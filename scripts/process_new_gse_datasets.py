#!/usr/bin/env python3
"""Process and integrate GSE156688 (TSS), GSE184402 (IpsA RNA-seq), GSE72453 (SigC RNA-seq), and GSE72452 (SigC ChIP-chip)."""

import csv
import gzip
import json
import math
import sys
import openpyxl
import xlrd
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))
from build_edge_evidence_features import read_gene_mapping, normalize_locus

DATASET_BASE = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/c_glutamicum_chipseq")
CHIPSEQ_CSV = ROOT / "data" / "reference" / "chipseq_regulations.csv"
TSS_JSON = ROOT / "data" / "reference" / "tss_promoter_annotations.json"
EXPRESSION_JSON = ROOT / "data" / "reference" / "tf_perturbation_expression.json"

cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product = read_gene_mapping(ROOT / "data" / "reference" / "gene_mapping.csv")

def parse_target_locus(raw_str):
    if not raw_str or str(raw_str).lower() in ["none", "n.d.", "nan", "-", "unknown locus tag"]:
        return None
    s = str(raw_str).strip().split()[0]
    mapped = normalize_locus(s, cg_to_cgl, cgl_to_cg, name_to_cg)
    return mapped if mapped else None

# --- 1. Process GSE156688 (TSS / 5' UTR Annotations) ---
def process_tss_dataset():
    xls_path = DATASET_BASE / "GSE156688/extracted/GSM4742245_TSS_detection.xls"
    if not xls_path.exists():
        print(f"TSS file not found: {xls_path}")
        return {}

    print("Processing GSE156688 (TSS detection)...")
    wb = xlrd.open_workbook(xls_path)
    s = wb.sheet_by_index(0)

    tss_map = {}
    for i in range(1, s.nrows):
        row = s.row_values(i)
        pos = row[0]
        strand = str(row[3])
        no_reads = row[4]
        tss_type = str(row[8])

        # Primary / Next downstream gene
        target_locus = parse_target_locus(row[11]) or parse_target_locus(row[24]) or parse_target_locus(row[17])
        if not target_locus:
            continue

        seq_70up = str(row[38]) if len(row) > 38 else ""
        leader_len = row[33] if len(row) > 33 and isinstance(row[33], (int, float)) else None

        if target_locus not in tss_map or (tss_type.lower() == "primary"):
            tss_map[target_locus] = {
                "locus": target_locus,
                "gene_name": gene_name.get(target_locus, target_locus),
                "tss_position": pos,
                "strand": strand,
                "reads": no_reads,
                "tss_type": tss_type,
                "leader_length": leader_len,
                "promoter_70bp_upstream": seq_70up,
            }

    print(f"Extracted TSS annotations for {len(tss_map)} C. glutamicum genes.")
    TSS_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(TSS_JSON, "w", encoding="utf-8") as f:
        json.dump(tss_map, f, indent=2, ensure_ascii=False)
    return tss_map

# --- 2. Process GSE184402 (IpsA RNA-seq FPKM) & GSE72453 (SigC RNA-seq) ---
def process_expression_perturbations():
    print("Processing RNA-seq perturbations (GSE184402 IpsA & GSE72453 SigC)...")
    results = {"ipsa_targets": {}, "sigc_targets": {}}

    # GSE184402 (IpsA)
    ipsa_path = DATASET_BASE / "GSE184402/GSE184402_all.fpkm.xls.gz"
    if ipsa_path.exists():
        with gzip.open(ipsa_path, "rt", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gname = row.get("GeneName") or row.get("GeneId")
                locus = parse_target_locus(gname) or parse_target_locus(row.get("GeneId"))
                if not locus:
                    continue

                # FPKM values
                try:
                    fpkm_cgs11 = (float(row.get("CGS11-1-2_FPKM", 0)) + float(row.get("CGS11-2-1_FPKM", 0)) + float(row.get("CGS11-3-1_FPKM", 0))) / 3.0
                    fpkm_cgs15 = (float(row.get("CGS15-1-2_FPKM", 0)) + float(row.get("CGS15-2-1_FPKM", 0)) + float(row.get("CGS15-3-2_FPKM", 0))) / 3.0
                    fpkm_cgs7  = (float(row.get("CGS7-1-2_FPKM", 0))  + float(row.get("CGS7-2-1_FPKM", 0))  + float(row.get("CGS7-3-1_FPKM", 0)))  / 3.0

                    log2fc = math.log2((fpkm_cgs15 + 0.1) / (fpkm_cgs11 + 0.1))
                    results["ipsa_targets"][locus] = {
                        "gene_name": gene_name.get(locus, locus),
                        "fpkm_wt": round(fpkm_cgs11, 2),
                        "fpkm_ipsa_mut": round(fpkm_cgs15, 2),
                        "log2fc": round(log2fc, 4),
                    }
                except Exception:
                    pass

    # GSE72453 (SigC RNA-seq)
    sigc_del_path = DATASET_BASE / "GSE72453/extracted/GSM1862830_sigCdel1_expo_1st.txt.gz"
    sigc_wt_path  = DATASET_BASE / "GSE72453/extracted/GSM1862829_wt_expo_1st.txt.gz"

    if sigc_del_path.exists() and sigc_wt_path.exists():
        wt_signals = {}
        with gzip.open(sigc_wt_path, "rt", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[0].startswith("cg"):
                    try:
                        wt_signals[parts[0]] = float(parts[1])
                    except ValueError:
                        pass

        with gzip.open(sigc_del_path, "rt", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[0].startswith("cg"):
                    locus = parse_target_locus(parts[0])
                    if locus and locus in wt_signals:
                        try:
                            val_del = float(parts[1])
                            val_wt = wt_signals[locus]
                            log2fc = math.log2((val_del + 1.0) / (val_wt + 1.0))
                            results["sigc_targets"][locus] = {
                                "gene_name": gene_name.get(locus, locus),
                                "signal_wt": round(val_wt, 2),
                                "signal_sigCdel": round(val_del, 2),
                                "log2fc": round(log2fc, 4),
                            }
                        except Exception:
                            pass

    EXPRESSION_JSON.parent.mkdir(parents=True, exist_ok=True)
    with open(EXPRESSION_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"Processed IpsA FPKM ({len(results['ipsa_targets'])} genes) & SigC RNA-seq ({len(results['sigc_targets'])} genes).")
    return results

# --- 3. Process GSE72452 (SigC ChIP-chip -> chipseq_regulations.csv) ---
def process_sigc_chipchip():
    print("Processing GSE72452 (SigC ChIP-chip)...")
    sigc_chip_path = DATASET_BASE / "GSE72452/extracted/GSM1862839_sigC_expo_1st.txt.gz"
    if not sigc_chip_path.exists():
        return 0

    sigc_locus = parse_target_locus("sigC") or "cg0309"
    sigc_targets = set()

    with gzip.open(sigc_chip_path, "rt", encoding="utf-8") as f:
        in_features = False
        for line in f:
            line_str = line.strip()
            if line_str.startswith("FEATURES"):
                in_features = True
                continue
            if in_features and line_str and not line_str.startswith("!"):
                parts = line_str.split("\t")
                if len(parts) >= 10:
                    probe_name = parts[5]
                    sys_name = parts[6]
                    try:
                        log_ratio = float(parts[9])
                        p_val = float(parts[11]) if len(parts) > 11 else 0.0
                        if log_ratio > 1.0 and p_val < 0.05:  # High-confidence ChIP binding peak
                            mapped = parse_target_locus(sys_name) or parse_target_locus(probe_name)
                            if mapped:
                                sigc_targets.add(mapped)
                    except Exception:
                        pass

    print(f"Extracted {len(sigc_targets)} SigC ChIP-chip direct binding targets.")

    # Read and update chipseq_regulations.csv
    existing_rows = []
    existing_keys = set()
    if CHIPSEQ_CSV.exists():
        with CHIPSEQ_CSV.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                existing_rows.append(r)
                existing_keys.add((r["TF_locusTag"], r["TG_locusTag"]))

    added_count = 0
    for tg_locus in sigc_targets:
        key = (sigc_locus, tg_locus)
        if key not in existing_keys:
            existing_keys.add(key)
            existing_rows.append({
                "TF_locusTag": sigc_locus,
                "TF_altLocusTag": "",
                "TF_name": "sigC",
                "TF_role": "A",
                "TG_locusTag": tg_locus,
                "TG_altLocusTag": "",
                "TG_name": gene_name.get(tg_locus, tg_locus),
                "Operon": "",
                "Binding_site": "",
                "Role": "A",
                "Is_sigma_factor": "yes",
                "Evidence": "ChIP-chip",
                "PMID": "26500293",
                "Source": "GEO_GSE72452_SigC",
                "evidence_score": "1.0",
                "confidence_label": "HIGH",
                "strain_note": "ATCC13032",
            })
            added_count += 1

    fieldnames = [
        "TF_locusTag", "TF_altLocusTag", "TF_name", "TF_role",
        "TG_locusTag", "TG_altLocusTag", "TG_name", "Operon",
        "Binding_site", "Role", "Is_sigma_factor", "Evidence",
        "PMID", "Source", "evidence_score", "confidence_label", "strain_note"
    ]
    with CHIPSEQ_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(existing_rows)

    print(f"Added {added_count} new SigC ChIP-chip interactions. Total in {CHIPSEQ_CSV}: {len(existing_rows)}")
    return added_count

def main():
    process_tss_dataset()
    process_expression_perturbations()
    process_sigc_chipchip()

if __name__ == "__main__":
    main()
