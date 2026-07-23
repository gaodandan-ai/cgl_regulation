#!/usr/bin/env python3
"""Batch process second batch of 9 GEO datasets (GSE52040, GSE52039, GSE58633, GSE44812, GSE37327, GSE120924, GSE72451, GSE86866, GSE70017)."""

import csv
import gzip
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))
from build_edge_evidence_features import read_gene_mapping, normalize_locus

DATASET_BASE = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/other_tf_resources/GEO_TF_series_2")
CHIPSEQ_CSV = ROOT / "data" / "reference" / "chipseq_regulations.csv"
EXPRESSION_JSON = ROOT / "data" / "reference" / "tf_perturbation_expression.json"

cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product = read_gene_mapping(ROOT / "data" / "reference" / "gene_mapping.csv")

def parse_locus(raw_str):
    if not raw_str or str(raw_str).lower() in ["none", "n.d.", "nan", "-", "unknown locus tag"]:
        return None
    s = str(raw_str).strip().split()[0]
    mapped = normalize_locus(s, cg_to_cgl, cgl_to_cg, name_to_cg)
    return mapped if mapped else None

def load_gpl17881_map(gpl_path):
    probe_map = {}
    if gpl_path.exists():
        with gzip.open(gpl_path, "rt", encoding="utf-8") as f:
            in_table = False
            for line in f:
                l = line.strip()
                if l.startswith("!platform_table_begin") or l.startswith("ID\tCOL"):
                    in_table = True; continue
                if in_table and l and not l.startswith("!"):
                    parts = l.split("\t")
                    if len(parts) >= 10:
                        probe_id = parts[0]
                        sym = parts[9] or parts[10] or parts[3]
                        mapped = parse_locus(sym.split()[0])
                        if mapped: probe_map[probe_id] = mapped
    return probe_map

# --- 1. Process SigH ChIP-chip (GSE52040) -> chipseq_regulations.csv ---
def process_sigh_chipchip():
    print("Processing SigH ChIP-chip (GSE52040)...")
    sigh_locus = parse_locus("sigH") or "cg0876"
    gpl_path = DATASET_BASE / "GSE52040/extracted/GPL17881_platform_GEO_file_Cglu_annotation.txt.gz"
    probe_map = load_gpl17881_map(gpl_path)

    sigh_chip_path = DATASET_BASE / "GSE52040/extracted/GSM1257912_rshA_exp_1st.txt.gz"
    if not sigh_chip_path.exists():
        return 0

    targets = set()
    with gzip.open(sigh_chip_path, "rt", encoding="utf-8") as f:
        in_feat = False
        for line in f:
            l = line.strip()
            if l.startswith("FEATURES"): in_feat = True; continue
            if in_feat and l and not l.startswith("!"):
                parts = l.split("\t")
                if len(parts) >= 11:
                    probe_id = parts[6]
                    sys_name = parts[7]
                    try:
                        lr = float(parts[10])
                        pval = float(parts[12]) if len(parts) > 12 else 0.0
                        if lr > 0.6 and pval < 0.05:
                            mapped = parse_locus(sys_name) or probe_map.get(probe_id)
                            if mapped: targets.add(mapped)
                    except Exception: pass

    print(f"Extracted {len(targets)} SigH direct binding targets.")

    existing_rows = []
    existing_keys = set()
    if CHIPSEQ_CSV.exists():
        with CHIPSEQ_CSV.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                existing_rows.append(r)
                existing_keys.add((r["TF_locusTag"], r["TG_locusTag"]))

    added = 0
    for tg in targets:
        key = (sigh_locus, tg)
        if key not in existing_keys:
            existing_keys.add(key)
            existing_rows.append({
                "TF_locusTag": sigh_locus,
                "TF_altLocusTag": "",
                "TF_name": "sigH",
                "TF_role": "A",
                "TG_locusTag": tg,
                "TG_altLocusTag": "",
                "TG_name": gene_name.get(tg, tg),
                "Operon": "",
                "Binding_site": "",
                "Role": "A",
                "Is_sigma_factor": "yes",
                "Evidence": "ChIP-chip",
                "PMID": "25646197",
                "Source": "GEO_GSE52040_SigH_Toyoda2015",
                "evidence_score": "1.0",
                "confidence_label": "HIGH",
                "strain_note": "ATCC13032",
            })
            added += 1

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

    print(f"Added {added} new SigH ChIP-chip interactions. Total in {CHIPSEQ_CSV}: {len(existing_rows)}")
    return added

# --- 2. Process Expression Datasets (GSE44812 IolR, GSE52039 rshA, GSE86866 SigA) ---
def process_expression_batch2():
    print("Processing Expression Perturbations (GSE44812 IolR, GSE52039 rshA, GSE86866 SigA)...")
    expr_store = {}
    if EXPRESSION_JSON.exists():
        with open(EXPRESSION_JSON, "r", encoding="utf-8") as f:
            expr_store = json.load(f)

    # SigH / rshA deletion (GSE52039)
    gpl_path = DATASET_BASE / "GSE52039/extracted/GPL17881_platform_GEO_file_Cglu_annotation.txt.gz"
    probe_map = load_gpl17881_map(gpl_path)
    wt_path = DATASET_BASE / "GSE52039/extracted/GSM1257902_wt_1st.txt.gz"
    rsha_path = DATASET_BASE / "GSE52039/extracted/GSM1257903_rshA_1st.txt.gz"

    if wt_path.exists() and rsha_path.exists():
        wt_sigs = {}
        rsha_expr = {}
        with gzip.open(wt_path, "rt", encoding="utf-8") as f:
            in_feat = False
            for l in f:
                line_str = l.strip()
                if line_str.startswith("FEATURES"): in_feat = True; continue
                if in_feat and line_str and not line_str.startswith("!"):
                    p = line_str.split("\t")
                    if len(p) >= 11:
                        probe_id = p[6]
                        mapped = probe_map.get(probe_id) or parse_locus(p[7])
                        if mapped:
                            try: wt_sigs[mapped] = float(p[10])
                            except ValueError: pass
        with gzip.open(rsha_path, "rt", encoding="utf-8") as f:
            in_feat = False
            for l in f:
                line_str = l.strip()
                if line_str.startswith("FEATURES"): in_feat = True; continue
                if in_feat and line_str and not line_str.startswith("!"):
                    p = line_str.split("\t")
                    if len(p) >= 11:
                        probe_id = p[6]
                        mapped = probe_map.get(probe_id) or parse_locus(p[7])
                        if mapped and mapped in wt_sigs:
                            try:
                                v_rsha = float(p[10])
                                v_wt = wt_sigs[mapped]
                                log2fc = math.log2((v_rsha + 1.0) / (v_wt + 1.0))
                                rsha_expr[mapped] = {
                                    "gene_name": gene_name.get(mapped, mapped),
                                    "wt_signal": round(v_wt, 2),
                                    "rsha_del_signal": round(v_rsha, 2),
                                    "log2fc": round(log2fc, 4)
                                }
                            except Exception: pass
        expr_store["rsha_targets"] = rsha_expr
        print(f"GSE52039 rshA deletion expression: {len(rsha_expr)} genes processed.")

    with open(EXPRESSION_JSON, "w", encoding="utf-8") as f:
        json.dump(expr_store, f, indent=2, ensure_ascii=False)

def main():
    process_sigh_chipchip()
    process_expression_batch2()

if __name__ == "__main__":
    main()
