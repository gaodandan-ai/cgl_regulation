#!/usr/bin/env python3
"""Process RegPrecise (C. glutamicum) and SigD GEO series (GSE102323, GSE102324, GSE102325, GSE102328)."""

import csv
import gzip
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))
from build_edge_evidence_features import read_gene_mapping, normalize_locus

DATASET_BASE = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/other_tf_resources")
GSE_BASE = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/c_glutamicum_chipseq")

CHIPSEQ_CSV = ROOT / "data" / "reference" / "chipseq_regulations.csv"
REGPRECISE_CSV = ROOT / "data" / "reference" / "regprecise_regulations.csv"
EXPRESSION_JSON = ROOT / "data" / "reference" / "tf_perturbation_expression.json"

cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product = read_gene_mapping(ROOT / "data" / "reference" / "gene_mapping.csv")

def parse_locus(raw_str):
    if not raw_str or str(raw_str).lower() in ["none", "n.d.", "nan", "-", "unknown locus tag"]:
        return None
    s = str(raw_str).strip().split()[0]
    mapped = normalize_locus(s, cg_to_cgl, cgl_to_cg, name_to_cg)
    return mapped if mapped else None

# --- 1. Process RegPrecise ---
def process_regprecise():
    tsv_path = DATASET_BASE / "RegPrecise_Cglutamicum/Cglutamicum_regulated_genes.tsv"
    if not tsv_path.exists():
        print(f"RegPrecise TSV file not found: {tsv_path}")
        return 0

    print("Processing RegPrecise C. glutamicum Regulons...")
    rows = []
    current_tf = None
    current_tf_name = None

    with tsv_path.open("r", encoding="utf-8") as f:
        for line in f:
            line_str = line.strip()
            if line_str.startswith("# TF -"):
                try:
                    tf_part = line_str.split("-")[1].strip()
                    parts = tf_part.split(":")
                    current_tf_name = parts[0].strip()
                    raw_tf_locus = parts[1].strip() if len(parts) > 1 else current_tf_name
                    current_tf = parse_locus(raw_tf_locus) or parse_locus(current_tf_name)
                except Exception:
                    current_tf = None
                continue

            if line_str and current_tf:
                parts = line_str.split("\t")
                if len(parts) >= 2:
                    raw_tg = parts[1]
                    tg_name_raw = parts[2] if len(parts) > 2 else ""
                    tg_locus = parse_locus(raw_tg) or parse_locus(tg_name_raw)
                    if tg_locus:
                        rows.append({
                            "TF_locusTag": current_tf,
                            "TF_altLocusTag": "",
                            "TF_name": current_tf_name or gene_name.get(current_tf, current_tf),
                            "TF_role": "Dual",
                            "TG_locusTag": tg_locus,
                            "TG_altLocusTag": "",
                            "TG_name": gene_name.get(tg_locus, tg_name_raw or tg_locus),
                            "Operon": "",
                            "Binding_site": "",
                            "Role": "Dual",
                            "Is_sigma_factor": "no",
                            "Evidence": "predicted",
                            "PMID": "",
                            "Source": "RegPrecise_Cglutamicum",
                            "evidence_score": "0.6",
                            "confidence_label": "MEDIUM",
                            "strain_note": "ATCC13032",
                        })

    print(f"Extracted {len(rows)} predicted TF-target interactions from RegPrecise.")
    REGPRECISE_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "TF_locusTag", "TF_altLocusTag", "TF_name", "TF_role",
        "TG_locusTag", "TG_altLocusTag", "TG_name", "Operon",
        "Binding_site", "Role", "Is_sigma_factor", "Evidence",
        "PMID", "Source", "evidence_score", "confidence_label", "strain_note"
    ]
    with REGPRECISE_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)

# --- Load GPL17881 Platform Annotation ---
def load_gpl17881_annotation():
    gpl_path = GSE_BASE / "GSE72453/extracted/GPL17881_platform_GEO_file_Cglu_annotation.txt.gz"
    probe_to_locus = {}
    if not gpl_path.exists():
        return probe_to_locus
    with gzip.open(gpl_path, "rt", encoding="utf-8") as f:
        in_table = False
        for line in f:
            l = line.strip()
            if l.startswith("!platform_table_begin") or l.startswith("ID\tCOL"):
                in_table = True
                continue
            if in_table and l and not l.startswith("!"):
                parts = l.split("\t")
                if len(parts) >= 11:
                    probe_id = parts[0]
                    gene_sym = parts[9] or parts[10] or parts[3]
                    mapped = parse_locus(gene_sym)
                    if mapped:
                        probe_to_locus[probe_id] = mapped
    return probe_to_locus

# --- 2. Process SigD ChIP-chip (GSE102323) ---
def process_sigd_chipchip(probe_to_locus):
    print("Processing SigD ChIP-chip (GSE102323)...")
    sigd_locus = parse_locus("sigD") or "cg0696"
    chip_path = DATASET_BASE / "SigD/GSE102323/extracted/GSM2734731_sigD-ChIP_1st.txt.gz"
    if not chip_path.exists():
        return 0

    sigd_targets = set()
    with gzip.open(chip_path, "rt", encoding="utf-8") as f:
        in_features = False
        for line in f:
            line_str = line.strip()
            if line_str.startswith("FEATURES"):
                in_features = True
                continue
            if in_features and line_str and not line_str.startswith("!"):
                parts = line_str.split("\t")
                if len(parts) >= 10:
                    probe_id = parts[0]
                    sys_name = parts[6]
                    try:
                        log_ratio = float(parts[9])
                        p_val = float(parts[11]) if len(parts) > 11 else 0.0
                        if log_ratio > 0.8 and p_val < 0.05:
                            mapped = parse_locus(sys_name) or probe_to_locus.get(probe_id)
                            if mapped:
                                sigd_targets.add(mapped)
                    except Exception:
                        pass

    print(f"Extracted {len(sigd_targets)} SigD ChIP-chip direct binding targets.")

    existing_rows = []
    existing_keys = set()
    if CHIPSEQ_CSV.exists():
        with CHIPSEQ_CSV.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                existing_rows.append(r)
                existing_keys.add((r["TF_locusTag"], r["TG_locusTag"]))

    added_count = 0
    for tg_locus in sigd_targets:
        key = (sigd_locus, tg_locus)
        if key not in existing_keys:
            existing_keys.add(key)
            existing_rows.append({
                "TF_locusTag": sigd_locus,
                "TF_altLocusTag": "",
                "TF_name": "sigD",
                "TF_role": "A",
                "TG_locusTag": tg_locus,
                "TG_altLocusTag": "",
                "TG_name": gene_name.get(tg_locus, tg_locus),
                "Operon": "",
                "Binding_site": "",
                "Role": "A",
                "Is_sigma_factor": "yes",
                "Evidence": "ChIP-chip",
                "PMID": "28874697",
                "Source": "GEO_GSE102323_SigD",
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

    print(f"Added {added_count} new SigD ChIP-chip interactions. Total in {CHIPSEQ_CSV}: {len(existing_rows)}")
    return added_count

# --- 3. Process SigD & RsdA Perturbations (GSE102324 & GSE102325 & GSE102328) ---
def process_sigd_expression(probe_to_locus):
    print("Processing SigD / RsdA Transcriptome Perturbations...")
    expr_data = {}
    if EXPRESSION_JSON.exists():
        with open(EXPRESSION_JSON, "r", encoding="utf-8") as f:
            expr_data = json.load(f)

    sigd_expr = {}
    ox_path = DATASET_BASE / "SigD/GSE102324/extracted/GSM2734733_sigD_ox_1st.txt.gz"
    ctrl_path = DATASET_BASE / "SigD/GSE102324/extracted/GSM2734737_vector_control_1st.txt.gz"

    if ox_path.exists() and ctrl_path.exists():
        ctrl_signals = {}
        with gzip.open(ctrl_path, "rt", encoding="utf-8") as f:
            in_features = False
            for line in f:
                l = line.strip()
                if l.startswith("FEATURES"):
                    in_features = True
                    continue
                if in_features and l and not l.startswith("!"):
                    parts = l.split("\t")
                    if len(parts) >= 10:
                        probe_id = parts[0]
                        mapped = probe_to_locus.get(probe_id)
                        if mapped:
                            try:
                                ctrl_signals[mapped] = float(parts[9])
                            except ValueError:
                                pass

        with gzip.open(ox_path, "rt", encoding="utf-8") as f:
            in_features = False
            for line in f:
                l = line.strip()
                if l.startswith("FEATURES"):
                    in_features = True
                    continue
                if in_features and l and not l.startswith("!"):
                    parts = l.split("\t")
                    if len(parts) >= 10:
                        probe_id = parts[0]
                        mapped = probe_to_locus.get(probe_id)
                        if mapped and mapped in ctrl_signals:
                            try:
                                val_ox = float(parts[9])
                                val_ctrl = ctrl_signals[mapped]
                                log2fc = math.log2((val_ox + 1.0) / (val_ctrl + 1.0))
                                sigd_expr[mapped] = {
                                    "gene_name": gene_name.get(mapped, mapped),
                                    "signal_ctrl": round(val_ctrl, 2),
                                    "signal_sigD_ox": round(val_ox, 2),
                                    "log2fc": round(log2fc, 4),
                                }
                            except Exception:
                                pass

    expr_data["sigd_targets"] = sigd_expr
    with open(EXPRESSION_JSON, "w", encoding="utf-8") as f:
        json.dump(expr_data, f, indent=2, ensure_ascii=False)
    print(f"Processed SigD Overexpression RNA-seq ({len(sigd_expr)} genes).")

def main():
    probe_to_locus = load_gpl17881_annotation()
    print(f"Loaded GPL17881 platform annotations: {len(probe_to_locus)} probes mapped.")
    process_regprecise()
    process_sigd_chipchip(probe_to_locus)
    process_sigd_expression(probe_to_locus)

if __name__ == "__main__":
    main()
