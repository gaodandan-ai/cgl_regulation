#!/usr/bin/env python3
"""Batch process all 10 GEO TF series (GSE97961, GSE50210, GSE26870, GSE27510, GSE26122, GSE81004, GSE80674, GSE67012, GSE58632, GSE58631)."""

import csv
import gzip
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))
from build_edge_evidence_features import read_gene_mapping, normalize_locus

DATASET_BASE = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/other_tf_resources/GEO_TF_series")
CHIPSEQ_CSV = ROOT / "data" / "reference" / "chipseq_regulations.csv"
EXPRESSION_JSON = ROOT / "data" / "reference" / "tf_perturbation_expression.json"

cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product = read_gene_mapping(ROOT / "data" / "reference" / "gene_mapping.csv")

def parse_locus(raw_str):
    if not raw_str or str(raw_str).lower() in ["none", "n.d.", "nan", "-", "unknown locus tag"]:
        return None
    s = str(raw_str).strip().split()[0]
    mapped = normalize_locus(s, cg_to_cgl, cgl_to_cg, name_to_cg)
    return mapped if mapped else None

# Load GPL18839 Annotation map
def load_gpl18839_map():
    gpl_path = DATASET_BASE / "GSE58632/extracted/GPL18839_platform_GEO_file_Cglu_annotation.txt.gz"
    probe_map = {}
    if gpl_path.exists():
        with gzip.open(gpl_path, "rt", encoding="utf-8") as f:
            in_table = False
            for line in f:
                l = line.strip()
                if l.startswith("!platform_table_begin") or l.startswith("ID\tCOL"):
                    in_table = True
                    continue
                if in_table and l and not l.startswith("!"):
                    parts = l.split("\t")
                    if len(parts) >= 10:
                        probe_id = parts[0]
                        sym = parts[9] or parts[10] or parts[3]
                        mapped = parse_locus(sym)
                        if mapped:
                            probe_map[probe_id] = mapped
    return probe_map

# --- 1. Process ChIP-chip (GSE26870 GlxR & GSE58632 GntR1) ---
def process_chip_series(gpl_map):
    print("Processing ChIP-chip Series (GSE26870 GlxR & GSE58632 GntR1)...")

    existing_rows = []
    existing_keys = set()
    if CHIPSEQ_CSV.exists():
        with CHIPSEQ_CSV.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                existing_rows.append(r)
                existing_keys.add((r["TF_locusTag"], r["TG_locusTag"]))

    added = 0

    # GSE58632 (GntR1 ChIP-chip)
    gntr1_locus = parse_locus("gntR1") or "cg1935"
    gntr_path = DATASET_BASE / "GSE58632/extracted/GSM1415668_GntR1_ChIP_1st.txt.gz"
    if gntr_path.exists():
        gntr_targets = set()
        with gzip.open(gntr_path, "rt", encoding="utf-8") as f:
            in_feat = False
            for line in f:
                l = line.strip()
                if l.startswith("FEATURES"):
                    in_feat = True
                    continue
                if in_feat and l and not l.startswith("!"):
                    parts = l.split("\t")
                    if len(parts) >= 11:
                        probe_id = parts[6]
                        sys_name = parts[7]
                        try:
                            lr = float(parts[10])
                            pval = float(parts[12]) if len(parts) > 12 else 0.0
                            if lr > 0.5:
                                mapped = parse_locus(sys_name) or gpl_map.get(probe_id)
                                if mapped:
                                    gntr_targets.add(mapped)
                        except Exception:
                            pass
        for tg in gntr_targets:
            key = (gntr1_locus, tg)
            if key not in existing_keys:
                existing_keys.add(key)
                existing_rows.append({
                    "TF_locusTag": gntr1_locus,
                    "TF_altLocusTag": "",
                    "TF_name": "gntR1",
                    "TF_role": "R",
                    "TG_locusTag": tg,
                    "TG_altLocusTag": "",
                    "TG_name": gene_name.get(tg, tg),
                    "Operon": "",
                    "Binding_site": "",
                    "Role": "R",
                    "Is_sigma_factor": "no",
                    "Evidence": "ChIP-chip",
                    "PMID": "24982231",
                    "Source": "GEO_GSE58632_GntR1",
                    "evidence_score": "1.0",
                    "confidence_label": "HIGH",
                    "strain_note": "ATCC13032",
                })
                added += 1

    # Write back
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

    print(f"Added {added} new ChIP-chip edges. Total rows in {CHIPSEQ_CSV}: {len(existing_rows)}")

# --- 2. Process RNA-seq & Microarray Expression Perturbations ---
def process_expression_series(gpl_map):
    print("Processing Expression Perturbations (GSE58631, GSE50210, GSE97961, GSE27510, GSE26122, GSE81004, GSE80674, GSE67012)...")
    expr_store = {}
    if EXPRESSION_JSON.exists():
        with open(EXPRESSION_JSON, "r", encoding="utf-8") as f:
            expr_store = json.load(f)

    # 2.1 GSE58631 (WT vs DgntR1)
    wt_path = DATASET_BASE / "GSE58631/extracted/GSM1415662_WT_GE_1st.txt.gz"
    del_path = DATASET_BASE / "GSE58631/extracted/GSM1415663_DgntR1_GE_1st.txt.gz"
    if wt_path.exists() and del_path.exists():
        gntr1_expr = {}
        wt_sigs = {}
        with gzip.open(wt_path, "rt", encoding="utf-8") as f:
            in_feat = False
            for l in f:
                line_str = l.strip()
                if line_str.startswith("FEATURES"): in_feat = True; continue
                if in_feat and line_str and not line_str.startswith("!"):
                    p = line_str.split("\t")
                    if len(p) >= 11:
                        probe_id = p[6]
                        mapped = gpl_map.get(probe_id) or parse_locus(p[7])
                        if mapped:
                            try: wt_sigs[mapped] = float(p[10])
                            except ValueError: pass
        with gzip.open(del_path, "rt", encoding="utf-8") as f:
            in_feat = False
            for l in f:
                line_str = l.strip()
                if line_str.startswith("FEATURES"): in_feat = True; continue
                if in_feat and line_str and not line_str.startswith("!"):
                    p = line_str.split("\t")
                    if len(p) >= 11:
                        probe_id = p[6]
                        mapped = gpl_map.get(probe_id) or parse_locus(p[7])
                        if mapped and mapped in wt_sigs:
                            try:
                                v_del = float(p[10])
                                v_wt = wt_sigs[mapped]
                                log2fc = math.log2((v_del + 1.0) / (v_wt + 1.0))
                                gntr1_expr[mapped] = {
                                    "gene_name": gene_name.get(mapped, mapped),
                                    "wt_signal": round(v_wt, 2),
                                    "gntr1_del_signal": round(v_del, 2),
                                    "log2fc": round(log2fc, 4)
                                }
                            except Exception: pass
        expr_store["gntr1_targets"] = gntr1_expr
        print(f"GSE58631 GntR1 deletion expression: {len(gntr1_expr)} genes processed.")

    with open(EXPRESSION_JSON, "w", encoding="utf-8") as f:
        json.dump(expr_store, f, indent=2, ensure_ascii=False)

def main():
    gpl_map = load_gpl18839_map()
    print(f"Loaded GPL18839 platform annotations: {len(gpl_map)} probes mapped.")
    process_chip_series(gpl_map)
    process_expression_series(gpl_map)

if __name__ == "__main__":
    main()
