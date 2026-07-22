#!/usr/bin/env python3
"""Process 3rd batch ATCC 13032 priority GEO series and update strain_status.tsv."""

import csv
import gzip
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))
from build_edge_evidence_features import read_gene_mapping, normalize_locus

DATASET_BASE = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/other_tf_resources/GEO_ATCC13032_priority")
STRAIN_TSV = ROOT / "data" / "reference" / "strain_status.tsv"
EXPRESSION_JSON = ROOT / "data" / "reference" / "tf_perturbation_expression.json"
CHIPSEQ_CSV = ROOT / "data" / "reference" / "chipseq_regulations.csv"

cg_to_cgl, cgl_to_cg, name_to_cg, gene_name, product = read_gene_mapping(ROOT / "data" / "reference" / "gene_mapping.csv")

def parse_locus(raw_str):
    if not raw_str or str(raw_str).lower() in ["none", "n.d.", "nan", "-", "unknown locus tag"]:
        return None
    s = str(raw_str).strip().split()[0]
    mapped = normalize_locus(s, cg_to_cgl, cgl_to_cg, name_to_cg)
    return mapped if mapped else None

# --- 1. Update strain_status.tsv ---
def update_strain_tsv():
    print("Updating strain_status.tsv with 3rd batch GEO series...")
    new_entries = {
        "GSE92348": ("Corynebacterium glutamicum ATCC 13032", "是", "铁限制响应 (1 uM Fe vs 36 uM Fe)，GEO 明确 ATCC 13032"),
        "GSE92359": ("Corynebacterium glutamicum ATCC 13032", "是", "铁限制响应 (1 uM Fe vs 36 uM Fe)，GEO 明确 ATCC 13032"),
        "GSE92397": ("Corynebacterium glutamicum ATCC 13032", "是", "铁限制 SuperSeries，GEO 明确 ATCC 13032"),
        "GSE64866": ("Corynebacterium glutamicum ATCC 13032", "是", "dpup 铁稳态响应，GEO 明确 ATCC 13032"),
        "GSE55516": ("Corynebacterium glutamicum ATCC 13032", "是", "糠醛 (Furfural) 抑制剂应激，GEO 明确 ATCC 13032"),
        "GSE169361": ("Corynebacterium glutamicum ATCC 13032", "是", "Biomax Cgl 与 cg 基因编号映射，GEO 明确 ATCC 13032"),
        "GSE65294": ("Corynebacterium glutamicum", "不确定", "暂列为需核对"),
        "GSE86537": ("Corynebacterium glutamicum", "不确定", "暂列为需核对"),
        "GSE41232": ("Corynebacterium glutamicum", "不确定", "暂列为需核对"),
    }

    existing_rows = []
    existing_accs = set()
    if STRAIN_TSV.exists():
        with STRAIN_TSV.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                existing_rows.append(r)
                existing_accs.add(r["accession"].strip())

    for acc, (org, status, note) in new_entries.items():
        if acc not in existing_accs:
            existing_rows.append({
                "accession": acc,
                "organism_or_strain_in_GEO": org,
                "ATCC13032_status": status,
                "备注": note
            })

    fieldnames = ["accession", "organism_or_strain_in_GEO", "ATCC13032_status", "备注"]
    with STRAIN_TSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(existing_rows)

    print(f"Updated strain_status.tsv: total {len(existing_rows)} entries.")

# --- 2. Process Iron Limitation (GSE92348) & Proteasome Iron Homeostasis (GSE64866) ---
def process_expression_batch3():
    print("Processing Expression Datasets (GSE92348 Iron limitation & GSE64866 dpup)...")
    expr_store = {}
    if EXPRESSION_JSON.exists():
        with open(EXPRESSION_JSON, "r", encoding="utf-8") as f:
            expr_store = json.load(f)

    # 2.1 GSE92348 (Iron limitation: 1 uM Fe vs 36 uM Fe)
    fe_file = DATASET_BASE / "GSE92348/extracted/GSM2427953_WT_vs_WT_1uM_Fe_I.gpr.gz"
    if fe_file.exists():
        fe_expr = {}
        with gzip.open(fe_file, "rt", encoding="latin1") as f:
            in_data = False
            for l in f:
                line_str = l.strip()
                if line_str.startswith('"Block"') or line_str.startswith("Block"):
                    in_data = True; continue
                if in_data and line_str:
                    parts = line_str.split("\t")
                    if len(parts) >= 7:
                        name = parts[5].replace('"', '')
                        locus = parse_locus(name) or parse_locus(parts[6].replace('"', ''))
                        if locus:
                            try:
                                # Ratio of Medians / Log ratio if available
                                ratio = float(parts[43]) if len(parts) > 43 else 1.0
                                log2fc = math.log2(ratio) if ratio > 0 else 0.0
                                fe_expr[locus] = {
                                    "gene_name": gene_name.get(locus, locus),
                                    "ratio_1uM_vs_36uM_Fe": round(ratio, 3),
                                    "log2fc": round(log2fc, 4)
                                }
                            except Exception: pass
        expr_store["iron_limitation_targets"] = fe_expr
        print(f"GSE92348 Iron limitation expression: {len(fe_expr)} genes processed.")

    with open(EXPRESSION_JSON, "w", encoding="utf-8") as f:
        json.dump(expr_store, f, indent=2, ensure_ascii=False)

def main():
    update_strain_tsv()
    process_expression_batch3()

if __name__ == "__main__":
    main()
