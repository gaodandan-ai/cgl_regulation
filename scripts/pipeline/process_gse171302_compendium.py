#!/usr/bin/env python3
"""Process GSE171302 Compendium (1214 samples), update strain_status.tsv and build quality-filtered index."""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))

STRAIN_TSV = ROOT / "data" / "reference" / "strain_status.tsv"
INDEX_CSV = ROOT / "data" / "reference" / "gse171302_metadata_index.csv"
GSE171302_DIR = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/other_tf_resources/GEO_ATCC13032_priority/GSE171302/extracted")

def update_strain_status():
    print("Registering GSE171302 in strain_status.tsv...")
    existing = []
    found = False
    if STRAIN_TSV.exists():
        with STRAIN_TSV.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                if r["accession"].strip() == "GSE171302":
                    r["organism_or_strain_in_GEO"] = "Corynebacterium glutamicum ATCC 13032"
                    r["ATCC13032_status"] = "是"
                    r["备注"] = "GEO 官方 ATCC 13032 表达谱 Compendium 全景库，关联 1,214 个样本；需按子系列和 QC 独立过滤使用"
                    found = True
                existing.append(r)

    if not found:
        existing.append({
            "accession": "GSE171302",
            "organism_or_strain_in_GEO": "Corynebacterium glutamicum ATCC 13032",
            "ATCC13032_status": "是",
            "备注": "GEO 官方 ATCC 13032 表达谱 Compendium 全景库，关联 1,214 个样本；需按子系列和 QC 独立过滤使用"
        })

    fieldnames = ["accession", "organism_or_strain_in_GEO", "ATCC13032_status", "备注"]
    with STRAIN_TSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(existing)

    print(f"Updated strain_status.tsv with GSE171302 entry ({len(existing)} total entries).")

def build_gse171302_index():
    print("Building GSE171302 sample quality and coverage index...")
    sample_files = list(GSE171302_DIR.glob("*.gz")) if GSE171302_DIR.exists() else []
    print(f"Found {len(sample_files)} sample files in GSE171302 directory.")

    records = []
    for p in sample_files:
        name = p.name
        gsm = name.split("_")[0]
        sz = p.stat().st_size
        qc_pass = "PASS" if sz > 500000 else "LOW_SIZE_WARN"
        records.append({
            "gsm_accession": gsm,
            "filename": name,
            "size_bytes": sz,
            "qc_status": qc_pass,
            "strain": "ATCC13032",
            "notes": "Compendium sub-sample"
        })

    INDEX_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["gsm_accession", "filename", "size_bytes", "qc_status", "strain", "notes"]
    with INDEX_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    print(f"Saved GSE171302 sample index to {INDEX_CSV} ({len(records)} samples processed).")

def main():
    update_strain_status()
    build_gse171302_index()

if __name__ == "__main__":
    main()
