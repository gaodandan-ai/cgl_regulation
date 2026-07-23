#!/usr/bin/env python3
"""Copy strain_status.tsv to data/reference and update strain_group in chipseq_regulations.csv."""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data_pipeline" / "scripts"))

SRC_TSV = Path(r"C:/Users/Tsuki/Documents/Codex/2026-07-21/interconnected-iron-and-heme-regulatory-networks/outputs/other_tf_resources/strain_status.tsv")
DEST_TSV = ROOT / "data" / "reference" / "strain_status.tsv"
CHIPSEQ_CSV = ROOT / "data" / "reference" / "chipseq_regulations.csv"

def main():
    # 1. Copy strain_status.tsv
    if SRC_TSV.exists():
        DEST_TSV.parent.mkdir(parents=True, exist_ok=True)
        content = SRC_TSV.read_text(encoding="utf-8")
        DEST_TSV.write_text(content, encoding="utf-8")
        print(f"Copied strain_status.tsv to {DEST_TSV}")

    # 2. Build GEO Accession to Strain Group map
    acc_map = {}
    if DEST_TSV.exists():
        with DEST_TSV.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                acc = r.get("accession", "").strip()
                status = r.get("ATCC13032_status", "").strip()
                if status == "是":
                    group = "ATCC13032"
                elif status == "否":
                    group = "Strain_R"
                else:
                    group = "Unspecified"
                if acc:
                    acc_map[acc] = group

    # 3. Update chipseq_regulations.csv with strain_group
    rows = []
    if CHIPSEQ_CSV.exists():
        with CHIPSEQ_CSV.open("r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for r in reader:
                src = r.get("Source", "")
                strain_note = r.get("strain_note", "")

                # Determine strain group
                group = "ATCC13032" # Default for literature ATCC13032
                for acc, grp in acc_map.items():
                    if acc in src:
                        group = grp
                        break

                if "cross-strain" in strain_note.lower():
                    group = "Strain_ATCC14067"
                elif "Strain_R" in src or "GSE58632" in src or "GSE52040" in src or "GSE10232" in src or "GSE72452" in src:
                    group = "Strain_R"

                r["strain_group"] = group
                rows.append(r)

    fieldnames = [
        "TF_locusTag", "TF_altLocusTag", "TF_name", "TF_role",
        "TG_locusTag", "TG_altLocusTag", "TG_name", "Operon",
        "Binding_site", "Role", "Is_sigma_factor", "Evidence",
        "PMID", "Source", "evidence_score", "confidence_label", "strain_note", "strain_group"
    ]
    with CHIPSEQ_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Count breakdown
    group_counts = {}
    for r in rows:
        g = r["strain_group"]
        group_counts[g] = group_counts.get(g, 0) + 1

    print(f"Updated {CHIPSEQ_CSV} with strain_group annotations:")
    for g, cnt in sorted(group_counts.items()):
        print(f"  - {g}: {cnt} interactions")

if __name__ == "__main__":
    main()
