#!/usr/bin/env python3
"""
process_chipseq_peaks_integration.py
=====================================
Processes and ingests lab_chipseq_peaks.csv (36.8k peak records), lab_chipseq_regulations.csv,
and chipseq_regulations.csv into the SQLite databases (data/cgl_regulation.db and data/reference/cgl_regulation.db).

Computes spatial promoter/TSS proximity confidence tiers and creates fast spatial/gene indexes.

Usage:
    python data_pipeline/scripts/process_chipseq_peaks_integration.py
"""

import os
import sys
import sqlite3
import pandas as pd
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parents[2]
REF_DIR = ROOT_DIR / "data" / "reference"
DATA_DB_PATH = ROOT_DIR / "data" / "cgl_regulation.db"
REF_DB_PATH = REF_DIR / "cgl_regulation.db"

PEAKS_CSV = REF_DIR / "lab_chipseq_peaks.csv"
LAB_REG_CSV = REF_DIR / "lab_chipseq_regulations.csv"
PUB_REG_CSV = REF_DIR / "chipseq_regulations.csv"


def classify_spatial_confidence(row):
    region = str(row.get("genomic_region_chip", "") or "").lower()
    context = str(row.get("tu_context", "") or "").lower()
    rel_tss = row.get("rel_pos_to_TSS")

    try:
        rel_tss_val = float(rel_tss)
    except (ValueError, TypeError):
        rel_tss_val = None

    if region == "promoter" or (rel_tss_val is not None and -350 <= rel_tss_val <= 50):
        return "PROMOTER_DIRECT"
    elif "upstream" in context or "promoter" in context:
        return "INTERGENIC_PROMOTER"
    elif region == "gene_body" or "inside_transcribed" in context:
        return "GENE_BODY_INTERNAL"
    else:
        return "UPSTREAM_DISTAL"


def process_and_ingest(db_path: Path):
    print(f"--- Ingesting ChIP-seq peaks and regulations into {db_path} ---")
    if not db_path.exists():
        print(f"Warning: Database {db_path} does not exist, creating new...")
        db_path.parent.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("PRAGMA journal_mode=WAL;")
    cursor.execute("PRAGMA synchronous=NORMAL;")

    # 1. Process ChIP-seq Peaks
    if PEAKS_CSV.exists():
        print(f"Reading {PEAKS_CSV}...")
        df_peaks = pd.read_csv(PEAKS_CSV)

        # Standardize column names to lowercase snake_case
        df_peaks.rename(columns={
            "TF_id": "tf_id",
            "TF_name": "tf_name",
            "rel_pos_to_TSS": "rel_pos_to_tss",
            "TSS_like": "tss_like",
            "TU_id": "tu_id",
        }, inplace=True)

        # Add spatial_confidence column
        df_peaks["spatial_confidence"] = df_peaks.apply(classify_spatial_confidence, axis=1)

        # Drop existing table if any & re-create
        cursor.execute("DROP TABLE IF EXISTS chipseq_peaks;")
        cursor.execute("""
            CREATE TABLE chipseq_peaks (
                peak_id TEXT PRIMARY KEY,
                tf_id TEXT,
                tf_name TEXT,
                chip_prefix TEXT,
                chrom TEXT,
                peak_start INT,
                peak_end INT,
                peak_center INT,
                peak_score REAL,
                peak_signal REAL,
                neglog10q REAL,
                strength_tier TEXT,
                nearest_gene_locus TEXT,
                nearest_gene_name TEXT,
                tu_id TEXT,
                gene_list TEXT,
                strand TEXT,
                tss_like INT,
                promoter_start INT,
                promoter_end INT,
                rel_pos_to_tss INT,
                abs_rel_pos INT,
                overlap_bp INT,
                mapping_status TEXT,
                tu_context TEXT,
                genomic_region_chip TEXT,
                spatial_confidence TEXT
            );
        """)

        # Ensure correct column ordering
        cols = [
            "peak_id", "tf_id", "tf_name", "chip_prefix", "chrom",
            "peak_start", "peak_end", "peak_center", "peak_score", "peak_signal",
            "neglog10q", "strength_tier", "nearest_gene_locus", "nearest_gene_name",
            "tu_id", "gene_list", "strand", "tss_like", "promoter_start",
            "promoter_end", "rel_pos_to_tss", "abs_rel_pos", "overlap_bp",
            "mapping_status", "tu_context", "genomic_region_chip", "spatial_confidence"
        ]

        # Filter cols that actually exist in dataframe
        avail_cols = [c for c in cols if c in df_peaks.columns]
        df_to_save = df_peaks[avail_cols].drop_duplicates(subset=["peak_id"]).reset_index(drop=True)
        df_to_save.to_sql("chipseq_peaks", conn, if_exists="append", index=False)

        # Build Indexes
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipseq_peaks_coord ON chipseq_peaks(peak_start, peak_end);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipseq_peaks_gene ON chipseq_peaks(nearest_gene_locus);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipseq_peaks_tf ON chipseq_peaks(tf_name);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipseq_peaks_tfid ON chipseq_peaks(tf_id);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipseq_peaks_region ON chipseq_peaks(genomic_region_chip);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipseq_peaks_conf ON chipseq_peaks(spatial_confidence);")

        print(f"Successfully loaded {len(df_to_save)} peak records into 'chipseq_peaks'.")
    else:
        print(f"Warning: {PEAKS_CSV} not found.")

    # 2. Process ChIP-seq Regulations (Lab + Public Literature)
    reg_rows = []
    if LAB_REG_CSV.exists():
        print(f"Reading {LAB_REG_CSV}...")
        df_lab_reg = pd.read_csv(LAB_REG_CSV)
        reg_rows.append(df_lab_reg)

    if PUB_REG_CSV.exists():
        print(f"Reading {PUB_REG_CSV}...")
        df_pub_reg = pd.read_csv(PUB_REG_CSV)
        reg_rows.append(df_pub_reg)

    if reg_rows:
        df_all_reg = pd.concat(reg_rows, ignore_index=True)
        # Rename columns to standard
        df_all_reg.rename(columns={
            "TF_locusTag": "tf_locus_tag",
            "TF_altLocusTag": "tf_alt_locus_tag",
            "TF_name": "tf_name",
            "TF_role": "tf_role",
            "TG_locusTag": "tg_locus_tag",
            "TG_altLocusTag": "tg_alt_locus_tag",
            "TG_name": "tg_name",
            "Operon": "operon",
            "Binding_site": "binding_site",
            "Role": "role",
            "Is_sigma_factor": "is_sigma_factor",
            "Evidence": "evidence",
            "PMID": "pmid",
            "Source": "source",
            "evidence_score": "evidence_score",
            "confidence_label": "confidence_label",
            "strain_note": "strain_note",
            "strain_group": "strain_group"
        }, inplace=True)

        cursor.execute("DROP TABLE IF EXISTS chipseq_regulations;")
        cursor.execute("""
            CREATE TABLE chipseq_regulations (
                tf_locus_tag TEXT,
                tf_alt_locus_tag TEXT,
                tf_name TEXT,
                tf_role TEXT,
                tg_locus_tag TEXT,
                tg_alt_locus_tag TEXT,
                tg_name TEXT,
                operon TEXT,
                binding_site TEXT,
                role TEXT,
                is_sigma_factor TEXT,
                evidence TEXT,
                pmid TEXT,
                source TEXT,
                evidence_score REAL,
                confidence_label TEXT,
                strain_note TEXT,
                strain_group TEXT
            );
        """)

        reg_cols = [
            "tf_locus_tag", "tf_alt_locus_tag", "tf_name", "tf_role",
            "tg_locus_tag", "tg_alt_locus_tag", "tg_name", "operon",
            "binding_site", "role", "is_sigma_factor", "evidence",
            "pmid", "source", "evidence_score", "confidence_label",
            "strain_note", "strain_group"
        ]
        avail_reg_cols = [c for c in reg_cols if c in df_all_reg.columns]
        df_reg_to_save = df_all_reg[avail_reg_cols].drop_duplicates().reset_index(drop=True)
        df_reg_to_save.to_sql("chipseq_regulations", conn, if_exists="append", index=False)

        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipreg_tf ON chipseq_regulations(tf_locus_tag);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipreg_tg ON chipseq_regulations(tg_locus_tag);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipreg_source ON chipseq_regulations(source);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_chipreg_tfname ON chipseq_regulations(tf_name);")

        print(f"Successfully loaded {len(df_reg_to_save)} regulation records into 'chipseq_regulations'.")

    conn.commit()
    conn.close()


def main():
    if DATA_DB_PATH.exists():
        process_and_ingest(DATA_DB_PATH)
    if REF_DB_PATH.exists():
        process_and_ingest(REF_DB_PATH)
    print("Done processing ChIP-seq integration.")


if __name__ == "__main__":
    main()
