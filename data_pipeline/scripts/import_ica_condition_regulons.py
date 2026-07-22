#!/usr/bin/env python3
"""
import_ica_condition_regulons.py
================================
1. Builds `imodulon_condition_activities` table in cgl_regulation.db containing
   ICA-inferred iModulon activity profiles across transcriptomic conditions
   (Glucose Exponential, Acetate Growth, Phenol Stress, Iron Limitation,
    Cold Shock, Heat Shock, Oxygen Deprivation, Acid Stress).

2. Builds `imodulon_regulon_overlaps` table mapping ICA co-expression modules to
   known TF regulons with F1-scores, Precision, Recall, and overlapping gene lists.
"""

import os
import sys
import json
import sqlite3
import numpy as np
import pandas as pd

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF_DIR = os.path.join(ROOT_DIR, "data", "reference")
DB_PATH = os.path.join(REF_DIR, "cgl_regulation.db")
IMOD_META_PATH = os.path.join(REF_DIR, "imodulon", "imodulon_metadata.json")
IMOD_WEIGHTS_PATH = os.path.join(REF_DIR, "imodulon", "imodulon_gene_weights.json")
DE_MATRIX_PATH = os.path.join(
    REF_DIR, "expression_compendium", "raw", "Filtered differential expression results.xlsx"
)

CONDITIONS_LIST = [
    "Glucose Exponential",
    "Glucose Stationary",
    "Acetate Growth",
    "Phenol Stress",
    "Iron Limitation",
    "Cold Shock (15C)",
    "Heat Shock (40C)",
    "Oxygen Deprivation",
    "Acid Stress (pH 5.5)"
]

def import_ica_condition_regulons():
    print(f"Building ICA Condition-Specific Regulon Engine in: {DB_PATH}")
    if not os.path.exists(IMOD_META_PATH):
        print(f"ERROR: {IMOD_META_PATH} not found.")
        sys.exit(1)

    with open(IMOD_META_PATH, "r", encoding="utf-8") as f:
        imod_meta = json.load(f)
    if not os.path.exists(IMOD_WEIGHTS_PATH) or not os.path.exists(DE_MATRIX_PATH):
        raise FileNotFoundError("iModulon weights or differential-expression workbook is missing")
    with open(IMOD_WEIGHTS_PATH, "r", encoding="utf-8") as f:
        weight_data = json.load(f)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. Create imodulon_condition_activities table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS imodulon_condition_activities (
            activity_id INTEGER PRIMARY KEY AUTOINCREMENT,
            imodulon_id TEXT,
            imodulon_name TEXT,
            linked_regulator TEXT,
            category TEXT,
            condition_name TEXT,
            activity_score REAL,
            is_significant BOOLEAN,
            sample_id TEXT,
            raw_activity REAL,
            activity_method TEXT,
            data_source TEXT,
            is_measured BOOLEAN,
            gene_overlap_count INTEGER
        )
    """)

    # Project the real filtered differential-expression matrix onto the
    # thresholded iModulon gene-weight matrix using least squares.  These are
    # reproducible inferred activities, not hand-written or measured values.
    act_rows = []
    overlap_rows = []
    de_frame = pd.read_excel(DE_MATRIX_PATH, header=6)
    condition_col, sample_col = de_frame.columns[:2]
    expression = de_frame.iloc[:, 2:].apply(pd.to_numeric, errors="coerce")
    keep = de_frame[condition_col].notna() & (expression.notna().sum(axis=1) > 0)
    de_frame = de_frame.loc[keep].reset_index(drop=True)
    expression = expression.loc[keep].reset_index(drop=True)

    imod_ids = [item.get("id") for item in imod_meta if isinstance(item, dict) and item.get("id") in weight_data]
    genes = sorted(set(expression.columns) & {
        gene for imod_id in imod_ids for gene in weight_data[imod_id].get("genes", {})
    })
    weight_matrix = np.zeros((len(genes), len(imod_ids)), dtype=float)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    for imod_index, imod_id in enumerate(imod_ids):
        for gene, weight in weight_data[imod_id].get("genes", {}).items():
            if gene in gene_index:
                weight_matrix[gene_index[gene], imod_index] = float(weight)
    expression_matrix = expression[genes].fillna(0.0).to_numpy(dtype=float).T
    raw_activities = np.linalg.pinv(weight_matrix, rcond=1e-6) @ expression_matrix
    means = raw_activities.mean(axis=1, keepdims=True)
    stds = raw_activities.std(axis=1, keepdims=True)
    z_activities = np.divide(raw_activities - means, stds, out=np.zeros_like(raw_activities), where=stds > 0)
    meta_by_id = {item.get("id"): item for item in imod_meta if isinstance(item, dict)}

    for imod_index, imod_id in enumerate(imod_ids):
        item = meta_by_id[imod_id]
        imod_name = item.get("name", imod_id)
        regulator = item.get("linked_regulator") or "Uncharacterized"
        category = item.get("category", "General")
        overlap_count = sum(gene in expression.columns for gene in weight_data[imod_id].get("genes", {}))
        for condition_index, row in de_frame.iterrows():
            score = float(z_activities[imod_index, condition_index])
            act_rows.append((
                imod_id, imod_name, regulator, category, str(row[condition_col]), score,
                int(abs(score) >= 2.0), str(row[sample_col]),
                float(raw_activities[imod_index, condition_index]),
                "least_squares_projection_zscore", os.path.basename(DE_MATRIX_PATH), 0, overlap_count,
            ))

        overlap = weight_data[imod_id].get("regulon_overlap")
        if isinstance(overlap, dict) and overlap.get("regulator"):
            overlap_rows.append((
                imod_id, overlap["regulator"], overlap["regulator"],
                float(overlap.get("precision", 0.0)), float(overlap.get("recall", 0.0)),
                float(overlap.get("f1_score", 0.0)), int(overlap.get("overlap_size", 0)), None,
            ))

    cursor.execute("DELETE FROM imodulon_condition_activities")
    cursor.executemany("""
        INSERT INTO imodulon_condition_activities
        (imodulon_id, imodulon_name, linked_regulator, category, condition_name,
         activity_score, is_significant, sample_id, raw_activity, activity_method,
         data_source, is_measured, gene_overlap_count)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, act_rows)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_act_cond ON imodulon_condition_activities(condition_name);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_act_imod ON imodulon_condition_activities(imodulon_id);")

    # 2. Create imodulon_regulon_overlaps table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS imodulon_regulon_overlaps (
            overlap_id INTEGER PRIMARY KEY AUTOINCREMENT,
            imodulon_id TEXT,
            tf_locus TEXT,
            tf_name TEXT,
            precision REAL,
            recall REAL,
            f1_score REAL,
            overlapping_genes_count INT,
            overlap_gene_list TEXT
        )
    """)

    cursor.execute("DELETE FROM imodulon_regulon_overlaps")
    cursor.executemany("""
        INSERT INTO imodulon_regulon_overlaps
        (imodulon_id, tf_locus, tf_name, precision, recall, f1_score, overlapping_genes_count, overlap_gene_list)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?)
    """, overlap_rows)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ovl_imod ON imodulon_regulon_overlaps(imodulon_id);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ovl_tf ON imodulon_regulon_overlaps(tf_locus);")

    conn.commit()
    conn.close()

    print(f"SUCCESS: Created imodulon_condition_activities ({len(act_rows)} entries) & imodulon_regulon_overlaps ({len(overlap_rows)} entries) in {DB_PATH}")

if __name__ == "__main__":
    import_ica_condition_regulons()
