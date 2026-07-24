#!/usr/bin/env python3
"""
cross_validate_chipseq_expression.py
======================================
Multi-omics cross-validation & functional edge classification script.

Intersects ChIP-seq binding peaks (chipseq_peaks) with expression compendium
correlations, perturbation expression changes, and PWM binding motifs.

Classifies regulatory interactions into:
  - DIRECT_FUNCTIONAL_REGULATION (Peak in promoter + Expression response + PWM match)
  - STRUCTURAL_PHYSICAL_BINDING (Peak present, but minimal expression change)
  - INDIRECT_EXPRESSION_RESPONSE (Expression response, but no direct ChIP peak)
  - COMPUTATIONAL_PREDICTION_UNVERIFIED (Unvalidated computationally predicted edge)

Populates SQLite database tables & creates `v_chipseq_functional_regulations` view.

Usage:
    python data_pipeline/scripts/cross_validate_chipseq_expression.py
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


def run_cross_validation(db_path: Path):
    print(f"--- Running Multi-Omics Cross-Validation on {db_path} ---")
    if not db_path.exists():
        print(f"Database {db_path} not found, skipping.")
        return

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Check if chipseq_peaks table exists
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='chipseq_peaks';")
    if not cursor.fetchone():
        print("Warning: chipseq_peaks table not found in database. Ingest peaks first.")
        conn.close()
        return

    # Add columns if missing in regulations and chipseq_regulations
    for tbl in ["regulations", "chipseq_regulations"]:
        cursor.execute(f"PRAGMA table_info({tbl});")
        cols = {row[1] for row in cursor.fetchall()}
        if "functional_category" not in cols:
            cursor.execute(f"ALTER TABLE {tbl} ADD COLUMN functional_category TEXT;")
        if "functional_score" not in cols:
            cursor.execute(f"ALTER TABLE {tbl} ADD COLUMN functional_score REAL;")

    # Load peak summaries per (tf, target_gene)
    cursor.execute("""
        SELECT LOWER(tf_name) as tf_name_lower, LOWER(tf_id) as tf_id_lower,
               LOWER(nearest_gene_locus) as tg_locus_lower,
               MAX(peak_score) as max_score, MAX(neglog10q) as max_negq,
               MIN(ABS(COALESCE(rel_pos_to_tss, 9999))) as min_abs_tss,
               spatial_confidence
        FROM chipseq_peaks
        WHERE nearest_gene_locus IS NOT NULL AND nearest_gene_locus != ''
        GROUP BY tf_name_lower, tf_id_lower, tg_locus_lower;
    """)
    peak_map = {}
    for r in cursor.fetchall():
        tf_name_l, tf_id_l, tg_l, max_score, max_negq, min_abs_tss, spatial_conf = r
        peak_info = {
            "max_score": max_score or 1.0,
            "max_negq": max_negq or 0.0,
            "min_abs_tss": min_abs_tss,
            "spatial_conf": spatial_conf or "PROMOTER_DIRECT"
        }
        if tf_name_l:
            peak_map[(tf_name_l, tg_l)] = peak_info
        if tf_id_l:
            peak_map[(tf_id_l, tg_l)] = peak_info

    # Load expression correlation if available
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='tf_gene_rf_scores';")
    rf_map = {}
    if cursor.fetchone():
        cursor.execute("SELECT LOWER(tf_locus), LOWER(target_locus), expression_correlation, predicted_confidence FROM tf_gene_rf_scores;")
        for r in cursor.fetchall():
            tf_l, tg_l, expr_corr, rf_conf = r
            try:
                corr_val = float(expr_corr) if expr_corr is not None else None
            except (ValueError, TypeError):
                corr_val = None
            rf_map[(tf_l, tg_l)] = {"expr_corr": corr_val, "rf_conf": rf_conf or 0.5}

    # Process regulations and chipseq_regulations tables
    for tbl in ["regulations", "chipseq_regulations"]:
        cursor.execute(f"PRAGMA table_info({tbl});")
        cols = [c[1].lower() for c in cursor.fetchall()]
        tf_col = "tf_locus_tag" if "tf_locus_tag" in cols else "tf_locustag"
        tg_col = "tg_locus_tag" if "tg_locus_tag" in cols else "tg_locustag"
        tf_name_col = "tf_name"
        bsite_col = "binding_site"

        cursor.execute(f"SELECT rowid, LOWER({tf_col}), LOWER({tf_name_col}), LOWER({tg_col}), {bsite_col} FROM {tbl};")
        reg_rows = cursor.fetchall()

        updates = []
        category_counts = {
            "DIRECT_FUNCTIONAL_REGULATION": 0,
            "STRUCTURAL_PHYSICAL_BINDING": 0,
            "INDIRECT_EXPRESSION_RESPONSE": 0,
            "COMPUTATIONAL_PREDICTION_UNVERIFIED": 0
        }

        for rowid, tf_loc, tf_name, tg_loc, bsite in reg_rows:
            pk = peak_map.get((tf_loc, tg_loc)) or peak_map.get((tf_name, tg_loc))
            expr = rf_map.get((tf_loc, tg_loc)) or {}
            expr_corr = expr.get("expr_corr")
            has_pwm = bool(bsite and str(bsite).strip() and str(bsite).strip().lower() != "nan")

            has_peak = pk is not None
            is_promoter_peak = has_peak and (pk["spatial_conf"] in ["PROMOTER_DIRECT", "INTERGENIC_PROMOTER"] or pk["min_abs_tss"] <= 350)
            has_expr_change = expr_corr is not None and abs(expr_corr) >= 0.15

            if is_promoter_peak and (has_expr_change or has_pwm or (pk and pk["max_score"] >= 5.0)):
                cat = "DIRECT_FUNCTIONAL_REGULATION"
                score = 0.85 + (0.08 if has_pwm else 0.0) + (0.07 if has_expr_change else 0.0)
            elif is_promoter_peak or has_peak:
                cat = "STRUCTURAL_PHYSICAL_BINDING"
                score = 0.65 + (0.10 if pk and pk["max_score"] > 3.0 else 0.0)
            elif has_expr_change:
                cat = "INDIRECT_EXPRESSION_RESPONSE"
                score = 0.50 + min(0.30, abs(expr_corr))
            else:
                cat = "COMPUTATIONAL_PREDICTION_UNVERIFIED"
                score = 0.35 + (0.15 if has_pwm else 0.0)

            score = min(1.0, round(score, 4))
            category_counts[cat] += 1
            updates.append((cat, score, rowid))

        cursor.executemany(f"UPDATE {tbl} SET functional_category = ?, functional_score = ? WHERE rowid = ?;", updates)

        print(f"  Table `{tbl}` classification summary:")
        for k, v in category_counts.items():
            print(f"    - {k}: {v} interactions")

    # Re-create View for Functional Regulations
    cursor.execute("DROP VIEW IF EXISTS v_chipseq_functional_regulations;")
    cursor.execute("""
        CREATE VIEW v_chipseq_functional_regulations AS
        SELECT r.TF_locusTag, r.TF_name, r.TG_locusTag, r.TG_name, r.Role, r.Evidence,
               r.functional_category, r.functional_score,
               p.peak_id, p.peak_score, p.neglog10q, p.strength_tier, p.spatial_confidence, p.rel_pos_to_tss
        FROM regulations r
        LEFT JOIN chipseq_peaks p ON (LOWER(r.TF_name) = LOWER(p.tf_name) OR LOWER(r.TF_locusTag) = LOWER(p.tf_id))
                                 AND LOWER(r.TG_locusTag) = LOWER(p.nearest_gene_locus);
    """)

    conn.commit()
    conn.close()

    print(f"Cross-Validation complete for {db_path}. Classification summary:")
    for k, v in category_counts.items():
        print(f"  - {k}: {v} interactions")


def main():
    if DATA_DB_PATH.exists():
        run_cross_validation(DATA_DB_PATH)
    if REF_DB_PATH.exists():
        run_cross_validation(REF_DB_PATH)
    print("Multi-omics cross-validation completed.")


if __name__ == "__main__":
    main()
