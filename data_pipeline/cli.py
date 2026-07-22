#!/usr/bin/env python3
"""
cli.py
======
Unified Command-Line Interface (CLI) for data pipeline operations, database building,
integrity auditing, and statistics reporting.

Usage:
    python -m data_pipeline.cli build-db
    python -m data_pipeline.cli audit-db
    python -m data_pipeline.cli stats
"""

import os
import sys
import argparse
import sqlite3
import shutil
from datetime import datetime

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF_DIR = os.path.join(ROOT_DIR, "data", "reference")
DB_PATH = os.path.join(REF_DIR, "cgl_regulation.db")

def cmd_build_db():
    """Build every database layer, restoring the previous DB on failure."""
    print("=== Data Pipeline: Building Complete SQLite Database ===")
    backup_path = None
    if os.path.exists(DB_PATH):
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup_path = f"{DB_PATH}.prebuild_{stamp}.bak"
        shutil.copy2(DB_PATH, backup_path)
        print(f"Safety backup: {backup_path}")
    try:
        from data_pipeline.scripts.build_sqlite_db import build_database
        from data_pipeline.scripts.enrich_database import enrich_database
        from data_pipeline.scripts.import_refseq_genome import import_refseq_genome
        from data_pipeline.scripts.import_all_omics import import_all_omics
        from data_pipeline.scripts.import_full_biocyc_kegg import import_full_biocyc_kegg
        from data_pipeline.scripts.process_collectf_tfbs import deduplicate_and_merge
        from data_pipeline.scripts.import_rfam_ncrna import import_rfam_ncrna
        from data_pipeline.scripts.import_ica_condition_regulons import import_ica_condition_regulons
        from data_pipeline.scripts.import_integrated_omics import import_integrated_omics
        from data_pipeline.scripts.build_condition_harmonization import build_condition_harmonization
        from data_pipeline.scripts.build_condition_regulatory_scores import build_condition_regulatory_scores
        from data_pipeline.scripts.build_oxygen_regulatory_scores import build_oxygen_regulatory_scores
        from data_pipeline.scripts.build_carbon_regulatory_scores import build_carbon_regulatory_scores
        from data_pipeline.scripts.build_nitrogen_regulatory_scores import build_nitrogen_regulatory_scores
        from data_pipeline.scripts.build_stress_regulatory_scores import build_stress_regulatory_scores
        from data_pipeline.scripts.build_intervention_priorities import build_intervention_priorities
        from data_pipeline.scripts.finalize_database import finalize_database

        build_database()
        enrich_database()
        import_refseq_genome()
        import_all_omics()
        import_full_biocyc_kegg()
        deduplicate_and_merge()
        import_rfam_ncrna()
        import_ica_condition_regulons()
        import_integrated_omics()
        build_condition_harmonization()
        finalize_database()
        build_condition_regulatory_scores()
        build_oxygen_regulatory_scores()
        build_carbon_regulatory_scores()
        build_nitrogen_regulatory_scores()
        build_stress_regulatory_scores()
        build_intervention_priorities()
    except BaseException:
        if backup_path and os.path.exists(backup_path):
            shutil.copy2(backup_path, DB_PATH)
            print(f"Build failed; restored database from {backup_path}")
        raise
    print("Complete database build finished successfully.")

def cmd_audit_db():
    print(f"=== Data Pipeline: Auditing Database ({DB_PATH}) ===")
    if not os.path.exists(DB_PATH):
        print("ERROR: Database file does not exist. Run build-db first.")
        sys.exit(1)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tables = [r[0] for r in cursor.fetchall()]

    print(f"Found {len(tables)} tables in database:\n" + "-" * 50)
    total_rows = 0
    for tbl in sorted(tables):
        cursor.execute(f"SELECT COUNT(*) FROM `{tbl}`")
        count = cursor.fetchone()[0]
        total_rows += count
        print(f"  - Table `{tbl:<24}` : {count:>7,d} rows")

    cursor.execute("SELECT COUNT(*) FROM sqlite_master WHERE type='index';")
    index_count = cursor.fetchone()[0]
    print("-" * 50)
    print(f"Total Rows: {total_rows:,d} | Indexes: {index_count} | Size: {os.path.getsize(DB_PATH)/(1024*1024):.2f} MB")
    conn.close()

def main():
    parser = argparse.ArgumentParser(description="Cgl Regulation Explorer Data Pipeline CLI")
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    subparsers.add_parser("build-db", help="Build SQLite cgl_regulation.db from reference files")
    subparsers.add_parser("audit-db", help="Audit database integrity and count table rows")
    subparsers.add_parser("stats", help="Show database summary statistics")

    args = parser.parse_args()
    if args.command == "build-db":
        cmd_build_db()
    elif args.command in ("audit-db", "stats"):
        cmd_audit_db()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
