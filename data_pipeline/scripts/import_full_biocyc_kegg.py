#!/usr/bin/env python3
"""
import_full_biocyc_kegg.py
==========================
Imports full 250+ BioCyc/KEGG pathways and 1,915 gene-pathway mappings
from data/reference/kegg_cache/kegg_cgl_cgb.json into cgl_regulation.db.

Creates:
  - biocyc_kegg_pathways (pathway_id, pathway_name, category, gene_count, source)
  - gene_pathway_mappings (gene_locus, pathway_id, pathway_name)
"""

import os
import sys
import json
import sqlite3
import re

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF_DIR = os.path.join(ROOT_DIR, "data", "reference")
DB_PATH = os.path.join(REF_DIR, "cgl_regulation.db")
KEGG_JSON_PATH = os.path.join(REF_DIR, "kegg_cache", "kegg_cgl_cgb.json")
KEGG_HIERARCHY_PATH = os.path.join(REF_DIR, "kegg_cache", "kegg_pathway_hierarchy.txt")


def load_kegg_categories():
    categories = {}
    current_group = "Unclassified"
    if not os.path.exists(KEGG_HIERARCHY_PATH):
        return categories
    with open(KEGG_HIERARCHY_PATH, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("B  "):
                current_group = line[3:].strip()
            elif line.startswith("C"):
                match = re.match(r"C\s+(\d{5})\s+", line)
                if match:
                    categories[match.group(1)] = current_group
    return categories

def import_full_biocyc_kegg():
    print(f"Importing full BioCyc/KEGG pathways into: {DB_PATH}")
    if not os.path.exists(KEGG_JSON_PATH):
        print(f"ERROR: {KEGG_JSON_PATH} not found.")
        sys.exit(1)

    with open(KEGG_JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)

    pathway_names = data.get("pathway_names", {})
    gene_to_pathways = data.get("gene_to_pathways", {})
    pathway_to_genes = data.get("pathway_to_genes", {})
    kegg_categories = load_kegg_categories()

    print(f"Loaded {len(pathway_names)} pathway names, {len(gene_to_pathways)} gene mappings, {len(pathway_to_genes)} pathway-gene sets.")

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. Update biocyc_kegg_pathways table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS biocyc_kegg_pathways (
            pathway_id TEXT PRIMARY KEY,
            pathway_name TEXT,
            category TEXT,
            gene_list TEXT,
            source TEXT
        )
    """)
    columns = {row[1] for row in cursor.execute("PRAGMA table_info(biocyc_kegg_pathways)")}
    if "category_source" not in columns:
        cursor.execute("ALTER TABLE biocyc_kegg_pathways ADD COLUMN category_source TEXT")

    pathway_rows = []
    for pid, pname in pathway_names.items():
        genes = pathway_to_genes.get(pid, [])
        gene_str = ",".join(genes)
        code_match = re.search(r"(\d{5})$", pid)
        category = kegg_categories.get(code_match.group(1), "Unclassified") if code_match else "Unclassified"
        pathway_rows.append((pid, pname, category, gene_str, "BioCyc/KEGG Full", "KEGG BRITE br08901"))

    cursor.executemany("""
        INSERT OR REPLACE INTO biocyc_kegg_pathways
        (pathway_id, pathway_name, category, gene_list, source, category_source)
        VALUES (?, ?, ?, ?, ?, ?)
    """, pathway_rows)

    # 2. Create gene_pathway_mappings relational table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gene_pathway_mappings (
            gene_locus TEXT,
            pathway_id TEXT,
            pathway_name TEXT,
            PRIMARY KEY (gene_locus, pathway_id)
        )
    """)

    g2p_rows = []
    for gene_locus, pids in gene_to_pathways.items():
        g_clean = gene_locus.lower().strip()
        for pid in pids:
            pname = pathway_names.get(pid, pid)
            g2p_rows.append((g_clean, pid, pname))

    cursor.executemany("""
        INSERT OR IGNORE INTO gene_pathway_mappings (gene_locus, pathway_id, pathway_name)
        VALUES (?, ?, ?)
    """, g2p_rows)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_g2p_gene ON gene_pathway_mappings(gene_locus);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_g2p_pathway ON gene_pathway_mappings(pathway_id);")

    conn.commit()
    conn.close()

    print(f"SUCCESS: Imported {len(pathway_rows)} BioCyc/KEGG pathways and {len(g2p_rows)} gene-pathway links into {DB_PATH}")

if __name__ == "__main__":
    import_full_biocyc_kegg()
