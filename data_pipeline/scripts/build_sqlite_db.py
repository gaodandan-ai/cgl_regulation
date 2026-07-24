#!/usr/bin/env python3
"""
build_sqlite_db.py
===================
Consolidates all static CSV, JSON, TSV reference files into a single, indexed
SQLite database: data/reference/cgl_regulation.db.

Usage:
    python data_pipeline/scripts/build_sqlite_db.py
"""

import os
import sys
import json
import sqlite3
import pandas as pd
from pathlib import Path

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF_DIR = os.path.join(ROOT_DIR, "data", "reference")
ANALYSIS_DIR = os.path.join(ROOT_DIR, "analysis_output")
DB_PATH = os.path.join(REF_DIR, "cgl_regulation.db")

def build_database():
    print(f"Building SQLite database at: {DB_PATH}")
    if os.path.exists(DB_PATH):
        os.remove(DB_PATH)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # Enable WAL mode for high performance
    cursor.execute("PRAGMA journal_mode=WAL;")
    cursor.execute("PRAGMA synchronous=NORMAL;")

    # 1. Gene Mappings & Canonical Locus Map
    print("[1/13] Processing gene_mapping.csv...")
    gm_path = os.path.join(REF_DIR, "gene_mapping.csv")
    if os.path.exists(gm_path):
        df_gm = pd.read_csv(gm_path)
        # Remove byte-for-byte duplicate mapping rows while preserving genuine
        # many-to-one aliases (for example insertion elements and legacy loci).
        df_gm = df_gm.drop_duplicates().reset_index(drop=True)
        df_gm.to_sql("gene_mappings", conn, if_exists="replace", index=False)

        # Build canonical locus map table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS canonical_locus_map (
                alias TEXT PRIMARY KEY,
                canonical_cg TEXT,
                canonical_cgl TEXT,
                gene_name TEXT,
                product TEXT
            )
        """)

        mappings = []
        for _, row in df_gm.iterrows():
            cgl = str(row.get("cgl_locus", "") or "").strip()
            cg = str(row.get("cg_locus", "") or "").strip()
            name = str(row.get("gene_name", "") or "").strip()
            prod = str(row.get("product", "") or "").strip()

            canonical_cg = cg if cg else (cgl if cgl else name)
            canonical_cgl = cgl if cgl else (cg if cg else name)

            if cg:
                mappings.append((cg, canonical_cg, canonical_cgl, name, prod))
                mappings.append((cg.lower(), canonical_cg, canonical_cgl, name, prod))
            if cgl:
                mappings.append((cgl, canonical_cg, canonical_cgl, name, prod))
                mappings.append((cgl.lower(), canonical_cg, canonical_cgl, name, prod))
            if name:
                mappings.append((name, canonical_cg, canonical_cgl, name, prod))
                mappings.append((name.lower(), canonical_cg, canonical_cgl, name, prod))

        cursor.executemany("""
            INSERT OR IGNORE INTO canonical_locus_map (alias, canonical_cg, canonical_cgl, gene_name, product)
            VALUES (?, ?, ?, ?, ?)
        """, mappings)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_canonical_alias ON canonical_locus_map(alias);")

    # 2. Regulations
    print("[2/13] Processing regulations.csv...")
    reg_path = os.path.join(REF_DIR, "regulations.csv")
    if os.path.exists(reg_path):
        df_reg = pd.read_csv(reg_path)
        df_reg.to_sql("regulations", conn, if_exists="replace", index=False)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_reg_tf ON regulations(TF_locusTag);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_reg_tg ON regulations(TG_locusTag);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_reg_tf_name ON regulations(TF_name);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_reg_tg_name ON regulations(TG_name);")

    # 3. Evidence Credibility
    print("[3/13] Processing tf_target_credibility.csv...")
    cred_path = os.path.join(ANALYSIS_DIR, "evidence_credibility", "tf_target_credibility.csv")
    if os.path.exists(cred_path):
        df_cred = pd.read_csv(cred_path)
        df_cred.to_sql("evidence_credibility", conn, if_exists="replace", index=False)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_cred_tf_tg ON evidence_credibility(tf, tg);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_cred_tf_name ON evidence_credibility(tf_name);")

    # 4. Essential Genes
    print("[4/13] Processing essential_genes.json...")
    eg_path = os.path.join(REF_DIR, "essential_genes.json")
    if os.path.exists(eg_path):
        with open(eg_path, "r", encoding="utf-8") as f:
            eg_data = json.load(f)
        eg_rows = []
        for k, v in eg_data.items():
            if isinstance(v, dict):
                eg_rows.append((k.lower(), v.get("symbol", ""), v.get("essentiality", "essential"), json.dumps(v)))
            else:
                eg_rows.append((k.lower(), "", str(v), json.dumps(v)))
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS essential_genes (
                locus_tag TEXT PRIMARY KEY,
                symbol TEXT,
                essentiality TEXT,
                details TEXT
            )
        """)
        cursor.executemany("INSERT OR IGNORE INTO essential_genes VALUES (?, ?, ?, ?)", eg_rows)

    # 5. Operons
    print("[5/13] Processing operons.csv...")
    op_path = os.path.join(REF_DIR, "operons.csv")
    if os.path.exists(op_path):
        op_rows = []
        with open(op_path, "r", encoding="utf-8") as f:
            lines = f.readlines()
            for line in lines[1:]: # skip header
                parts = [p.strip() for p in line.strip().split(",") if p.strip()]
                if len(parts) >= 2:
                    op_id = parts[0].lstrip(">")
                    orientation = parts[1]
                    genes = ",".join(parts[2:])
                    op_rows.append((op_id, orientation, genes))

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS operons (
                operon_id TEXT PRIMARY KEY,
                orientation TEXT,
                genes TEXT
            )
        """)
        cursor.executemany("INSERT OR IGNORE INTO operons VALUES (?, ?, ?)", op_rows)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_operon_id ON operons(operon_id);")

    # 6. BRENDA kcat
    print("[6/13] Processing brenda_kcat_mappings.json...")
    bk_path = os.path.join(REF_DIR, "brenda_kcat_mappings.json")
    if os.path.exists(bk_path):
        with open(bk_path, "r", encoding="utf-8") as f:
            bk_data = json.load(f)
        bk_rows = []
        for rxn_id, info in bk_data.items():
            if rxn_id.startswith("_"):
                continue
            kcat = info.get("kcat") or info.get("kcat_max") or 0.0 if isinstance(info, dict) else 0.0
            bk_rows.append((rxn_id, info.get("ec", "") if isinstance(info, dict) else "", float(kcat), info.get("substrate", "") if isinstance(info, dict) else "", json.dumps(info)))
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS brenda_kcat (
                reaction_id TEXT PRIMARY KEY,
                ec_number TEXT,
                kcat_val REAL,
                substrate TEXT,
                details TEXT
            )
        """)
        cursor.executemany("INSERT OR IGNORE INTO brenda_kcat VALUES (?, ?, ?, ?, ?)", bk_rows)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_brenda_rxn ON brenda_kcat(reaction_id);")

    # 7. STRING Interactions
    print("[7/13] Processing string_interactions.json...")
    si_path = os.path.join(REF_DIR, "string_interactions.json")
    if os.path.exists(si_path):
        with open(si_path, "r", encoding="utf-8") as f:
            si_data = json.load(f)
        si_rows = []
        for gene_a, targets in si_data.items():
            if gene_a.startswith("_"):
                continue
            if isinstance(targets, dict):
                for gene_b, score in targets.items():
                    if gene_b.startswith("_"):
                        continue
                    score_val = score.get("score", score) if isinstance(score, dict) else score
                    try:
                        si_rows.append((gene_a.lower(), gene_b.lower(), float(score_val)))
                    except (ValueError, TypeError):
                        pass
            elif isinstance(targets, list):
                for item in targets:
                    if isinstance(item, dict):
                        gb = item.get("partner", "") or item.get("stringId_B", "") or item.get("gene_b", "")
                        sc = item.get("score", 0.0)
                        if gb:
                            try:
                                si_rows.append((gene_a.lower(), str(gb).lower(), float(sc)))
                            except (ValueError, TypeError):
                                pass

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS string_interactions (
                gene_a TEXT,
                gene_b TEXT,
                score REAL
            )
        """)
        cursor.executemany("INSERT INTO string_interactions VALUES (?, ?, ?)", si_rows)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_string_ga ON string_interactions(gene_a);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_string_gb ON string_interactions(gene_b);")

    # 8. Abasy Roles
    print("[8/13] Processing abasy_roles.json...")
    ab_path = os.path.join(REF_DIR, "abasy_roles.json")
    if os.path.exists(ab_path):
        with open(ab_path, "r", encoding="utf-8") as f:
            ab_data = json.load(f)
        ab_rows = []
        for locus, role in ab_data.items():
            role_str = role.get("systemic_role", str(role)) if isinstance(role, dict) else str(role)
            ab_rows.append((locus.lower(), role_str, json.dumps(role) if isinstance(role, dict) else role_str))
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS abasy_roles (
                locus_tag TEXT PRIMARY KEY,
                systemic_role TEXT,
                details TEXT
            )
        """)
        cursor.executemany("INSERT OR IGNORE INTO abasy_roles VALUES (?, ?, ?)", ab_rows)

    # 9. Rhea & ChEBI Mappings
    print("[9/13] Processing rhea & chebi mappings...")
    rhea_path = os.path.join(REF_DIR, "rhea_mappings.json")
    if os.path.exists(rhea_path):
        with open(rhea_path, "r", encoding="utf-8") as f:
            rhea_data = json.load(f)
        rhea_rows = [(k, json.dumps(v)) for k, v in rhea_data.items()]
        cursor.execute("CREATE TABLE IF NOT EXISTS rhea_mappings (rxn_id TEXT PRIMARY KEY, links TEXT)")
        cursor.executemany("INSERT OR IGNORE INTO rhea_mappings VALUES (?, ?)", rhea_rows)

    chebi_path = os.path.join(REF_DIR, "chebi_mappings.json")
    if os.path.exists(chebi_path):
        with open(chebi_path, "r", encoding="utf-8") as f:
            chebi_data = json.load(f)
        chebi_rows = [(k, json.dumps(v)) for k, v in chebi_data.items()]
        cursor.execute("CREATE TABLE IF NOT EXISTS chebi_mappings (met_id TEXT PRIMARY KEY, details TEXT)")
        cursor.executemany("INSERT OR IGNORE INTO chebi_mappings VALUES (?, ?)", chebi_rows)

    # 10. COG Annotations
    print("[10/13] Processing cog_annotations.json...")
    cog_path = os.path.join(REF_DIR, "cog_annotations.json")
    if os.path.exists(cog_path):
        with open(cog_path, "r", encoding="utf-8") as f:
            cog_data = json.load(f)
        cog_rows = []
        for locus, info in cog_data.items():
            if isinstance(info, dict):
                cog_rows.append((locus.lower(), info.get("cog_id", ""), info.get("category", ""), info.get("description", "")))
            else:
                cog_rows.append((locus.lower(), "", str(info), ""))
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS cog_annotations (
                locus_tag TEXT PRIMARY KEY,
                cog_id TEXT,
                category TEXT,
                description TEXT
            )
        """)
        cursor.executemany("INSERT OR IGNORE INTO cog_annotations VALUES (?, ?, ?, ?)", cog_rows)

    # 11. Network Centrality
    print("[11/13] Processing network_centrality.json...")
    nc_path = os.path.join(REF_DIR, "network_centrality.json")
    if os.path.exists(nc_path):
        with open(nc_path, "r", encoding="utf-8") as f:
            nc_data = json.load(f)
        nc_rows = []
        # network_centrality.json stores per-gene records under ``nodes``.
        # Older code iterated the top-level metadata object and silently
        # created bogus rows named ``_meta`` and ``nodes``.
        node_data = nc_data.get("nodes", nc_data) if isinstance(nc_data, dict) else {}
        for locus, info in node_data.items():
            if isinstance(info, dict):
                in_degree = int(info.get("in_degree", 0) or 0)
                out_degree = int(info.get("out_degree", 0) or 0)
                nc_rows.append((
                    locus.lower(),
                    int(info.get("degree", in_degree + out_degree) or 0),
                    in_degree,
                    out_degree,
                    float(info.get("betweenness", 0.0)),
                    float(info.get("closeness", 0.0)),
                    float(info.get("pagerank", 0.0))
                ))
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS network_centrality (
                locus_tag TEXT PRIMARY KEY,
                degree INT,
                in_degree INT,
                out_degree INT,
                betweenness REAL,
                closeness REAL,
                pagerank REAL
            )
        """)
        cursor.executemany("INSERT OR IGNORE INTO network_centrality VALUES (?, ?, ?, ?, ?, ?, ?)", nc_rows)

    # 12. Full-Text Search (FTS5) for Literature RAG
    print("[12/13] Setting up Literature FTS5 index...")
    lit_path = os.path.join(REF_DIR, "literature_cache.json")
    cursor.execute("""
        CREATE VIRTUAL TABLE IF NOT EXISTS literature_fts USING fts5(
            doc_id,
            gene_locus,
            title,
            abstract,
            pmid
        );
    """)

    if os.path.exists(lit_path):
        with open(lit_path, "r", encoding="utf-8") as f:
            lit_data = json.load(f)
        lit_rows = []
        if isinstance(lit_data, dict):
            files_dict = lit_data.get("files", {})
            if files_dict and isinstance(files_dict, dict):
                for filename, fmeta in files_dict.items():
                    chunks = fmeta.get("chunks", [])
                    for c_idx, chunk in enumerate(chunks):
                        text = chunk.get("text", "") if isinstance(chunk, dict) else str(chunk)
                        lit_rows.append((
                            f"{filename}_{c_idx}",
                            filename,
                            filename,
                            text,
                            ""
                        ))
            for locus, articles in lit_data.items():
                if locus == "files":
                    continue
                if isinstance(articles, list):
                    for idx, art in enumerate(articles):
                        if isinstance(art, dict):
                            lit_rows.append((
                                f"{locus}_{idx}",
                                locus,
                                art.get("title", ""),
                                art.get("abstract", art.get("text", "")),
                                str(art.get("pmid", ""))
                            ))
                elif isinstance(articles, dict):
                    lit_rows.append((
                        locus,
                        locus,
                        articles.get("title", ""),
                        articles.get("abstract", articles.get("text", "")),
                        str(articles.get("pmid", ""))
                    ))
        elif isinstance(lit_data, list):
            for idx, art in enumerate(lit_data):
                locus = art.get("locus_tag", art.get("gene", "unknown"))
                lit_rows.append((
                    f"doc_{idx}",
                    locus,
                    art.get("title", ""),
                    art.get("abstract", ""),
                    str(art.get("pmid", ""))
                ))

        cursor.executemany("INSERT INTO literature_fts (doc_id, gene_locus, title, abstract, pmid) VALUES (?, ?, ?, ?, ?)", lit_rows)

    # 13. ChIP-seq Peaks & Regulations Ingestion
    print("[13/14] Processing ChIP-seq peaks and regulations...")
    try:
        from data_pipeline.scripts.process_chipseq_peaks_integration import process_and_ingest
    except ImportError:
        from process_chipseq_peaks_integration import process_and_ingest
    conn.commit()
    conn.close()
    process_and_ingest(Path(DB_PATH))

    # Re-open connection for final optimization
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 14. Final Commit & Optimization
    print("[14/14] Committing changes & optimizing DB...")
    conn.commit()
    cursor.execute("VACUUM;")
    cursor.execute("ANALYZE;")
    conn.close()

    db_size_mb = os.path.getsize(DB_PATH) / (1024 * 1024)
    print(f"SUCCESS: Created {DB_PATH} ({db_size_mb:.2f} MB)")

if __name__ == "__main__":
    build_database()
