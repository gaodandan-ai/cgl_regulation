#!/usr/bin/env python3
"""
import_all_omics.py
===================
Imports the remaining 5 omics and machine learning datasets into data/reference/cgl_regulation.db:
  1. imodulons & imodulon_gene_weights (ICA Co-Expression Modules)
  2. tf_gene_rf_scores (Random Forest ML Edge Confidence Predictions)
  3. reaction_thermodynamics (eQuilibrATOR Gibbs Free Energy & MW)
  4. tf_hierarchy_rankings (TF 3-Tier Pyramid Hierarchy Rankings)
  5. network_rewired_edges (Cross-strain Evolutionary Rewired Network Edges)

Usage:
    python data_pipeline/scripts/import_all_omics.py
"""

import os
import sys
import json
import sqlite3
import pandas as pd

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF_DIR = os.path.join(ROOT_DIR, "data", "reference")
ANALYSIS_DIR = os.path.join(ROOT_DIR, "analysis_output")
DB_PATH = os.path.join(REF_DIR, "cgl_regulation.db")

def import_all_omics():
    print(f"Importing all omics & ML datasets into: {DB_PATH}")
    if not os.path.exists(DB_PATH):
        print("ERROR: cgl_regulation.db not found. Build DB first.")
        sys.exit(1)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. iModulon Co-Expression Modules & Gene Weights
    print("[1/5] Importing iModulons & Gene Weights...")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS imodulons (
            imodulon_id TEXT PRIMARY KEY,
            name TEXT,
            size INT,
            explanation TEXT,
            details TEXT
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS imodulon_gene_weights (
            imodulon_id TEXT,
            gene_locus TEXT,
            weight REAL,
            PRIMARY KEY (imodulon_id, gene_locus)
        )
    """)

    imod_meta_path = os.path.join(REF_DIR, "imodulon", "imodulon_metadata.json")
    if os.path.exists(imod_meta_path):
        with open(imod_meta_path, "r", encoding="utf-8") as f:
            imod_meta = json.load(f)
        imod_rows = []
        if isinstance(imod_meta, list):
            for v in imod_meta:
                if isinstance(v, dict):
                    imod_id = v.get("id") or v.get("name", "iM_unknown")
                    imod_rows.append((
                        imod_id,
                        v.get("name", imod_id),
                        int(v.get("gene_count", v.get("size", 0))),
                        v.get("description", v.get("explanation", "")),
                        json.dumps(v)
                    ))
        elif isinstance(imod_meta, dict):
            for k, v in imod_meta.items():
                if isinstance(v, dict):
                    imod_rows.append((
                        k,
                        v.get("name", k),
                        int(v.get("size", 0)),
                        v.get("explanation", ""),
                        json.dumps(v)
                    ))
                else:
                    imod_rows.append((k, k, 0, str(v), json.dumps(v)))

        cursor.executemany("INSERT OR IGNORE INTO imodulons VALUES (?, ?, ?, ?, ?)", imod_rows)

    imod_weights_path = os.path.join(REF_DIR, "imodulon", "imodulon_gene_weights.json")
    if os.path.exists(imod_weights_path):
        with open(imod_weights_path, "r", encoding="utf-8") as f:
            imod_weights = json.load(f)
        w_rows = []
        for imod_id, imod_obj in imod_weights.items():
            if isinstance(imod_obj, dict):
                genes_dict = imod_obj.get("genes", imod_obj)
                if isinstance(genes_dict, dict):
                    for gene, weight in genes_dict.items():
                        try:
                            w_rows.append((imod_id, gene.lower(), float(weight)))
                        except (ValueError, TypeError):
                            pass
        cursor.executemany("INSERT OR IGNORE INTO imodulon_gene_weights VALUES (?, ?, ?)", w_rows)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_imod_gene ON imodulon_gene_weights(gene_locus);")

    # 2. Random Forest ML Edge Confidence Predictions
    print("[2/5] Importing Random Forest Edge Scores...")
    rf_scores_path = os.path.join(REF_DIR, "edge_confidence", "tf_gene_edge_scores.csv")
    if os.path.exists(rf_scores_path):
        df_rf = pd.read_csv(rf_scores_path)
        df_rf.to_sql("tf_gene_rf_scores", conn, if_exists="replace", index=False)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_rf_tf ON tf_gene_rf_scores(tf_locus);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_rf_tg ON tf_gene_rf_scores(target_locus);")

    # 3. Reaction Thermodynamics & Gibbs Free Energy
    print("[3/5] Importing Reaction Thermodynamics & MW...")
    cursor.execute("DROP TABLE IF EXISTS reaction_thermodynamics")
    cursor.execute("""
        CREATE TABLE reaction_thermodynamics (
            reaction_id TEXT PRIMARY KEY,
            dgr_prime_o REAL,
            dgr_uncertainty REAL,
            mw REAL,
            kcat REAL,
            dgr_status TEXT NOT NULL CHECK(dgr_status IN ('available', 'missing')),
            missing_reason TEXT,
            thermo_source TEXT,
            thermo_confidence TEXT,
            direction_locked TEXT,
            in_model INTEGER NOT NULL CHECK(in_model IN (0, 1)),
            enzyme_parameter_source TEXT,
            details TEXT
        )
    """)

    thermo_path = os.path.join(REF_DIR, "thermo_dgr_data.json")
    thermo_rows = []
    if os.path.exists(thermo_path):
        with open(thermo_path, "r", encoding="utf-8") as f:
            thermo_data = json.load(f)
        # thermo_dgr_data.json stores reaction records under ``reactions``.
        # Iterating the top-level object used to create one bogus aggregate
        # row called ``reactions``.
        reactions = thermo_data.get("reactions", thermo_data) if isinstance(thermo_data, dict) else {}
        for rxn_id, info in reactions.items():
            if rxn_id.startswith("_"):
                continue
            if isinstance(info, dict):
                dgr = info.get("dgr_prime_0", info.get("dgr_prime_o", info.get("dgr")))
                unc = info.get("uncertainty", info.get("dgr_uncertainty"))
                if unc is None and info.get("dgr_prime_min") is not None and info.get("dgr_prime_max") is not None:
                    unc = (float(info["dgr_prime_max"]) - float(info["dgr_prime_min"])) / 2.0
                mw = info.get("mw")
                kcat = info.get("kcat")
                status = "available" if dgr is not None else "missing"
                thermo_rows.append((
                    rxn_id,
                    float(dgr) if dgr is not None else None,
                    float(unc) if unc is not None else None,
                    float(mw) if mw is not None else None,
                    float(kcat) if kcat is not None else None,
                    status,
                    info.get("note") if status == "missing" else None,
                    info.get("source"),
                    info.get("confidence"),
                    info.get("direction_locked"),
                    int(bool(info.get("in_model", True))),
                    "thermodynamics JSON" if mw is not None or kcat is not None else None,
                    json.dumps(info),
                ))

    cursor.executemany("""
        INSERT OR REPLACE INTO reaction_thermodynamics
        (reaction_id, dgr_prime_o, dgr_uncertainty, mw, kcat, dgr_status,
         missing_reason, thermo_source, thermo_confidence, direction_locked,
         in_model, enzyme_parameter_source, details)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, thermo_rows)
    # Link enzyme-constrained model parameters when available. This project
    # file contains reaction-level mean molecular mass and kcat values.
    enzyme_path = os.path.join(
        REF_DIR, "model", "ecCGL1-main", "ecCGL1-main", "iCW773_get_data", "reaction_kcat_MW.csv"
    )
    if os.path.exists(enzyme_path):
        df_enzyme = pd.read_csv(enzyme_path)
        reaction_col = df_enzyme.columns[0]
        enzyme_updates = []
        for _, row in df_enzyme.iterrows():
            rxn_id = str(row.get(reaction_col, "") or "").strip()
            if not rxn_id:
                continue
            kcat = row.get("kcat")
            mw = row.get("MW")
            enzyme_updates.append((
                float(mw) if pd.notnull(mw) else None,
                float(kcat) if pd.notnull(kcat) else None,
                rxn_id,
            ))
        cursor.executemany("""
            UPDATE reaction_thermodynamics
            SET mw=?, kcat=?, enzyme_parameter_source='ecCGL1 reaction_kcat_MW.csv'
            WHERE reaction_id=?
        """, enzyme_updates)
    # Fill remaining exact-ID kcat values from the consolidated BRENDA/DLKcat layer.
    cursor.execute("""
        UPDATE reaction_thermodynamics
        SET kcat=(SELECT b.kcat_val FROM brenda_kcat b
                  WHERE b.reaction_id=reaction_thermodynamics.reaction_id),
            enzyme_parameter_source='BRENDA/DLKcat exact reaction ID'
        WHERE kcat IS NULL AND EXISTS (
            SELECT 1 FROM brenda_kcat b WHERE b.reaction_id=reaction_thermodynamics.reaction_id
        )
    """)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_thermo_rxn ON reaction_thermodynamics(reaction_id);")

    # 4. TF 3-Tier Hierarchy Rankings
    print("[4/5] Importing TF Hierarchy Rankings...")
    tf_hier_path = os.path.join(ANALYSIS_DIR, "tf_hierarchy", "tf_hierarchy_rankings.csv")
    if os.path.exists(tf_hier_path):
        df_tf = pd.read_csv(tf_hier_path)
        df_tf.to_sql("tf_hierarchy_rankings", conn, if_exists="replace", index=False)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_tfh_locus ON tf_hierarchy_rankings(locus);")

    # 5. Cross-Strain Evolutionary Rewired Network Edges
    print("[5/5] Importing Cross-Strain Rewired Network Edges...")
    rewired_path = os.path.join(ANALYSIS_DIR, "cross_strain", "rewired_edges.csv")
    if os.path.exists(rewired_path):
        df_rew = pd.read_csv(rewired_path)
        df_rew.to_sql("network_rewired_edges", conn, if_exists="replace", index=False)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_rew_tf ON network_rewired_edges(tf_cgl);")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_rew_tg ON network_rewired_edges(tg_cgl);")

    conn.commit()
    conn.close()

    db_size_mb = os.path.getsize(DB_PATH) / (1024 * 1024)
    print(f"SUCCESS: Imported all omics & ML datasets into {DB_PATH} ({db_size_mb:.2f} MB)")

if __name__ == "__main__":
    import_all_omics()
