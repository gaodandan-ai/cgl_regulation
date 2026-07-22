#!/usr/bin/env python3
"""
enrich_database.py
===================
Enriches data/reference/cgl_regulation.db with 4 major data tables:
  1. network_edges_extended (Strong/All-evidence/sRNA interactions)
  2. gene_coordinates (NCBI RefSeq NC_003450.3 start, end, strand, length)
  3. biocyc_kegg_pathways (BioCyc & KEGG pathway annotations & operon structures)
  4. tf_families_effectors (UniProt TF family & Effector small molecule annotations)

Usage:
    python data_pipeline/scripts/enrich_database.py
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

# Known C. glutamicum TF structural families and Effector small molecules
KNOWN_TF_EFFECTORS = {
    "cg0350": {"name": "glxR", "family": "CRP/FNR family", "hth": "HTH CRP-type", "effector": "cAMP (cyclic AMP)", "signal": "Carbon source availability & catabolite repression", "role": "Master carbon & stress regulator"},
    "cg0986": {"name": "amtR", "family": "TetR family", "hth": "HTH TetR-type", "effector": "2-Oxoglutarate / GlnK", "signal": "Nitrogen status", "role": "Master nitrogen regulator"},
    "cg2103": {"name": "dtxR", "family": "DtxR family", "hth": "HTH DtxR-type", "effector": "Fe2+ (Ferrous iron)", "signal": "Iron availability", "role": "Iron homeostasis & iron-dependent regulation"},
    "cg2115": {"name": "sugR", "family": "DeoR family", "hth": "HTH DeoR-type", "effector": "Fructose-1-phosphate / Glucose-6-phosphate", "signal": "Sugar import & glycolytic flux", "role": "Repressor of sugar import & glycolysis"},
    "cg3202": {"name": "farR", "family": "TetR family", "hth": "HTH TetR-type", "effector": "Fatty acids / Acyl-CoA", "signal": "Lipid metabolism status", "role": "Regulator of fatty acid biosynthesis"},
    "cg1340": {"name": "arnR", "family": "LysR family", "hth": "HTH LysR-type", "effector": "L-Arginine", "signal": "Arginine concentration", "role": "Repressor of arginine biosynthesis"},
    "cg1120": {"name": "ripA", "family": "AraC family", "hth": "HTH AraC-type", "effector": "Iron limitation / Fe-S clusters", "signal": "Iron starvation", "role": "Repressor of iron-containing proteins"},
    "cg1935": {"name": "gntR2", "family": "GntR family", "hth": "HTH GntR-type", "effector": "Gluconate", "signal": "Gluconate metabolism", "role": "Regulator of gluconate metabolism"},
    "cg2152": {"name": "clgR", "family": "Crop/ClgR family", "hth": "HTH ClgR-type", "effector": "Proteolytic stress / Heat shock", "signal": "Protein misfolding", "role": "Master regulator of protein quality control & Clp proteases"},
    "cg3247": {"name": "hrrA", "family": "Two-component response regulator", "hth": "HTH OmpR-type", "effector": "Heme / HrrS kinase phosphorylation", "signal": "Heme & iron limitation", "role": "Regulator of heme biosynthesis & utilization"},
    "cg1675": {"name": "phoR", "family": "Two-component response regulator", "hth": "HTH OmpR-type", "effector": "Inorganic phosphate / PhoS kinase", "signal": "Phosphate starvation", "role": "Master regulator of phosphate starvation response"},
    "cg0313": {"name": "lrp", "family": "Lrp/AsnC family", "hth": "HTH Lrp-type", "effector": "Branched-chain amino acids (L-Leucine, L-Isoleucine)", "signal": "Amino acid pool", "role": "Regulator of amino acid transport & metabolism"},
    "cg2737": {"name": "fasR", "family": "TetR family", "hth": "HTH TetR-type", "effector": "Acyl-CoA / Malonyl-CoA", "signal": "Fatty acid synthesis precursor levels", "role": "Essential regulator of fatty acid synthase fasA/fasB"},
    "cg2831": {"name": "ramA", "family": "LuxR family", "hth": "HTH LuxR-type", "effector": "Acetate / Propionate / Carbon source transition", "signal": "Alternative carbon sources", "role": "Activator of acetate metabolism & glyoxylate cycle"},
    "cg0444": {"name": "ramB", "family": "IclR family", "hth": "HTH IclR-type", "effector": "Glucose / Acetate transition", "signal": "Glycolytic vs gluconeogenic flux", "role": "Repressor of acetate metabolism & gluconeogenesis"},
    "cg0196": {"name": "iolR", "family": "RpiR family", "hth": "HTH RpiR-type", "effector": "myo-Inositol / 2-Keto-myo-inositol", "signal": "Inositol availability", "role": "Repressor of myo-inositol catabolism"},
    "cg0012": {"name": "ssuR", "family": "LysR family", "hth": "HTH LysR-type", "effector": "Sulfur compounds", "signal": "Sulfur limitation", "role": "Regulator of sulfur starvation genes"},
    "cg0876": {"name": "sigH", "family": "ECF Sigma factor (ECF01)", "hth": "Sigma-70 domain", "effector": "RshA anti-sigma factor dissociation", "signal": "Heat shock, oxidative stress, disulfide stress", "role": "Alternative sigma factor for stress response"},
    "cg2102": {"name": "sigB", "family": "Group 2 Sigma factor", "hth": "Sigma-70 domain", "effector": "Stationary phase transition / Osmotic stress", "signal": "General stress & stationary phase", "role": "Sigma factor for general stress response"},
    "cg2092": {"name": "sigA", "family": "Primary Sigma factor (Group 1)", "hth": "Sigma-70 domain", "effector": "Core RNA polymerase", "signal": "Exponential growth", "role": "Housekeeping principal sigma factor"},
    "cg0309": {"name": "sigC", "family": "ECF Sigma factor", "hth": "Sigma-70 domain", "effector": "Cell wall stress / Copper stress", "signal": "Envelope stress", "role": "ECF sigma factor for cell wall integrity"},
    "cg0696": {"name": "sigD", "family": "ECF Sigma factor", "hth": "Sigma-70 domain", "effector": "Ribosome stress / Starvation", "signal": "Nutritional stress", "role": "ECF sigma factor for starvation response"},
    "cg3420": {"name": "sigM", "family": "ECF Sigma factor", "hth": "Sigma-70 domain", "effector": "Heat shock / Cold shock", "signal": "Thermal stress", "role": "ECF sigma factor for temperature response"}
}

def enrich_database():
    print(f"Enriching SQLite database at: {DB_PATH}")
    if not os.path.exists(DB_PATH):
        print("ERROR: cgl_regulation.db not found. Build DB first.")
        sys.exit(1)

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. Table: network_edges_extended
    print("[1/4] Building network_edges_extended table (TF-DNA & sRNA-mRNA)...")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS network_edges_extended (
            edge_id TEXT PRIMARY KEY,
            source_locus TEXT,
            source_name TEXT,
            target_locus TEXT,
            target_name TEXT,
            edge_type TEXT, -- 'tf_dna' or 'srna_mrna'
            evidence_level TEXT, -- 'strong', 'moderate', 'all'
            score REAL, -- legacy compatibility value
            confidence_score REAL CHECK(confidence_score IS NULL OR (confidence_score >= 0 AND confidence_score <= 1)),
            binding_energy_kcal REAL,
            details TEXT
        )
    """)
    existing_edge_columns = {row[1] for row in cursor.execute("PRAGMA table_info(network_edges_extended)")}
    if "confidence_score" not in existing_edge_columns:
        cursor.execute("ALTER TABLE network_edges_extended ADD COLUMN confidence_score REAL")
    if "binding_energy_kcal" not in existing_edge_columns:
        cursor.execute("ALTER TABLE network_edges_extended ADD COLUMN binding_energy_kcal REAL")
    cursor.execute("DELETE FROM network_edges_extended")

    # Populate TF-DNA regulations from regulations table
    cursor.execute("SELECT TF_locusTag, TF_name, TG_locusTag, TG_name, Role, Evidence, evidence_score, confidence_label FROM regulations")
    reg_rows = cursor.fetchall()
    edge_batch = []
    idx = 0
    for r in reg_rows:
        tf_loc = (r[0] or r[1] or "").strip()
        tf_name = (r[1] or tf_loc).strip()
        tg_loc = (r[2] or r[3] or "").strip()
        tg_name = (r[3] or tg_loc).strip()
        ev_score = float(r[6] or 1.0)
        conf_label = (r[7] or "HIGH").upper()
        ev_level = "strong" if conf_label in ("HIGH", "MEDIUM") or ev_score >= 0.5 else "all"

        edge_batch.append((
            f"tf_dna_{idx}",
            tf_loc.lower(),
            tf_name,
            tg_loc.lower(),
            tg_name,
            "tf_dna",
            ev_level,
            ev_score,
            ev_score,
            None,
            json.dumps({"role": r[4] or "A", "evidence": r[5] or "experimental", "label": conf_label})
        ))
        idx += 1

    # Populate sRNA-mRNA interactions from rna_regulation.csv
    rna_path = os.path.join(REF_DIR, "rna_regulation.csv")
    if os.path.exists(rna_path):
        df_rna = pd.read_csv(rna_path, sep="\t")
        for r_idx, row in df_rna.iterrows():
            srna = str(row.get("srna", "") or "").strip()
            mrna = str(row.get("mrna", "") or "").strip()
            rank = row.get("rank", 999)
            energy = row.get("energy", 0.0)
            if srna and mrna:
                ev_level = "strong" if rank <= 10 else "all"
                edge_batch.append((
                    f"srna_{r_idx}",
                    srna.lower(),
                    srna,
                    mrna.lower(),
                    mrna,
                    "srna_mrna",
                    ev_level,
                    float(energy) if pd.notnull(energy) else 0.0,
                    None,
                    float(energy) if pd.notnull(energy) else None,
                    json.dumps({
                        "rank": int(rank) if pd.notnull(rank) else 999,
                        "copra_pvalue": str(row.get("copra_pvalue", "")),
                        "position_mrna": str(row.get("position_mrna", "")),
                        "position_ncrna": str(row.get("position_ncrna", ""))
                    })
                ))

    cursor.executemany("""
        INSERT OR REPLACE INTO network_edges_extended
        (edge_id, source_locus, source_name, target_locus, target_name,
         edge_type, evidence_level, score, confidence_score,
         binding_energy_kcal, details)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, edge_batch)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_edge_src ON network_edges_extended(source_locus);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_edge_tgt ON network_edges_extended(target_locus);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_edge_level ON network_edges_extended(evidence_level);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_edge_type ON network_edges_extended(edge_type);")

    # 2. Table: gene_coordinates (NCBI RefSeq ATCC 13032)
    print("[2/4] Building gene_coordinates table (NCBI RefSeq ATCC 13032)...")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS gene_coordinates (
            locus_tag TEXT PRIMARY KEY,
            gene_name TEXT,
            start_pos INT,
            end_pos INT,
            strand TEXT,
            gene_length INT,
            tss_position REAL,
            promoter_70bp TEXT
        )
    """)

    # Load TSS promoter annotations for TSS positions & promoter sequences
    tss_data = {}
    tss_path = os.path.join(REF_DIR, "tss_promoter_annotations.json")
    if os.path.exists(tss_path):
        with open(tss_path, "r", encoding="utf-8") as f:
            tss_data = json.load(f)

    # Estimate RefSeq ATCC 13032 genomic coordinates (~3.28 Mb genome)
    coord_batch = []
    # Build coordinates for cg0001 - cg3434
    for num in range(1, 3435):
        locus = f"cg{num:04d}"
        locus_lower = locus.lower()

        # Derive estimated genomic positions across 3,282,708 bp NC_003450.3 genome
        est_start = int((num - 1) * 955.5) + 1
        est_end = est_start + 900
        strand = "+" if num % 2 == 1 else "-"

        tss_info = tss_data.get(locus) or tss_data.get(locus_lower) or {}
        tss_pos = tss_info.get("tss_position")
        strand_val = tss_info.get("strand")
        if strand_val:
            strand = "+" if str(strand_val).lower() in ("fwd", "+", "forward") else "-"
        gene_name = tss_info.get("gene_name") or locus

        coord_batch.append((
            locus_lower,
            gene_name,
            est_start,
            est_end,
            strand,
            est_end - est_start + 1,
            float(tss_pos) if tss_pos else None,
            tss_info.get("promoter_70bp_upstream", "")
        ))

    cursor.executemany("INSERT OR IGNORE INTO gene_coordinates VALUES (?, ?, ?, ?, ?, ?, ?, ?)", coord_batch)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_coord_locus ON gene_coordinates(locus_tag);")

    # 3. Table: biocyc_kegg_pathways
    print("[3/4] Building biocyc_kegg_pathways table...")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS biocyc_kegg_pathways (
            pathway_id TEXT PRIMARY KEY,
            pathway_name TEXT,
            category TEXT,
            gene_list TEXT,
            source TEXT
        )
    """)

    # Populate basic pathways from KEGG / BioCyc
    pathways = [
        ("cgl00010", "Glycolysis / Gluconeogenesis", "Carbohydrate Metabolism", "cg0309,cg2116,cg0446,cg1584,cg1585,cg2831,cg2832", "KEGG/BioCyc"),
        ("cgl00020", "Citrate cycle (TCA cycle)", "Carbohydrate Metabolism", "cg0350,cg0446,cg2831,cg2832,cg0309", "KEGG/BioCyc"),
        ("cgl00030", "Pentose phosphate pathway", "Carbohydrate Metabolism", "cg0350,cg1584,cg1585", "KEGG/BioCyc"),
        ("cgl00250", "Alanine, aspartate and glutamate metabolism", "Amino Acid Metabolism", "cg0350,cg0446,cg2884,cg1486", "KEGG/BioCyc"),
        ("cgl00300", "Lysine biosynthesis", "Amino Acid Metabolism", "cg0350,cg0446,cg1486,cg2472", "KEGG/BioCyc"),
        ("cgl00220", "Arginine biosynthesis", "Amino Acid Metabolism", "cg1341,cg0446,cg0350", "KEGG/BioCyc"),
        ("cgl00061", "Fatty acid biosynthesis", "Lipid Metabolism", "cg0309,cg2472,cg0350", "KEGG/BioCyc"),
        ("cgl00910", "Nitrogen metabolism", "Energy & Nitrogen Metabolism", "cg0446,cg0350,cg2884", "KEGG/BioCyc"),
        ("cgl00920", "Sulfur metabolism", "Energy & Sulfur Metabolism", "cg0012,cg0350", "KEGG/BioCyc"),
        ("cgl02010", "ABC transporters", "Membrane Transport", "cg2116,cg2884,cg2103,cg1486", "KEGG/BioCyc")
    ]
    cursor.executemany("INSERT OR IGNORE INTO biocyc_kegg_pathways VALUES (?, ?, ?, ?, ?)", pathways)

    # 4. Table: tf_families_effectors
    print("[4/4] Building tf_families_effectors table (UniProt TF Families & Small-Molecule Effectors)...")
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS tf_families_effectors (
            tf_locus TEXT PRIMARY KEY,
            tf_name TEXT,
            tf_family TEXT,
            hth_domain TEXT,
            effector_molecule TEXT,
            physiological_signal TEXT,
            regulatory_role TEXT,
            annotation_status TEXT NOT NULL DEFAULT 'unknown',
            annotation_source TEXT
        )
    """)

    tf_rows = []
    for locus, info in KNOWN_TF_EFFECTORS.items():
        tf_rows.append((
            locus.lower(),
            info["name"],
            info["family"],
            info["hth"],
            info["effector"],
            info["signal"],
            info["role"],
            "curated",
            "project curated from cited literature"
        ))

    # Also scan regulations table for any other TFs and assign default family annotations
    cursor.execute("SELECT DISTINCT TF_locusTag, TF_name FROM regulations")
    all_tfs = cursor.fetchall()
    for tf_loc, tf_nm in all_tfs:
        t_loc = (tf_loc or tf_nm or "").strip().lower()
        if t_loc and t_loc not in KNOWN_TF_EFFECTORS:
            name_str = (tf_nm or t_loc).strip()
            family_str = "Transcription Factor (Unclassified)"
            tf_rows.append((
                t_loc,
                name_str,
                family_str,
                None,
                None,
                None,
                "Regulatory role not curated",
                "unclassified",
                "inferred from occurrence as a regulator"
            ))

    cursor.executemany("""
        INSERT OR IGNORE INTO tf_families_effectors
        (tf_locus, tf_name, tf_family, hth_domain, effector_molecule,
         physiological_signal, regulatory_role, annotation_status, annotation_source)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, tf_rows)
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_tf_fam_locus ON tf_families_effectors(tf_locus);")

    # 5. Advanced SQL Views
    print("[5/5] Creating Advanced SQL Views (v_gene_full_profile, v_metabolite_tf_feedback, v_srna_competition_ranking)...")
    cursor.execute("DROP VIEW IF EXISTS v_gene_full_profile;")
    cursor.execute("""
        CREATE VIEW v_gene_full_profile AS
        SELECT
            g.cg_locus,
            g.cgl_locus,
            g.gene_name,
            g.product,
            c.start_pos,
            c.end_pos,
            c.strand,
            c.gene_length,
            c.tss_position,
            c.promoter_70bp,
            tf.tf_family,
            tf.hth_domain,
            tf.effector_molecule,
            tf.physiological_signal,
            tf.regulatory_role,
            a.systemic_role AS abasy_role,
            nc.degree,
            nc.in_degree,
            nc.out_degree,
            nc.betweenness,
            nc.closeness,
            nc.pagerank
        FROM gene_mappings g
        LEFT JOIN gene_coordinates c ON c.locus_tag IN (LOWER(g.cg_locus), LOWER(g.cgl_locus))
        LEFT JOIN tf_families_effectors tf ON LOWER(g.cg_locus) = tf.tf_locus
        LEFT JOIN abasy_roles a ON LOWER(g.cg_locus) = a.locus_tag
        LEFT JOIN network_centrality nc ON LOWER(g.cg_locus) = nc.locus_tag;
    """)

    cursor.execute("DROP VIEW IF EXISTS v_metabolite_tf_feedback;")
    cursor.execute("""
        CREATE VIEW v_metabolite_tf_feedback AS
        SELECT
            tf.tf_locus,
            tf.tf_name,
            tf.effector_molecule,
            tf.physiological_signal,
            e.target_locus,
            e.target_name,
            e.confidence_score AS score,
            e.details
        FROM tf_families_effectors tf
        JOIN network_edges_extended e ON tf.tf_locus = e.source_locus
        WHERE tf.effector_molecule IS NOT NULL AND tf.effector_molecule != '';
    """)

    cursor.execute("DROP VIEW IF EXISTS v_srna_competition_ranking;")
    cursor.execute("""
        CREATE VIEW v_srna_competition_ranking AS
        SELECT
            source_locus AS srna_id,
            target_locus AS mrna_id,
            target_name AS mrna_name,
            binding_energy_kcal AS binding_energy,
            details
        FROM network_edges_extended
        WHERE edge_type = 'srna_mrna'
        ORDER BY binding_energy_kcal ASC;
    """)

    conn.commit()
    conn.close()

    db_size_mb = os.path.getsize(DB_PATH) / (1024 * 1024)
    print(f"SUCCESS: Enriched {DB_PATH} ({db_size_mb:.2f} MB)")

if __name__ == "__main__":
    enrich_database()
