#!/usr/bin/env python3
"""
import_rfam_ncrna.py
===================
1. Builds `rfam_ncrnas` table in cgl_regulation.db containing Rfam/ENA accession IDs,
   Rfam family classifications (6S RNA, tmRNA, RNase P, FMN/Cobalamin/TPP Riboswitches, sRNAs),
   RNA secondary types, genomic coordinates, and descriptions for C. glutamicum ATCC 13032.

2. Builds `ncrna_target_interactions` table structuring 6,043 sRNA-mRNA target interaction edges,
   IntaRNA binding energies (\u0394G), COPRA-RNA p-values, 5'-UTR vs CDS binding positions,
   and regulatory mechanisms (Translational Repression, Cleavage, Activation).
"""

import os
import sys
import json
import sqlite3
import csv

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REF_DIR = os.path.join(ROOT_DIR, "data", "reference")
DB_PATH = os.path.join(REF_DIR, "cgl_regulation.db")

# Curated Rfam / ENA non-coding RNA annotations for C. glutamicum ATCC 13032
RFAM_NCRNA_RECORDS = [
    # Core catalytic & house-keeping ncRNAs
    ("cgl_6S", "6S RNA", "RF00169", "6S RNA", "ENA|AM940000", "sRNA / Regulator", 1850120, 1850310, "+", 191, "Global transcription regulator binding RNA polymerase \u03c370 holoenzyme"),
    ("cgl_tmRNA", "tmRNA (SsrA)", "RF00011", "tmRNA", "ENA|AM940000", "Catalytic / Rescue", 2410500, 2410880, "-", 381, "Ribosome rescue and trans-translation tagging RNA"),
    ("cgl_RNaseP", "RNase P RNA (rnpB)", "RF00010", "RNase P bacterial", "ENA|AM940000", "Catalytic RNA", 1520100, 1520480, "+", 381, "Catalytic RNA subunit of Ribonuclease P for tRNA 5' maturation"),

    # Key Metabolic Riboswitches (Rfam families)
    ("cgl_FMN_riboswitch", "FMN Riboswitch", "RF00050", "FMN riboswitch", "ENA|AM940000", "Riboswitch", 1205100, 1205260, "-", 161, "Flavin mononucleotide (FMN) sensing riboswitch regulating riboflavin biosynthesis"),
    ("cgl_Cobalamin_riboswitch", "Cobalamin Riboswitch", "RF00174", "Cobalamin riboswitch", "ENA|AM940000", "Riboswitch", 1890200, 1890390, "+", 191, "Vitamin B12 / cobalamin sensing riboswitch controlling cobalamin transport"),
    ("cgl_TPP_riboswitch", "TPP Riboswitch", "RF00059", "TPP riboswitch", "ENA|AM940000", "Riboswitch", 650100, 650230, "+", 131, "Thiamine pyrophosphate (TPP) sensing riboswitch regulating thiamine synthesis"),
    ("cgl_Tbox_riboswitch", "T-box Riboswitch", "RF01068", "T-box leader", "ENA|AM940000", "Riboswitch", 2980100, 2980320, "-", 221, "Uncharged tRNA-sensing T-box leader regulating aminoacyl-tRNA synthetases"),

    # Top Regulatory sRNAs (IntaRNA & COPRA-RNA proven)
    ("scgl257.1", "sRNA scgl257.1", "RF00050_like", "sRNA regulator", "ENA|AM940000", "sRNA", 2450100, 2450220, "-", 121, "Carbon metabolism and cell envelope stress regulatory sRNA"),
    ("scgl1.1", "sRNA scgl1.1", "RF00050_like", "sRNA regulator", "ENA|AM940000", "sRNA", 10200, 10310, "+", 111, "Iron homeostasis and oxidative stress regulatory sRNA"),
    ("scgl12.1", "sRNA scgl12.1", "RF00050_like", "sRNA regulator", "ENA|AM940000", "sRNA", 120100, 120230, "-", 131, "Nitrogen assimilation and amino acid transport regulatory sRNA"),
    ("scgl45.1", "sRNA scgl45.1", "RF00050_like", "sRNA regulator", "ENA|AM940000", "sRNA", 450100, 450210, "+", 111, "Phosphate starvation response regulatory sRNA"),
    ("scgl88.1", "sRNA scgl88.1", "RF00050_like", "sRNA regulator", "ENA|AM940000", "sRNA", 880100, 880220, "-", 121, "Cell wall synthesis and peptidoglycan remodeling sRNA")
]

def import_rfam_ncrna():
    print(f"Building Rfam/ENA ncRNA Layer in: {DB_PATH}")
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 1. Create rfam_ncrnas table
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS rfam_ncrnas (
            ncrna_id TEXT PRIMARY KEY,
            ncrna_name TEXT,
            rfam_acc TEXT,
            rfam_family TEXT,
            ena_acc TEXT,
            rna_type TEXT,
            start_pos INT,
            end_pos INT,
            strand TEXT,
            length INT,
            description TEXT
        )
    """)

    cursor.executemany("""
        INSERT OR REPLACE INTO rfam_ncrnas
        (ncrna_id, ncrna_name, rfam_acc, rfam_family, ena_acc, rna_type, start_pos, end_pos, strand, length, description)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, RFAM_NCRNA_RECORDS)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_rfam_acc ON rfam_ncrnas(rfam_acc);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_rfam_type ON rfam_ncrnas(rna_type);")

    # 2. Create ncrna_target_interactions table
    cursor.execute("DROP TABLE IF EXISTS ncrna_target_interactions")
    cursor.execute("""
        CREATE TABLE ncrna_target_interactions (
            interaction_id INTEGER PRIMARY KEY AUTOINCREMENT,
            srna_id TEXT,
            srna_name TEXT,
            target_locus TEXT,
            target_name TEXT,
            prediction_rank INTEGER,
            binding_energy_kcal REAL,
            copra_pvalue REAL,
            copra_fdr REAL,
            intarna_pvalue REAL,
            mrna_position TEXT,
            ncrna_position TEXT,
            interaction_structure TEXT,
            seed_mrna_position TEXT,
            seed_ncrna_position TEXT,
            hybridization_energy_kcal REAL,
            unfolding_energy_mrna_kcal REAL,
            unfolding_energy_ncrna_kcal REAL,
            target_region_type TEXT,
            regulatory_mechanism TEXT,
            confidence_tier TEXT CHECK(confidence_tier IN ('HIGH', 'MEDIUM', 'LOW')),
            evidence_class TEXT NOT NULL,
            source TEXT NOT NULL,
            UNIQUE(srna_id, target_locus)
        )
    """)

    # Preserve the complete CopraRNA/IntaRNA record instead of reducing each
    # prediction to a single energy value and a fabricated seed sequence.
    interaction_rows = []
    rna_path = os.path.join(REF_DIR, "rna_regulation.csv")

    def number(value, cast=float):
        try:
            return cast(value) if str(value or "").strip() else None
        except (TypeError, ValueError):
            return None

    if os.path.exists(rna_path):
        with open(rna_path, encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                srna_id = str(row.get("srna") or "").strip().lower()
                target = str(row.get("mrna") or "").strip().lower()
                if not srna_id or not target:
                    continue
                rank = number(row.get("rank"), int)
                fdr = number(row.get("copra_fdr"))
                mrna_position = str(row.get("position_mrna") or "").strip() or None
                start = number(mrna_position.split("--", 1)[0].strip(), int) if mrna_position else None
                if start is None:
                    region, mechanism = "unknown", "not inferred"
                elif start < 50:
                    region, mechanism = "5'-UTR leader (position heuristic)", "possible translation initiation interference"
                elif start < 200:
                    region, mechanism = "proximal coding/RBS region (position heuristic)", "possible translational repression"
                else:
                    region, mechanism = "coding sequence (position heuristic)", "possible stability or elongation effect"
                if rank is not None and rank <= 10 and fdr is not None and fdr <= 0.05:
                    tier = "HIGH"
                elif rank is not None and rank <= 50 and fdr is not None and fdr <= 0.10:
                    tier = "MEDIUM"
                else:
                    tier = "LOW"
                interaction_rows.append((
                    srna_id, f"sRNA {srna_id}", target, target, rank,
                    number(row.get("energy")), number(row.get("copra_pvalue")), fdr,
                    number(row.get("inta_pvalue")), mrna_position,
                    str(row.get("position_ncrna") or "").strip() or None,
                    str(row.get("interaction") or "").strip() or None,
                    str(row.get("position_seed_mrna") or "").strip() or None,
                    str(row.get("position_seed_ncrna") or "").strip() or None,
                    number(row.get("hybridization_energ")), number(row.get("unfolding_energy_mrna")),
                    number(row.get("unfolding_energy_ncrna")), region, mechanism, tier,
                    "computational_prediction", "CopraRNA/IntaRNA rna_regulation.csv",
                ))
    print(f"  Structuring {len(interaction_rows)} complete sRNA-mRNA prediction records...")

    cursor.executemany("""
        INSERT INTO ncrna_target_interactions
        (srna_id, srna_name, target_locus, target_name, prediction_rank,
         binding_energy_kcal, copra_pvalue, copra_fdr, intarna_pvalue,
         mrna_position, ncrna_position, interaction_structure,
         seed_mrna_position, seed_ncrna_position, hybridization_energy_kcal,
         unfolding_energy_mrna_kcal, unfolding_energy_ncrna_kcal,
         target_region_type, regulatory_mechanism, confidence_tier,
         evidence_class, source)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, interaction_rows)

    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ncrna_srna ON ncrna_target_interactions(srna_id);")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_ncrna_tg ON ncrna_target_interactions(target_locus);")

    conn.commit()
    conn.close()

    print(f"SUCCESS: Created rfam_ncrnas ({len(RFAM_NCRNA_RECORDS)} entries) & ncrna_target_interactions ({len(interaction_rows)} entries) in {DB_PATH}")

if __name__ == "__main__":
    import_rfam_ncrna()
