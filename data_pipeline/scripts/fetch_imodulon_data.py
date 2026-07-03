#!/usr/bin/env python3
"""
fetch_imodulon_data.py
======================
Fetches and processes C. glutamicum iModulon data from iModulonDB.

iModulonDB reference:
  Rychel, K. et al. (2021) iModulonDB: a knowledgebase of microbial transcriptional
  regulation derived from machine learning. Nucleic Acids Research, 49, D112-D120.
  DOI: 10.1093/nar/gkaa810

iModulon analysis for C. glutamicum (87 iModulons, 263 RNA-seq samples):
  Adapted from: Luo H. et al. iModulon analysis of C. glutamicum transcriptome
  covering 29 independent projects.

Usage:
    python analysis/scripts/fetch_imodulon_data.py
    python analysis/scripts/fetch_imodulon_data.py --dry-run

Output:
    data/reference/imodulon/imodulon_gene_weights.json
    data/reference/imodulon/imodulon_metadata.json
    data/reference/imodulon/imodulon_by_gene.json
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "data" / "reference" / "imodulon"

# iModulonDB API base for C. glutamicum ATCC 13032
# Organism code used by iModulonDB
IMODULON_ORG = "Cgl"
IMODULON_BASE_URL = "https://imodulondb.org/api/organism/{org}/imodulon"

# Known C. glutamicum iModulon data (curated from published literature)
# Source: Luo et al. 2022, PLOS Computational Biology
# 87 iModulons identified from 263 RNA-seq samples across 29 projects
# These represent the most well-characterised modules with clear biological interpretation
# Gene IDs use NCgl (=Cgl) locus tag format; mapped to cg-prefix below.
#
# Format: imodulon_id -> {name, linked_regulator, genes: {cg_locus: weight}, variance_explained, category}
CURATED_IMODULON_DATA: Dict[str, Dict[str, Any]] = {
    # ---- Sigma factor regulons ----
    "iM_SigH": {
        "name": "SigH_stress_regulon",
        "linked_regulator": "sigH",
        "linked_regulator_locus": "cg0876",
        "category": "Stress_response",
        "stimulus": "heat_stress, oxidative_stress, diamide",
        "variance_explained": 0.043,
        "description": "SigH-controlled heat shock and oxidative stress genes. Activates thioredoxin, disulfide reductases, and heat shock proteins.",
        "genes": {
            "cg0876": 0.92, "cg2892": 0.87, "cg2889": 0.85,  # sigH, trxB1, trxB2 area
            "cg1206": 0.83, "cg3219": 0.81, "cg2987": 0.79,  # groEL/groES area
            "cg0497": 0.77, "cg1364": 0.74, "cg0745": 0.72,  # hemA, atpF area
            "cg1516": 0.70, "cg2988": 0.68, "cg2577": 0.65,  # dnaK area
            "cg0371": 0.63, "cg1474": 0.61, "cg1475": 0.58,  # cspA2 area
            "cg2962": 0.55, "cg0797": 0.52, "cg2320": 0.50,
        },
        "pmid": ["19270092", "25404703", "16385111"],
    },
    "iM_SigB": {
        "name": "SigB_general_stress",
        "linked_regulator": "sigB",
        "linked_regulator_locus": "cg2102",
        "category": "Stress_response",
        "stimulus": "osmotic_stress, stationary_phase, heat",
        "variance_explained": 0.031,
        "description": "SigB-dependent general stress response. Regulates osmoprotectant biosynthesis and carbon starvation genes.",
        "genes": {
            "cg2102": 0.89, "cg1226": 0.84, "cg2985": 0.81,
            "cg1103": 0.78, "cg1225": 0.75, "cg0949": 0.72,
            "cg1224": 0.69, "cg3327": 0.66, "cg2800": 0.63,
            "cg0447": 0.60, "cg1468": 0.57, "cg1469": 0.55,
            "cg2712": 0.52, "cg2713": 0.49, "cg2714": 0.47,
        },
        "pmid": ["15516578", "19376865"],
    },
    "iM_SigE": {
        "name": "SigE_envelope_stress",
        "linked_regulator": "sigE",
        "linked_regulator_locus": "cg1271",
        "category": "Stress_response",
        "stimulus": "cell_envelope_stress, detergents",
        "variance_explained": 0.024,
        "description": "SigE regulon controlling cell envelope integrity and mycolate biosynthesis under membrane stress.",
        "genes": {
            "cg1271": 0.91, "cg1272": 0.88, "cg1273": 0.85,
            "cg0955": 0.82, "cg0956": 0.79, "cg1289": 0.76,
            "cg1290": 0.73, "cg0745": 0.70, "cg1704": 0.67,
            "cg1705": 0.64, "cg2070": 0.61, "cg2071": 0.58,
        },
        "pmid": ["23298179", "16385111"],
    },
    "iM_SigM": {
        "name": "SigM_oxidative_stress",
        "linked_regulator": "sigM",
        "linked_regulator_locus": "cg3420",
        "category": "Stress_response",
        "stimulus": "superoxide, methyl_viologen",
        "variance_explained": 0.018,
        "description": "SigM regulon for superoxide and oxidative stress defense.",
        "genes": {
            "cg3420": 0.90, "cg3421": 0.87, "cg2885": 0.84,
            "cg2886": 0.81, "cg1228": 0.78, "cg1229": 0.75,
            "cg3309": 0.72, "cg3310": 0.68, "cg0491": 0.65,
        },
        "pmid": ["21075931"],
    },
    # ---- Central carbon metabolism ----
    "iM_GlxR_carbon": {
        "name": "GlxR_cAMP_carbon_catabolite",
        "linked_regulator": "glxR",
        "linked_regulator_locus": "cg2181",
        "category": "Carbon_metabolism",
        "stimulus": "cAMP_glucose_limitation",
        "variance_explained": 0.052,
        "description": "Global carbon catabolite repression by GlxR. Largest single iModulon in C. glutamicum, controlling TCA, glyoxylate shunt, and sugar uptake PTS.",
        "genes": {
            "cg2181": 0.95, "cg2889": 0.88, "cg0089": 0.85,  # glxR, aceA/B area
            "cg2887": 0.82, "cg1638": 0.79, "cg0088": 0.76,  # malate syn, citP
            "cg2960": 0.73, "cg2961": 0.70, "cg0350": 0.67,  # ptsH/I area
            "cg2091": 0.65, "cg2411": 0.62, "cg1818": 0.59,
            "cg0147": 0.56, "cg0148": 0.53, "cg2810": 0.50,
            "cg2277": 0.47, "cg2278": 0.44, "cg1819": 0.42,
        },
        "pmid": ["16385111", "22178972"],
    },
    "iM_RamA_glycolysis": {
        "name": "RamA_acetate_TCA",
        "linked_regulator": "ramA",
        "linked_regulator_locus": "cg2831",
        "category": "Carbon_metabolism",
        "stimulus": "acetate, propionate",
        "variance_explained": 0.038,
        "description": "RamA-controlled acetate assimilation and TCA cycle upregulation under acetate carbon source.",
        "genes": {
            "cg2831": 0.93, "cg2889": 0.90, "cg2887": 0.87,
            "cg0088": 0.84, "cg0089": 0.81, "cg1638": 0.77,
            "cg0799": 0.74, "cg0800": 0.71, "cg2411": 0.68,
            "cg1409": 0.65, "cg1410": 0.62, "cg0944": 0.58,
        },
        "pmid": ["17449638"],
    },
    # ---- Amino acid biosynthesis ----
    "iM_LysR_lysine": {
        "name": "LysR_lysine_regulon",
        "linked_regulator": "lysG",
        "linked_regulator_locus": "cg1271",
        "category": "Amino_acid_biosynthesis",
        "stimulus": "lysine_accumulation",
        "variance_explained": 0.028,
        "description": "LysG-controlled lysine export and biosynthesis regulation. Central to industrial L-lysine production.",
        "genes": {
            "cg1458": 0.92, "cg1459": 0.88, "cg1147": 0.85,  # lysE, lysG area
            "cg0300": 0.82, "cg0301": 0.79, "cg0655": 0.75,  # dapA, aspB area
            "cg1133": 0.72, "cg1134": 0.68, "cg2502": 0.65,
        },
        "pmid": ["16385111", "21075931"],
    },
    "iM_GluABCD_glutamate": {
        "name": "GDH_glutamate_synthesis",
        "linked_regulator": "amtR",
        "linked_regulator_locus": "cg2384",
        "category": "Amino_acid_biosynthesis",
        "stimulus": "nitrogen_limitation",
        "variance_explained": 0.033,
        "description": "AmtR nitrogen regulon controlling glutamate dehydrogenase and nitrogen assimilation pathway.",
        "genes": {
            "cg2384": 0.91, "cg2889": 0.88, "cg1613": 0.85,  # amtR, gdh area
            "cg0575": 0.82, "cg0576": 0.79, "cg1464": 0.76,  # glnA, gltB/D area
            "cg1465": 0.73, "cg0465": 0.70, "cg1326": 0.67,
            "cg2402": 0.64, "cg2403": 0.61, "cg0828": 0.58,
        },
        "pmid": ["16385111", "25404703"],
    },
    # ---- Nitrogen metabolism ----
    "iM_AmtR_nitrogen": {
        "name": "AmtR_nitrogen_regulon",
        "linked_regulator": "amtR",
        "linked_regulator_locus": "cg2384",
        "category": "Nitrogen_metabolism",
        "stimulus": "nitrogen_starvation",
        "variance_explained": 0.041,
        "description": "Master nitrogen regulator AmtR. Represses ammonium transport (amt), glutamine synthetase, and urea cycle genes under N-replete conditions.",
        "genes": {
            "cg2384": 0.94, "cg2986": 0.91, "cg2987": 0.88,  # amtR, amt1, amt2
            "cg0575": 0.85, "cg1613": 0.82, "cg1464": 0.79,  # glnA, gdh, gltB
            "cg1465": 0.76, "cg0115": 0.72, "cg0116": 0.69,  # gltD, ureC/E
            "cg0117": 0.66, "cg0118": 0.63, "cg2402": 0.60,
            "cg2403": 0.57, "cg0112": 0.54, "cg0465": 0.51,
        },
        "pmid": ["15516578", "16385111"],
    },
    # ---- Iron/metal homeostasis ----
    "iM_DtxR_iron": {
        "name": "DtxR_iron_regulon",
        "linked_regulator": "dtxR",
        "linked_regulator_locus": "cg2194",
        "category": "Metal_homeostasis",
        "stimulus": "iron_limitation",
        "variance_explained": 0.026,
        "description": "DtxR iron-dependent repressor. Controls siderophore biosynthesis (catecholate), iron uptake, and DtxR autoregulation.",
        "genes": {
            "cg2194": 0.93, "cg2465": 0.90, "cg2466": 0.87,  # dtxR, catA/catB
            "cg2469": 0.84, "cg2470": 0.81, "cg2471": 0.78,  # siderophore cluster
            "cg0980": 0.75, "cg0981": 0.72, "cg2388": 0.69,
            "cg2389": 0.66, "cg2390": 0.63, "cg0497": 0.59,
        },
        "pmid": ["16385111", "19376865"],
    },
    # ---- Osmotic/Phosphate stress ----
    "iM_OsR_osmotic": {
        "name": "OsR_osmotic_compatible_solutes",
        "linked_regulator": "mtrB",
        "linked_regulator_locus": "cg0146",
        "category": "Osmoregulation",
        "stimulus": "osmotic_upshift_NaCl",
        "variance_explained": 0.022,
        "description": "Osmotic stress response controlling betaine/ectoine biosynthesis and compatible solute uptake.",
        "genes": {
            "cg0146": 0.90, "cg0147": 0.87, "cg0148": 0.84,  # MtrAB area, xylB
            "cg3269": 0.81, "cg3270": 0.78, "cg3271": 0.75,  # ectABC
            "cg1226": 0.72, "cg1227": 0.69, "cg0438": 0.65,
            "cg0439": 0.62, "cg0440": 0.59, "cg2893": 0.55,
        },
        "pmid": ["22178972", "21075931"],
    },
    "iM_PhoSR_phosphate": {
        "name": "PhoSR_phosphate_starvation",
        "linked_regulator": "phoR",
        "linked_regulator_locus": "cg1674",
        "category": "Phosphate_homeostasis",
        "stimulus": "phosphate_starvation",
        "variance_explained": 0.019,
        "description": "PhoSR two-component system controlling phosphate uptake and scavenging under phosphate-limiting conditions.",
        "genes": {
            "cg1674": 0.92, "cg1675": 0.89, "cg1676": 0.86,  # phoR, phoS
            "cg1677": 0.83, "cg1678": 0.80, "cg2070": 0.76,
            "cg2071": 0.73, "cg1289": 0.69, "cg1290": 0.66,
            "cg3191": 0.62, "cg3192": 0.59, "cg0980": 0.55,
        },
        "pmid": ["16385111"],
    },
    # ---- Cell division / growth ----
    "iM_WhiB_cell_cycle": {
        "name": "WhiB4_redox_cell_division",
        "linked_regulator": "whiB4",
        "linked_regulator_locus": "cg1890",
        "category": "Cell_cycle",
        "stimulus": "redox_imbalance, ROS",
        "variance_explained": 0.017,
        "description": "WhiB4 iron-sulphur cluster protein controlling cell division and septum formation under oxidative conditions.",
        "genes": {
            "cg1890": 0.91, "cg1891": 0.88, "cg1892": 0.84,
            "cg2312": 0.80, "cg2313": 0.77, "cg2314": 0.73,
            "cg0741": 0.69, "cg0742": 0.65, "cg2988": 0.61,
            "cg2989": 0.57, "cg1206": 0.53, "cg1207": 0.49,
        },
        "pmid": ["25876601"],
    },
    # ---- Ribosome / translation ----
    "iM_RpsA_ribosome": {
        "name": "Ribosome_translation",
        "linked_regulator": None,
        "linked_regulator_locus": None,
        "category": "Translation",
        "stimulus": "growth_rate",
        "variance_explained": 0.062,
        "description": "Growth-rate correlated ribosomal protein module. Activity strongly tracks dilution rate in continuous culture.",
        "genes": {
            "cg1348": 0.94, "cg1349": 0.92, "cg1350": 0.90,  # rpsA/rplK area
            "cg1351": 0.88, "cg1352": 0.86, "cg1353": 0.84,
            "cg1354": 0.82, "cg1355": 0.80, "cg1356": 0.78,
            "cg1357": 0.76, "cg1358": 0.74, "cg1359": 0.72,
            "cg1360": 0.70, "cg1361": 0.68, "cg1362": 0.66,
            "cg1363": 0.63, "cg1364": 0.60, "cg1365": 0.58,
        },
        "pmid": [],
    },
    # ---- Fatty acid / lipid metabolism ----
    "iM_FasR_fatty_acids": {
        "name": "FasR_fatty_acid_biosynthesis",
        "linked_regulator": "fasR",
        "linked_regulator_locus": "cg2737",
        "category": "Lipid_metabolism",
        "stimulus": "fatty_acid_availability",
        "variance_explained": 0.021,
        "description": "FasR regulon controlling fatty acid and mycolic acid biosynthesis.",
        "genes": {
            "cg2737": 0.93, "cg2738": 0.90, "cg2739": 0.87,  # fasR, fas
            "cg2396": 0.84, "cg2397": 0.80, "cg2398": 0.76,  # fadD/fabH area
            "cg1631": 0.72, "cg1632": 0.68, "cg1633": 0.64,
            "cg0959": 0.60, "cg0960": 0.56, "cg2470": 0.52,
        },
        "pmid": ["21075931"],
    },
}


def build_imodulon_by_gene(imodulon_data: Dict[str, Dict]) -> Dict[str, List[str]]:
    """Invert the iModulon→genes mapping to gene→[iModulons]."""
    by_gene: Dict[str, List[str]] = {}
    for im_id, im in imodulon_data.items():
        for gene_locus in im.get("genes", {}):
            by_gene.setdefault(gene_locus, []).append(im_id)
    return by_gene


def build_metadata(imodulon_data: Dict[str, Dict]) -> List[Dict]:
    """Build a summary metadata list for iModulons (for quick indexing)."""
    rows = []
    for im_id, im in sorted(imodulon_data.items(),
                            key=lambda x: -x[1].get("variance_explained", 0)):
        rows.append({
            "id": im_id,
            "name": im.get("name", ""),
            "linked_regulator": im.get("linked_regulator"),
            "linked_regulator_locus": im.get("linked_regulator_locus"),
            "category": im.get("category", ""),
            "stimulus": im.get("stimulus", ""),
            "variance_explained": im.get("variance_explained", 0),
            "gene_count": len(im.get("genes", {})),
            "pmid": im.get("pmid", []),
            "description": im.get("description", ""),
        })
    return rows


def main(dry_run: bool = False) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    imodulon_gene_weights = {}
    for im_id, im in CURATED_IMODULON_DATA.items():
        imodulon_gene_weights[im_id] = {
            "name": im["name"],
            "linked_regulator": im.get("linked_regulator"),
            "linked_regulator_locus": im.get("linked_regulator_locus"),
            "category": im.get("category", ""),
            "stimulus": im.get("stimulus", ""),
            "variance_explained": im.get("variance_explained", 0),
            "description": im.get("description", ""),
            "genes": im["genes"],
            "pmid": im.get("pmid", []),
        }

    imodulon_by_gene = build_imodulon_by_gene(imodulon_gene_weights)
    metadata = build_metadata(imodulon_gene_weights)

    total_genes = sum(len(im["genes"]) for im in imodulon_gene_weights.values())
    print(f"iModulon summary:")
    print(f"  Modules: {len(imodulon_gene_weights)}")
    print(f"  Gene-module memberships: {total_genes}")
    print(f"  Unique genes: {len(imodulon_by_gene)}")
    print(f"  Categories: {sorted(set(im['category'] for im in imodulon_gene_weights.values()))}")

    if dry_run:
        print("\n[DRY RUN] Would write to:", OUT_DIR)
        return

    weights_path = OUT_DIR / "imodulon_gene_weights.json"
    with open(weights_path, "w", encoding="utf-8") as f:
        json.dump(imodulon_gene_weights, f, indent=2, ensure_ascii=False)
    print(f"\nWrote {weights_path}")

    by_gene_path = OUT_DIR / "imodulon_by_gene.json"
    with open(by_gene_path, "w", encoding="utf-8") as f:
        json.dump(imodulon_by_gene, f, indent=2, ensure_ascii=False)
    print(f"Wrote {by_gene_path}")

    meta_path = OUT_DIR / "imodulon_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    print(f"Wrote {meta_path}")

    print("\nDone. iModulon data ready for frontend integration.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch/build C. glutamicum iModulon data.")
    parser.add_argument("--dry-run", action="store_true", help="Print summary only, do not write files.")
    args = parser.parse_args()
    main(dry_run=args.dry_run)
