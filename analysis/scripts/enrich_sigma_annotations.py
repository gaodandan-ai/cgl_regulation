#!/usr/bin/env python3
"""
enrich_sigma_annotations.py
============================
Builds a curated JSON annotation file for the alternative sigma factors
in Corynebacterium glutamicum ATCC 13032 / DSM 20300.

Sigma factors present in regulations.csv (Is_sigma_factor=yes):
  - sigH  (cg0876): Heat shock / oxidative stress
  - sigB  (cg2102): General stress / osmotic stress
  - sigE  (cg1271): Cell envelope stress
  - sigM  (cg3420): Oxidative / superoxide stress
  - cspA2 (cg0371): Cold shock protein (RNA chaperone, not canonical sigma)

Literature sources:
  [1] Ehira S et al. (2009) J Bacteriol 191:4274-4281  — sigH regulon definition
  [2] Schröder J et al. (2003) J Bacteriol 185:4679-4689 — sigH heat shock
  [3] Kim TH et al. (2005) J Bacteriol 187:4565-4572 — sigE envelope stress
  [4] Busche T et al. (2012) J Bacteriol 194:940-952 — sigB general stress
  [5] Frunzke J et al. (2008) Mol Microbiol 67:305-321 — sigM oxidative
  [6] Larisch C et al. (2018) Front Microbiol 9:2774 — sigH/sigE overlap & TSS
  [7] Musiol-Kroll EM et al. (2023) — promoter consensus updated
  [8] CoryneRegNet 7 (2020) Scientific Data 7:142
  [9] Niebisch A & Bott M (2001) Microbiology — heat stress transcriptomics
  [10] Gaigalat L et al. (2007) BMC Genomics — SigB regulon ChIP/expression

Usage:
    python analysis/scripts/enrich_sigma_annotations.py
    python analysis/scripts/enrich_sigma_annotations.py --dry-run

Output:
    data/reference/sigma_factor_annotations.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
DATA_REF = ROOT / "data" / "reference"
REGULATIONS_CSV = DATA_REF / "regulations.csv"

# =============================================================================
# Curated sigma factor annotations for C. glutamicum ATCC 13032
# 
# Fields:
#   locus            : primary cg locus tag
#   alt_locus        : Cgl locus tag
#   gene_name        : standard gene name
#   sigma_class      : ECF sigma or Group 2/3 classification
#   stimulus         : primary activating condition(s)
#   anti_sigma       : anti-sigma factor (if known)
#   anti_sigma_locus : locus of anti-sigma factor
#   consensus_minus35: -35 promoter motif (or None if ECF uses single motif)
#   consensus_minus10: -10 promoter motif
#   spacer_bp        : typical spacer length between -35 and -10 (bp)
#   extended_minus10 : TGn motif at -12/-13 (Actinobacteria feature)
#   heat_stress_activation  : whether activated by temperature upshift
#   tm_activation_degC      : approximate temperature threshold for activation
#   targets_count    : number of target genes in regulations.csv (Is_sigma_factor=yes)
#   key_targets      : curated list of most important direct targets with loci
#   regulon_description: biological summary
#   overlap_with     : other sigma factors with overlapping target sets
#   pmid             : key PMIDs for regulon characterization
# =============================================================================
SIGMA_ANNOTATIONS: Dict[str, Dict[str, Any]] = {
    "sigH": {
        "locus": "cg0876",
        "alt_locus": "Cgl0797",
        "gene_name": "sigH",
        "sigma_class": "ECF_sigma",
        "sigma_group": "σ^H (ECF Group 8)",
        "stimulus": [
            "heat_stress (≥37°C)",
            "oxidative_stress (H₂O₂, diamide, menadione)",
            "disulfide_stress",
            "thiol_depletion",
        ],
        "anti_sigma": "rshA",
        "anti_sigma_locus": "cg0877",
        "anti_sigma_mechanism": "Disulfide bond formation between RshA Cys residues releases SigH under oxidative stress",
        "consensus_minus35": None,
        "consensus_minus10": "TGAC-N4-GTCAA",
        "ecf_consensus": "TGAACC-N17-GCTTGA",
        "ecf_minus35_half": "TGAACC",
        "ecf_minus10_half": "GCTTGA",
        "spacer_bp": 17,
        "extended_minus10": None,
        "heat_stress_activation": True,
        "tm_activation_degC": 37.0,
        "targets_count": 142,
        "key_targets": [
            {"locus": "cg2892", "name": "trxB1", "function": "thioredoxin reductase 1"},
            {"locus": "cg2889", "name": "trxC",  "function": "thioredoxin"},
            {"locus": "cg3219", "name": "groEL1", "function": "chaperonin GroEL"},
            {"locus": "cg1206", "name": "groEL2", "function": "chaperonin GroEL2"},
            {"locus": "cg2987", "name": "groES",  "function": "co-chaperonin GroES"},
            {"locus": "cg0497", "name": "hemA",   "function": "glutamyl-tRNA reductase"},
            {"locus": "cg1364", "name": "atpF",   "function": "ATP synthase subunit b"},
            {"locus": "cg2577", "name": "dnaK",   "function": "Hsp70 chaperone DnaK"},
            {"locus": "cg2985", "name": "dnaJ",   "function": "Hsp40 co-chaperone DnaJ"},
            {"locus": "cg1516", "name": "clpB",   "function": "disaggregase ClpB"},
            {"locus": "cg0371", "name": "cspA2",  "function": "cold shock protein / RNA chaperone"},
            {"locus": "cg2962", "name": "msrA",   "function": "methionine sulphoxide reductase"},
        ],
        "regulon_description": (
            "SigH is the primary heat shock and oxidative stress sigma factor in C. glutamicum. "
            "It controls the largest stress regulon (~142 targets), including thioredoxin systems, "
            "chaperones (GroEL/ES, DnaK/J, ClpB), and detoxification enzymes. "
            "SigH is sequestered by the anti-sigma RshA under reducing conditions; "
            "oxidative stress forms an RshA disulfide bond releasing active SigH. "
            "Under heat shock, SigH-dependent transcription is rapidly activated to induce chaperones and heat tolerance genes as the proteome is remodelled."
        ),
        "overlap_with": ["sigE"],
        "promoter_example": {
            "gene": "trxB1 (cg2892)",
            "sequence": "...TGAACC-N17-GCTTGA...",
            "distance_to_tss": -35,
        },
        "pmid": ["19270092", "25404703", "16385111", "23298179", "12791148", "17449638"],
    },

    "sigB": {
        "locus": "cg2102",
        "alt_locus": "Cgl1926",
        "gene_name": "sigB",
        "sigma_class": "Group_2_sigma",
        "sigma_group": "σ^B (non-essential Group 2)",
        "stimulus": [
            "osmotic_upshift (NaCl, sucrose)",
            "heat_stress (mild, ≥37°C)",
            "stationary_phase_entry",
            "carbon_starvation",
            "ethanol_stress",
        ],
        "anti_sigma": "rsbW",
        "anti_sigma_locus": "cg2101",
        "anti_sigma_mechanism": "RsbW kinase phosphorylates anti-anti-sigma RsbV, sequestering SigB under non-stress conditions",
        "consensus_minus35": "TTGACA",
        "consensus_minus10": "TATAAT",
        "ecf_consensus": None,
        "spacer_bp": 17,
        "extended_minus10": "TGN",
        "heat_stress_activation": True,
        "tm_activation_degC": 38.0,
        "targets_count": 19,
        "key_targets": [
            {"locus": "cg1226", "name": "cg1226", "function": "general stress protein"},
            {"locus": "cg1103", "name": "cg1103", "function": "stress-responsive protein"},
            {"locus": "cg0949", "name": "cg0949", "function": "osmoprotectant-related"},
            {"locus": "cg3327", "name": "cg3327", "function": "betaine biosynthesis"},
            {"locus": "cg2800", "name": "cg2800", "function": "compatible solute transport"},
        ],
        "regulon_description": (
            "SigB controls the general stress response and osmotic adaptation in C. glutamicum. "
            "It regulates a smaller regulon (~19 targets in CoryneRegNet) than SigH, "
            "but contributes broadly to tolerance of multiple stress conditions. "
            "SigB activity increases in stationary phase and under osmotic upshift, "
            "complementing SigH-mediated specific heat shock response."
        ),
        "overlap_with": ["sigH"],
        "pmid": ["15516578", "19376865", "16385111"],
    },

    "sigE": {
        "locus": "cg1271",
        "alt_locus": "Cgl1168",
        "gene_name": "sigE",
        "sigma_class": "ECF_sigma",
        "sigma_group": "σ^E (ECF Group 10, envelope stress)",
        "stimulus": [
            "cell_envelope_stress",
            "SDS",
            "triton_X-100 (detergents)",
            "vancomycin",
            "penicillin (beta-lactams)",
            "EDTA (chelation disrupting outer layer)",
        ],
        "anti_sigma": "rseA",
        "anti_sigma_locus": "cg1272",
        "anti_sigma_mechanism": "RseA (anti-sigma) tethers SigE at the membrane; membrane stress triggers RseA proteolysis releasing SigE",
        "consensus_minus35": None,
        "consensus_minus10": "TGAC-N4-GTCAA",
        "ecf_consensus": "GAACTT-N14-16-GTCAA",
        "ecf_minus35_half": "GAACTT",
        "ecf_minus10_half": "GTCAA",
        "spacer_bp": 15,
        "extended_minus10": None,
        "heat_stress_activation": True,
        "tm_activation_degC": 40.0,
        "targets_count": 1,
        "key_targets": [
            {"locus": "cg1271", "name": "sigE",  "function": "autoregulation"},
            {"locus": "cg1272", "name": "rseA",  "function": "anti-sigma factor"},
            {"locus": "cg0955", "name": "cg0955","function": "mycolyltransferase-related"},
            {"locus": "cg0956", "name": "cg0956","function": "cell wall remodelling"},
            {"locus": "cg1289", "name": "murI",  "function": "glutamate racemase (PG biosynthesis)"},
        ],
        "regulon_description": (
            "SigE controls cell envelope stress response in C. glutamicum. "
            "Although only 1 target appears in CoryneRegNet (due to limited experimental data), "
            "transcriptomic studies reveal ~30+ SigE-controlled genes including mycolyltransferases, "
            "cell wall biosynthesis enzymes, and the anti-sigma rseA itself (autoregulatory loop). "
            "SigE and SigH share overlapping target promoters, with SigE showing preference for "
            "membrane-related genes. Heat stress at 40°C activates SigE through membrane perturbation."
        ),
        "overlap_with": ["sigH"],
        "pmid": ["23298179", "16385111", "12791148"],
    },

    "sigM": {
        "locus": "cg3420",
        "alt_locus": "Cgl3198",
        "gene_name": "sigM",
        "sigma_class": "ECF_sigma",
        "sigma_group": "σ^M (ECF, superoxide/oxidative stress)",
        "stimulus": [
            "superoxide_stress (methyl viologen, paraquat)",
            "menadione",
            "diamide",
            "oxidative_burst",
        ],
        "anti_sigma": "rsmA",
        "anti_sigma_locus": "cg3421",
        "anti_sigma_mechanism": "RsmA anti-sigma; mechanism analogous to SigH/RshA redox switch",
        "consensus_minus35": None,
        "ecf_consensus": "CGAAAC-N16-18-CTGTCA",
        "ecf_minus35_half": "CGAAAC",
        "ecf_minus10_half": "CTGTCA",
        "consensus_minus10": None,
        "spacer_bp": 17,
        "extended_minus10": None,
        "heat_stress_activation": False,
        "tm_activation_degC": None,
        "targets_count": 17,
        "key_targets": [
            {"locus": "cg3420", "name": "sigM",  "function": "autoregulation"},
            {"locus": "cg3421", "name": "rsmA",  "function": "anti-sigma"},
            {"locus": "cg2885", "name": "cg2885","function": "thioredoxin-related"},
            {"locus": "cg2886", "name": "cg2886","function": "oxidoreductase"},
            {"locus": "cg1228", "name": "cg1228","function": "superoxide dismutase-related"},
        ],
        "regulon_description": (
            "SigM is the third ECF sigma factor in C. glutamicum, specialised for superoxide "
            "and severe oxidative stress response. Its regulon (~17 targets) partially overlaps "
            "with SigH, but SigM is the primary responder to superoxide-generating agents like "
            "methyl viologen. Not strongly heat-activated, but contributes to oxidative defence "
            "during prolonged heat stress when ROS generation increases."
        ),
        "overlap_with": ["sigH"],
        "pmid": ["21075931", "16385111"],
    },
}


def count_targets_from_csv(regulations_path: Path) -> Dict[str, int]:
    """Count sigma factor targets from regulations.csv (Is_sigma_factor=yes)."""
    counts: Dict[str, int] = {}
    if not regulations_path.exists():
        return counts
    with open(regulations_path, encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            if row.get("Is_sigma_factor", "").strip().lower() == "yes":
                name = row.get("TF_name", "").strip().lower()
                if name:
                    counts[name] = counts.get(name, 0) + 1
    return counts


def main(dry_run: bool = False) -> None:
    # Cross-check targets counts from CSV
    csv_counts = count_targets_from_csv(REGULATIONS_CSV)
    print("Sigma factor target counts from regulations.csv:")
    for name, count in sorted(csv_counts.items()):
        print(f"  {name}: {count} targets")

    # Update annotations with live counts
    for sigma_id, ann in SIGMA_ANNOTATIONS.items():
        name = ann["gene_name"]
        if name in csv_counts:
            ann["targets_count"] = csv_counts[name]
            print(f"  Updated {sigma_id} targets_count → {csv_counts[name]}")

    # Summary
    print(f"\nSigma factor annotations prepared: {len(SIGMA_ANNOTATIONS)}")
    for sigma_id, ann in SIGMA_ANNOTATIONS.items():
        heat = "[heat-active]" if ann["heat_stress_activation"] else ""
        print(f"  {ann['locus']} ({sigma_id}): class={ann['sigma_class']}, "
              f"targets={ann['targets_count']}, {heat}")

    if dry_run:
        print("\n[DRY RUN] Would write data/reference/sigma_factor_annotations.json")
        return

    out_path = DATA_REF / "sigma_factor_annotations.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(SIGMA_ANNOTATIONS, f, indent=2, ensure_ascii=False)
    print(f"\nWrote {out_path}")
    print("Done. Sigma factor annotations ready for frontend integration.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build sigma factor annotation file for C. glutamicum.")
    parser.add_argument("--dry-run", action="store_true", help="Print summary, do not write files.")
    args = parser.parse_args()
    main(dry_run=args.dry_run)
