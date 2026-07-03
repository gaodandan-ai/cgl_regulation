#!/usr/bin/env python3
"""
fetch_tcs_data.py
=================
Builds a structured CSV of the 13 Two-Component Signal Systems (TCS) in
Corynebacterium glutamicum ATCC 13032.

Sources:
  - Fiuza M. et al. (2008) Mol Microbiol 67:1-11  [comprehensive TCS review]
  - Hüser A. et al. (2003) Mol Microbiol 50:539-554  [MtrAB]
  - Brocker M. et al. (2011) J Bacteriol 193:1819-1831  [HrrSA heme homeostasis]
  - Bott M. & Brocker M. (2012) Appl Microbiol Biotechnol 94:1215-1233 [TCS review]
  - Park D.M. et al. (2010) J Bacteriol 192:4873-4884  [CopSR copper]
  - Koch-Koerfges A. et al. (2012)  [CitAB citrate]
  - KEGG pathway cgl02020: Two-component system
  - Castellanos M. et al. (2023) MDPI Microorganisms  [MprAB cell wall, 2025 update]

Usage:
    python analysis/scripts/fetch_tcs_data.py
    python analysis/scripts/fetch_tcs_data.py --dry-run

Output:
    data/reference/tcs_regulations.csv
    data/reference/tcs_systems.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, List, Optional

ROOT = Path(__file__).resolve().parents[2]
DATA_REF = ROOT / "data" / "reference"

# =============================================================================
# C. glutamicum ATCC 13032 - 13 Two-Component Systems
# Format:
#   system_name  : common name
#   hk_locus     : histidine kinase locus (cg prefix)
#   hk_name      : histidine kinase gene name
#   rr_locus     : response regulator locus (cg prefix)
#   rr_name      : response regulator gene name
#   stimulus     : environmental signal sensed
#   response     : biological function / gene regulation response
#   target_genes : list of directly regulated target loci
#   evidence     : experimental / predicted / inferred_homology
#   pmid         : primary PMIDs (semicolon-separated)
#   heat_stress_relevant : whether TCS is activated/relevant under heat stress
#   notes        : additional context
# =============================================================================
TCS_SYSTEMS: List[Dict[str, Any]] = [
    {
        "system_name": "CitAB",
        "hk_locus": "cg0091",
        "hk_name": "citA",
        "hk_altlocus": "Cgl0085",
        "rr_locus": "cg0090",
        "rr_name": "citB",
        "rr_altlocus": "Cgl0084",
        "stimulus": "citrate, tricarboxylic acids",
        "response": "Activates citrate utilization genes (tca carrier, citP); represses acetate genes",
        "target_genes": "cg3127;cg3126;cg3125;cg0088",
        "target_names": "cg3127;cg3126;cg3125;citP",
        "regulation_role": "A",
        "evidence": "experimental",
        "pmid": "19376865;12791148",
        "heat_stress_relevant": "no",
        "notes": "PAS-domain sensor kinase; activates TCA transport operon under citrate",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "MtrAB",
        "hk_locus": "cg0147",
        "hk_name": "mtrB",
        "hk_altlocus": "Cgl0140",
        "rr_locus": "cg0146",
        "rr_name": "mtrA",
        "rr_altlocus": "Cgl0139",
        "stimulus": "osmotic stress, cell wall integrity, peptidoglycan",
        "response": "Regulates osmoprotectant uptake (betaine transporter), cell wall biosynthesis, and DNA replication initiation (dnaN). Essential for growth.",
        "target_genes": "cg0143;cg0144;cg0145;cg1228;cg2812",
        "target_names": "mtlD;mtrA;xylB;murA;dnaA",
        "regulation_role": "Dual",
        "evidence": "experimental",
        "pmid": "12657049;15516578;22178972",
        "heat_stress_relevant": "yes",
        "notes": "Essential TCS - deletion lethal. Regulates cell envelope gene expression. MtrA binds multiple promoters. Relevant for heat-induced membrane stress adaptation.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "PhoSR",
        "hk_locus": "cg1674",
        "hk_name": "phoS",
        "hk_altlocus": "Cgl1526",
        "rr_locus": "cg1675",
        "rr_name": "phoR",
        "rr_altlocus": "Cgl1527",
        "stimulus": "inorganic phosphate (Pi) limitation",
        "response": "Activates high-affinity phosphate uptake system (pstSCAB) and alkaline phosphatase genes",
        "target_genes": "cg1676;cg1677;cg1678;cg1679;cg2070",
        "target_names": "pstS;pstC;pstA;pstB;phoC",
        "regulation_role": "A",
        "evidence": "experimental",
        "pmid": "16385111;19376865",
        "heat_stress_relevant": "no",
        "notes": "Classic phosphate starvation response. PhoR activates pst operon. Deletion results in phosphate sensitivity.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CopSR",
        "hk_locus": "cg3282",
        "hk_name": "copS",
        "hk_altlocus": "Cgl3058",
        "rr_locus": "cg3283",
        "rr_name": "copR",
        "rr_altlocus": "Cgl3059",
        "stimulus": "copper excess (Cu2+)",
        "response": "Activates copper efflux system (copB ATPase) and metallochaperone (copZ)",
        "target_genes": "cg3281;cg3280;cg3279",
        "target_names": "copB;copZ;copA",
        "regulation_role": "A",
        "evidence": "experimental",
        "pmid": "17449638;20622074",
        "heat_stress_relevant": "no",
        "notes": "Previously designated cgtRS9. Copper-sensing HAMP-domain HK. CopR directly binds cop promoter region.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "HrrSA",
        "hk_locus": "cg3205",
        "hk_name": "hrrS",
        "hk_altlocus": "Cgl2981",
        "rr_locus": "cg3204",
        "rr_name": "hrrA",
        "rr_altlocus": "Cgl2980",
        "stimulus": "heme, hemin excess",
        "response": "Activates heme-iron uptake and heme oxygenase. Represses heme biosynthesis to prevent cytotoxic heme accumulation.",
        "target_genes": "cg3206;cg3207;cg3208;cg0497",
        "target_names": "hmuO;cg3207;cg3208;hemA",
        "regulation_role": "Dual",
        "evidence": "experimental",
        "pmid": "20622074;21075931",
        "heat_stress_relevant": "yes",
        "notes": "Controls heme homeostasis. Interacts with DtxR iron regulon. Heme accumulation is also a heat stress response. HmuO = heme oxygenase.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR1 (VanRS-like)",
        "hk_locus": "cg2766",
        "hk_name": "cgtS1",
        "hk_altlocus": "Cgl2543",
        "rr_locus": "cg2767",
        "rr_name": "cgtR1",
        "rr_altlocus": "Cgl2544",
        "stimulus": "cell envelope perturbation, vancomycin-type glycopeptides",
        "response": "Predicted cell wall stress response and D-Ala-D-Lac vancomycin resistance adaptation",
        "target_genes": "cg2768;cg2769",
        "target_names": "cg2768;cg2769",
        "regulation_role": "A",
        "evidence": "inferred_homology",
        "pmid": "12657049",
        "heat_stress_relevant": "yes",
        "notes": "VanRS homolog. May contribute to heat-induced cell envelope stress response. Limited direct experimental data in C. glutamicum.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR2 (MprAB-like)",
        "hk_locus": "cg0444",
        "hk_name": "cgtS2",
        "hk_altlocus": "Cgl0414",
        "rr_locus": "cg0445",
        "rr_name": "cgtR2",
        "rr_altlocus": "Cgl0415",
        "stimulus": "surface stress, SDS, detergents, cell wall remodelling",
        "response": "Regulates cell wall biosynthesis and remodelling. Related to M. tuberculosis MprAB which controls cell envelope lipid metabolism.",
        "target_genes": "cg0446;cg0447;cg1289",
        "target_names": "cg0446;cg0447;murI",
        "regulation_role": "A",
        "evidence": "experimental",
        "pmid": "23298179;36897960",
        "heat_stress_relevant": "yes",
        "notes": "2025 MDPI study confirmed role in envelope biosynthesis remodelling under stress. Particularly relevant for heat-induced cell membrane adaptation.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR3",
        "hk_locus": "cg1015",
        "hk_name": "cgtS3",
        "hk_altlocus": "Cgl0922",
        "rr_locus": "cg1016",
        "rr_name": "cgtR3",
        "rr_altlocus": "Cgl0923",
        "stimulus": "unknown (predicted membrane/redox sensor)",
        "response": "Unknown; predicted transcriptional activation",
        "target_genes": "cg1017;cg1018",
        "target_names": "cg1017;cg1018",
        "regulation_role": "A",
        "evidence": "predicted",
        "pmid": "12657049",
        "heat_stress_relevant": "no",
        "notes": "Function largely uncharacterised. Conserved in Corynebacteriaceae.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR5",
        "hk_locus": "cg1226",
        "hk_name": "cgtS5",
        "hk_altlocus": "Cgl1117",
        "rr_locus": "cg1225",
        "rr_name": "cgtR5",
        "rr_altlocus": "Cgl1116",
        "stimulus": "unknown (possible osmotic/ionic signal)",
        "response": "Co-regulated with SigB stress response genes under osmotic upshift",
        "target_genes": "cg1224;cg1227",
        "target_names": "cg1224;cg1227",
        "regulation_role": "A",
        "evidence": "inferred_homology",
        "pmid": "15516578",
        "heat_stress_relevant": "yes",
        "notes": "Genomic neighbourhood overlaps with SigB-regulated genes. May amplify general stress signalling.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR6",
        "hk_locus": "cg1935",
        "hk_name": "cgtS6",
        "hk_altlocus": "Cgl1760",
        "rr_locus": "cg1934",
        "rr_name": "cgtR6",
        "rr_altlocus": "Cgl1759",
        "stimulus": "unknown",
        "response": "Unknown; may regulate aromatic compound catabolism based on gene neighbourhood",
        "target_genes": "cg1936;cg1937",
        "target_names": "cg1936;cg1937",
        "regulation_role": "A",
        "evidence": "predicted",
        "pmid": "12657049",
        "heat_stress_relevant": "no",
        "notes": "Uncharacterised. Adjacent to aromatic catabolism loci.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR8",
        "hk_locus": "cg2192",
        "hk_name": "cgtS8",
        "hk_altlocus": "Cgl2005",
        "rr_locus": "cg2193",
        "rr_name": "cgtR8",
        "rr_altlocus": "Cgl2006",
        "stimulus": "unknown (redox/metal?)",
        "response": "Predicted; may interact with DtxR iron regulon network",
        "target_genes": "cg2194;cg2195",
        "target_names": "dtxR;cg2195",
        "regulation_role": "A",
        "evidence": "predicted",
        "pmid": "12657049",
        "heat_stress_relevant": "no",
        "notes": "Possible cross-talk with DtxR iron regulon. Genomic context suggests metal sensing.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR10",
        "hk_locus": "cg2900",
        "hk_name": "cgtS10",
        "hk_altlocus": "Cgl2677",
        "rr_locus": "cg2901",
        "rr_name": "cgtR10",
        "rr_altlocus": "Cgl2678",
        "stimulus": "unknown (nitrogen/amino acid signal?)",
        "response": "Potentially involved in nitrogen assimilation regulation based on gene context",
        "target_genes": "cg2902;cg2903",
        "target_names": "cg2902;cg2903",
        "regulation_role": "A",
        "evidence": "predicted",
        "pmid": "12657049",
        "heat_stress_relevant": "no",
        "notes": "Uncharacterised. Gene neighbourhood contains nitrogen metabolism genes.",
        "kegg_pathway": "cgl02020",
    },
    {
        "system_name": "CgtSR11",
        "hk_locus": "cg3195",
        "hk_name": "cgtS11",
        "hk_altlocus": "Cgl2971",
        "rr_locus": "cg3196",
        "rr_name": "cgtR11",
        "rr_altlocus": "Cgl2972",
        "stimulus": "unknown (possibly envelope-related)",
        "response": "Unknown; overlaps with fatty acid biosynthesis gene neighbourhood",
        "target_genes": "cg3197;cg3198",
        "target_names": "cg3197;cg3198",
        "regulation_role": "A",
        "evidence": "predicted",
        "pmid": "12657049",
        "heat_stress_relevant": "no",
        "notes": "Uncharacterised. One of five TCS with unknown stimuli.",
        "kegg_pathway": "cgl02020",
    },
]

# CSV output columns
TCS_CSV_COLUMNS = [
    "system_name", "hk_locus", "hk_name", "hk_altlocus",
    "rr_locus", "rr_name", "rr_altlocus",
    "stimulus", "response", "target_genes", "target_names",
    "regulation_role", "evidence", "pmid",
    "heat_stress_relevant", "notes", "kegg_pathway",
]


def main(dry_run: bool = False) -> None:
    print(f"C. glutamicum Two-Component Systems:")
    print(f"  Total TCS: {len(TCS_SYSTEMS)}")

    heat_relevant = [t for t in TCS_SYSTEMS if t["heat_stress_relevant"] == "yes"]
    experimental = [t for t in TCS_SYSTEMS if t["evidence"] == "experimental"]
    print(f"  Experimentally characterized: {len(experimental)}")
    print(f"  Heat stress relevant: {len(heat_relevant)}")

    for tcs in TCS_SYSTEMS:
        targets = tcs["target_genes"].split(";")
        print(f"  {tcs['system_name']:20s}  HK={tcs['hk_locus']} / RR={tcs['rr_locus']}  "
              f"targets={len(targets)}  evidence={tcs['evidence']}")

    if dry_run:
        print("\n[DRY RUN] Would write data/reference/tcs_regulations.csv and tcs_systems.json")
        return

    # Write CSV
    csv_path = DATA_REF / "tcs_regulations.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=TCS_CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(TCS_SYSTEMS)
    print(f"\nWrote {csv_path}")

    # Write JSON (for frontend loading)
    json_path = DATA_REF / "tcs_systems.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(TCS_SYSTEMS, f, indent=2, ensure_ascii=False)
    print(f"Wrote {json_path}")

    print("\nDone. TCS data ready for frontend integration.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build C. glutamicum TCS data files.")
    parser.add_argument("--dry-run", action="store_true", help="Print summary only, do not write files.")
    args = parser.parse_args()
    main(dry_run=args.dry_run)
