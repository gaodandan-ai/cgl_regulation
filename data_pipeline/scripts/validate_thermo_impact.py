"""
scripts/validate_thermo_impact.py
===================================
Validates whether thermodynamic directionality pruning
actually changes ecFBA growth/flux predictions.

Runs:
  1. FBA on original model (no pruning)
  2. FBA on thermo-pruned model
  3. Single-gene knockout (KO) ecFBA on both
  4. Reports genes whose growth-impact changes after pruning

Output: data/reference/thermo_validation_report.json
"""

import os, json, warnings, logging
warnings.filterwarnings('ignore')
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("validate_thermo")

import cobra
from cobra.flux_analysis import single_gene_deletion

ROOT_DIR   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MODEL_PATH = os.path.join(ROOT_DIR, "backend", "models", "iCW773.xml")
THERMO_PATH= os.path.join(ROOT_DIR, "data", "reference", "thermo_dgr_data.json")
OUT_PATH   = os.path.join(ROOT_DIR, "data", "reference", "thermo_validation_report.json")

EPSILON_FLUX    = 1e-6     # below this = effectively zero growth
GROWTH_CHANGE_THRESH = 0.01  # 1% change is considered significant

def apply_thermo_pruning(model, thermo_data, epsilon_kJ=1.0):
    """Apply direction locks from thermo data to a cobra model (in-place on a copy)."""
    m = model.copy()
    rxns = thermo_data.get("reactions", {})
    n_locked = 0
    locked_details = []
    for rxn_id, info in rxns.items():
        if info.get("direction_locked") not in ("forward", "reverse"):
            continue
        if not info.get("in_model", False):
            continue
        if rxn_id not in m.reactions:
            continue
        rxn = m.reactions.get_by_id(rxn_id)
        direction = info["direction_locked"]
        if direction == "forward" and rxn.lower_bound < 0:
            rxn.lower_bound = 0
            n_locked += 1
            locked_details.append({"id": rxn_id, "lock": "forward->lb=0"})
        elif direction == "reverse" and rxn.upper_bound > 0:
            rxn.upper_bound = 0
            n_locked += 1
            locked_details.append({"id": rxn_id, "lock": "reverse->ub=0"})
    log.info(f"Applied {n_locked} thermodynamic direction locks")
    return m, n_locked, locked_details


def run_gene_ko_comparison(model_orig, model_pruned, gene_list=None, limit=150):
    """Compare single-gene KO growth between original and pruned model."""
    # Limit to first N genes to keep runtime manageable
    genes = gene_list or [g.id for g in model_orig.genes]
    genes = genes[:limit]

    log.info(f"Running single-gene KO on {len(genes)} genes (original model)...")
    wt_growth = model_orig.slim_optimize()
    wt_pruned = model_pruned.slim_optimize()

    results = []
    changed = []

    for gid in genes:
        # Original model KO
        with model_orig as m:
            m.genes.get_by_id(gid).knock_out()
            gr_orig = m.slim_optimize()
            if gr_orig is None or gr_orig < EPSILON_FLUX:
                essential_orig = True
            else:
                essential_orig = False

        # Pruned model KO
        if gid not in [g.id for g in model_pruned.genes]:
            essential_pruned = essential_orig
            gr_pruned = gr_orig
        else:
            with model_pruned as m:
                m.genes.get_by_id(gid).knock_out()
                gr_pruned = m.slim_optimize()
                if gr_pruned is None or gr_pruned < EPSILON_FLUX:
                    essential_pruned = True
                else:
                    essential_pruned = False

        gr_orig_val  = gr_orig  if gr_orig  is not None else 0.0
        gr_prune_val = gr_pruned if gr_pruned is not None else 0.0

        changed_essentiality = essential_orig != essential_pruned
        rel_change = abs(gr_orig_val - gr_prune_val) / max(wt_growth, 1e-9)
        growth_changed = rel_change > GROWTH_CHANGE_THRESH

        record = {
            "gene":              gid,
            "growth_orig":       round(gr_orig_val,  6),
            "growth_pruned":     round(gr_prune_val, 6),
            "essential_orig":    essential_orig,
            "essential_pruned":  essential_pruned,
            "essentiality_changed": changed_essentiality,
            "growth_rel_change": round(rel_change, 4),
            "growth_changed":    growth_changed or changed_essentiality,
        }
        results.append(record)
        if record["growth_changed"]:
            changed.append(record)
            log.info(f"  CHANGED: {gid}  orig={gr_orig_val:.4f}  pruned={gr_prune_val:.4f}  "
                     f"essential: {essential_orig}->{essential_pruned}")

    return results, changed, wt_growth, wt_pruned


def main():
    log.info("Loading model...")
    model = cobra.io.read_sbml_model(MODEL_PATH)
    log.info(f"Model: {len(model.reactions)} reactions, {len(model.genes)} genes")

    log.info("Loading thermodynamic data...")
    with open(THERMO_PATH, encoding="utf-8") as f:
        thermo_data = json.load(f)

    log.info("Applying thermodynamic pruning...")
    pruned_model, n_locked, locked_details = apply_thermo_pruning(model, thermo_data)

    log.info("Running WT flux analysis...")
    wt_orig   = model.slim_optimize()
    wt_pruned = pruned_model.slim_optimize()
    log.info(f"  WT growth  ORIGINAL: {wt_orig:.4f}")
    log.info(f"  WT growth  PRUNED  : {wt_pruned:.4f}")

    # Run gene KO comparison (first 200 genes for speed)
    all_genes = [g.id for g in model.genes]
    results, changed, wt_orig_val, wt_pruned_val = run_gene_ko_comparison(
        model, pruned_model, all_genes, limit=200
    )

    report = {
        "summary": {
            "wt_growth_original": round(wt_orig_val, 6),
            "wt_growth_pruned":   round(wt_pruned_val, 6),
            "n_reactions_locked": n_locked,
            "n_genes_tested":     len(results),
            "n_genes_changed":    len(changed),
            "n_newly_essential":  sum(1 for r in changed if r["essential_pruned"] and not r["essential_orig"]),
            "n_rescued":          sum(1 for r in changed if not r["essential_pruned"] and r["essential_orig"]),
        },
        "locked_reactions": locked_details,
        "changed_genes":    changed,
        "all_results":      results,
    }

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    log.info(f"Report saved: {OUT_PATH}")

    # Print summary
    s = report["summary"]
    print("\n" + "="*60)
    print(f"WT growth (original): {s['wt_growth_original']:.4f}")
    print(f"WT growth (pruned)  : {s['wt_growth_pruned']:.4f}")
    print(f"Reactions locked    : {s['n_reactions_locked']}")
    print(f"Genes tested        : {s['n_genes_tested']}")
    print(f"Genes with changed prediction: {s['n_genes_changed']}")
    print(f"  Newly essential   : {s['n_newly_essential']}")
    print(f"  Rescued (no longer essential): {s['n_rescued']}")
    print("="*60)

if __name__ == "__main__":
    main()
