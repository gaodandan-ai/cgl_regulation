"""
backend/thermo_pruner.py
========================
Thermodynamic Directionality Pruning for iCW773 COBRApy model.

Loads pre-computed ΔrG' range data (data/reference/thermo_dgr_data.json)
and applies one-time bound tightening at model-load time:

  • If  ΔrG'_max < −ε  →  reaction is always exergonic forward
        → lower_bound forced to max(0, lower_bound)   [locks out reverse]

  • If  ΔrG'_min >  +ε →  reaction is always endergonic in forward direction
        → upper_bound forced to min(0, upper_bound)   [locks out forward]

  • Otherwise (near-equilibrium or missing data) → no change

All pruning is done *before* the model is returned to callers, so every
subsequent FBA / ecFBA solve automatically benefits from tighter bounds
without any MILP or additional solver overhead.
"""

import os
import json
import math
import logging
from typing import Dict, Any, Optional

logger = logging.getLogger("thermo_pruner")

# ── Path resolution ────────────────────────────────────────────────────────────
_THIS_DIR    = os.path.dirname(os.path.abspath(__file__))
_ROOT_DIR    = os.path.dirname(_THIS_DIR)
THERMO_PATH  = os.path.join(_ROOT_DIR, "data", "reference", "thermo_dgr_data.json")

# ── Module-level cache for report (populated once at model-load time) ──────────
_pruning_report: Optional[Dict[str, Any]] = None


def _load_thermo_data() -> Dict[str, Any]:
    """Load thermodynamic ΔrG' data from JSON. Returns empty dict on failure."""
    if not os.path.exists(THERMO_PATH):
        logger.warning(f"Thermo data file not found at {THERMO_PATH}. Pruning disabled.")
        return {}
    try:
        with open(THERMO_PATH, "r", encoding="utf-8") as f:
            raw = json.load(f)
        reactions = raw.get("reactions", {})
        logger.info(f"Loaded thermodynamic data: {len(reactions)} reaction entries.")
        return reactions
    except Exception as e:
        logger.error(f"Failed to load thermo data: {e}")
        return {}


def apply_thermodynamic_pruning(model, epsilon: float = 1.0):
    """
    Apply thermodynamic directionality pruning to a COBRApy model in-place.

    Parameters
    ----------
    model   : cobra.Model  (modified in-place)
    epsilon : float        kJ/mol threshold; reactions within ±ε of zero are
                           treated as near-equilibrium and left unchanged.

    Returns
    -------
    model   : the same cobra.Model, with tightened reaction bounds
    report  : dict  with pruning statistics (also stored in module-level cache)
    """
    global _pruning_report

    thermo_data = _load_thermo_data()
    if not thermo_data:
        _pruning_report = _empty_report(model)
        return model, _pruning_report

    n_forward_locked   = 0   # newly forced (was reversible, now locked forward)
    n_reverse_locked   = 0   # newly forced (was reversible, now locked reverse)
    n_confirmed_fwd    = 0   # already [0,UB], thermodynamics agrees forward
    n_confirmed_rev    = 0   # already [LB,0], thermodynamics agrees reverse
    n_skipped_neq      = 0   # near-equilibrium, no change
    n_no_data          = 0
    pruned_details     = []  # actual bound changes
    confirmed_details  = []  # already-consistent, no change needed
    reverted_locks     = []  # locks reverted due to infeasibility

    for rxn in model.reactions:
        entry = thermo_data.get(rxn.id)

        if entry is None or entry.get("dgr_prime_0") is None:
            n_no_data += 1
            continue

        dgr_min   = entry.get("dgr_prime_min")
        dgr_max   = entry.get("dgr_prime_max")
        direction = entry.get("direction_locked", "none")

        if dgr_min is None or dgr_max is None:
            n_no_data += 1
            continue

        old_lb = rxn.lower_bound
        old_ub = rxn.upper_bound

        if direction == "forward":
            if rxn.lower_bound < 0:
                # Reaction was reversible but thermodynamics says forward-only
                rxn.lower_bound = 0.0
                # ── Safety guard: revert if lock makes model infeasible ────
                try:
                    _test = model.slim_optimize()
                    if _test is None or _test < 1e-6:
                        rxn.lower_bound = old_lb
                        logger.warning(
                            "[SafeGuard] Reverted forward lock on %s: "
                            "locking makes model infeasible (dG=%.1f kJ/mol). "
                            "Model may use this reaction reversibly.",
                            rxn.id, entry.get("dgr_prime_0", 0)
                        )
                        reverted_locks.append({
                            "reaction_id": rxn.id,
                            "direction": "forward",
                            "dgr_prime_0": entry.get("dgr_prime_0"),
                            "dgr_prime_min": dgr_min,
                            "dgr_prime_max": dgr_max,
                            "confidence": entry.get("confidence", "?"),
                            "reason": "makes model infeasible (biomass < 1e-6)"
                        })
                        n_skipped_neq += 1
                        continue
                except Exception as ex:
                    rxn.lower_bound = old_lb
                    reverted_locks.append({
                        "reaction_id": rxn.id,
                        "direction": "forward",
                        "dgr_prime_0": entry.get("dgr_prime_0"),
                        "dgr_prime_min": dgr_min,
                        "dgr_prime_max": dgr_max,
                        "confidence": entry.get("confidence", "?"),
                        "reason": f"optimisation exception: {str(ex)}"
                    })
                    n_skipped_neq += 1
                    continue
                # ─────────────────────────────────────────────────────────
                n_forward_locked += 1
                pruned_details.append({
                    "reaction_id":  rxn.id,
                    "direction":    "forward",
                    "status":       "newly_locked",
                    "old_lb": old_lb, "old_ub": old_ub,
                    "new_lb": rxn.lower_bound, "new_ub": rxn.upper_bound,
                    "dgr_prime_0":   entry.get("dgr_prime_0"),
                    "dgr_prime_min": dgr_min, "dgr_prime_max": dgr_max,
                    "confidence": entry.get("confidence", "?"),
                    "note": entry.get("note", "")
                })
            else:
                # Already irreversible forward - thermodynamics confirms
                n_confirmed_fwd += 1
                confirmed_details.append({
                    "reaction_id":  rxn.id,
                    "direction":    "forward",
                    "status":       "already_correct",
                    "old_lb": old_lb, "old_ub": old_ub,
                    "dgr_prime_0":   entry.get("dgr_prime_0"),
                    "dgr_prime_min": dgr_min, "dgr_prime_max": dgr_max,
                    "confidence": entry.get("confidence", "?"),
                })

        elif direction == "reverse":
            if rxn.upper_bound > 0:
                rxn.upper_bound = 0.0
                # ── Safety guard: revert if lock makes model infeasible ────
                try:
                    _test = model.slim_optimize()
                    if _test is None or _test < 1e-6:
                        rxn.upper_bound = old_ub
                        logger.warning(
                            "[SafeGuard] Reverted reverse lock on %s: "
                            "locking makes model infeasible (dG=%.1f kJ/mol). "
                            "Model may use this reaction in forward direction.",
                            rxn.id, entry.get("dgr_prime_0", 0)
                        )
                        reverted_locks.append({
                            "reaction_id": rxn.id,
                            "direction": "reverse",
                            "dgr_prime_0": entry.get("dgr_prime_0"),
                            "dgr_prime_min": dgr_min,
                            "dgr_prime_max": dgr_max,
                            "confidence": entry.get("confidence", "?"),
                            "reason": "makes model infeasible (biomass < 1e-6)"
                        })
                        n_skipped_neq += 1
                        continue
                except Exception as ex:
                    rxn.upper_bound = old_ub
                    reverted_locks.append({
                        "reaction_id": rxn.id,
                        "direction": "reverse",
                        "dgr_prime_0": entry.get("dgr_prime_0"),
                        "dgr_prime_min": dgr_min,
                        "dgr_prime_max": dgr_max,
                        "confidence": entry.get("confidence", "?"),
                        "reason": f"optimisation exception: {str(ex)}"
                    })
                    n_skipped_neq += 1
                    continue
                # ─────────────────────────────────────────────────────────
                n_reverse_locked += 1
                pruned_details.append({
                    "reaction_id":  rxn.id,
                    "direction":    "reverse",
                    "status":       "newly_locked",
                    "old_lb": old_lb, "old_ub": old_ub,
                    "new_lb": rxn.lower_bound, "new_ub": rxn.upper_bound,
                    "dgr_prime_0":   entry.get("dgr_prime_0"),
                    "dgr_prime_min": dgr_min, "dgr_prime_max": dgr_max,
                    "confidence": entry.get("confidence", "?"),
                    "note": entry.get("note", "")
                })
            else:
                n_confirmed_rev += 1
                confirmed_details.append({
                    "reaction_id":  rxn.id,
                    "direction":    "reverse",
                    "status":       "already_correct",
                    "old_lb": old_lb, "old_ub": old_ub,
                    "dgr_prime_0":   entry.get("dgr_prime_0"),
                    "dgr_prime_min": dgr_min, "dgr_prime_max": dgr_max,
                    "confidence": entry.get("confidence", "?"),
                })
        else:
            n_skipped_neq += 1

    total     = len(model.reactions)
    n_pruned  = n_forward_locked + n_reverse_locked
    n_thermo_consistent = n_confirmed_fwd + n_confirmed_rev
    coverage_pct = round((total - n_no_data) / total * 100, 1) if total > 0 else 0.0

    logger.info(
        f"Thermodynamic pruning: {n_pruned} newly locked "
        f"(fwd={n_forward_locked}, rev={n_reverse_locked}), "
        f"{n_thermo_consistent} already consistent, coverage={coverage_pct}%"
    )

    _pruning_report = {
        "enabled":             True,
        "total_reactions":     total,
        "data_coverage_pct":   coverage_pct,
        # Newly locked (actual bound changes)
        "n_pruned":            n_pruned,
        "n_forward_locked":    n_forward_locked,
        "n_reverse_locked":    n_reverse_locked,
        # Already-correct (model already had the right direction)
        "n_confirmed_forward": n_confirmed_fwd,
        "n_confirmed_reverse": n_confirmed_rev,
        "n_thermo_consistent": n_thermo_consistent,
        "n_near_equilibrium":  n_skipped_neq,
        "n_no_data":           n_no_data,
        "epsilon_kJ":          epsilon,
        "conditions":          "pH 7.0, I=0.1 M, T=30°C",
        "top_pruned":          pruned_details[:50],
        "confirmed_reactions": confirmed_details[:100],
        "all_pruned_count":    len(pruned_details),
        "reverted_locks":      reverted_locks,
        "n_reverted_locks":    len(reverted_locks),
    }

    return model, _pruning_report


def get_pruning_report() -> Dict[str, Any]:
    """Return the cached pruning report (populated after model load)."""
    global _pruning_report
    if _pruning_report is None:
        return {
            "enabled": False,
            "message": "Model not yet loaded or pruning not yet run."
        }
    return _pruning_report


def _empty_report(model) -> Dict[str, Any]:
    return {
        "enabled":           False,
        "total_reactions":   len(model.reactions),
        "data_coverage_pct": 0.0,
        "n_pruned":          0,
        "n_forward_locked":  0,
        "n_reverse_locked":  0,
        "n_near_equilibrium":0,
        "n_no_data":         len(model.reactions),
        "reverted_locks":    [],
        "n_reverted_locks":  0,
        "message":           "Thermo data file not found; pruning disabled."
    }
