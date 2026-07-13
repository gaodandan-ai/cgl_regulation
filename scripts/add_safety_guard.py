"""scripts/add_safety_guard.py — inserts FBA feasibility guard into thermo_pruner.py"""

OLD_BLOCK = (
    "        if direction == \"forward\":\n"
    "            if rxn.lower_bound < 0:\n"
    "                # Reaction was reversible but thermodynamics says forward-only\n"
    "                rxn.lower_bound = 0.0\n"
    "                n_forward_locked += 1\n"
    "                pruned_details.append({\n"
    "                    \"reaction_id\":  rxn.id,\n"
    "                    \"direction\":    \"forward\",\n"
    "                    \"status\":       \"newly_locked\",\n"
    "                    \"old_lb\": old_lb, \"old_ub\": old_ub,\n"
    "                    \"new_lb\": rxn.lower_bound, \"new_ub\": rxn.upper_bound,\n"
    "                    \"dgr_prime_0\":   entry.get(\"dgr_prime_0\"),\n"
    "                    \"dgr_prime_min\": dgr_min, \"dgr_prime_max\": dgr_max,\n"
    "                    \"confidence\": entry.get(\"confidence\", \"?\"),\n"
    "                    \"note\": entry.get(\"note\", \"\")\n"
    "                })\n"
)

NEW_BLOCK = (
    "        if direction == \"forward\":\n"
    "            if rxn.lower_bound < 0:\n"
    "                # Reaction was reversible but thermodynamics says forward-only\n"
    "                rxn.lower_bound = 0.0\n"
    "                # ── Safety guard: revert if lock makes model infeasible ────\n"
    "                try:\n"
    "                    _test = model.slim_optimize()\n"
    "                    if _test is None or _test < 1e-6:\n"
    "                        rxn.lower_bound = old_lb\n"
    "                        logger.warning(\n"
    "                            \"[SafeGuard] Reverted forward lock on %s: \"\n"
    "                            \"locking makes model infeasible (dG=%.1f kJ/mol). \"\n"
    "                            \"Model may use this reaction reversibly.\",\n"
    "                            rxn.id, entry.get(\"dgr_prime_0\", 0)\n"
    "                        )\n"
    "                        n_skipped_neq += 1\n"
    "                        continue\n"
    "                except Exception:\n"
    "                    rxn.lower_bound = old_lb\n"
    "                    n_skipped_neq += 1\n"
    "                    continue\n"
    "                # ─────────────────────────────────────────────────────────\n"
    "                n_forward_locked += 1\n"
    "                pruned_details.append({\n"
    "                    \"reaction_id\":  rxn.id,\n"
    "                    \"direction\":    \"forward\",\n"
    "                    \"status\":       \"newly_locked\",\n"
    "                    \"old_lb\": old_lb, \"old_ub\": old_ub,\n"
    "                    \"new_lb\": rxn.lower_bound, \"new_ub\": rxn.upper_bound,\n"
    "                    \"dgr_prime_0\":   entry.get(\"dgr_prime_0\"),\n"
    "                    \"dgr_prime_min\": dgr_min, \"dgr_prime_max\": dgr_max,\n"
    "                    \"confidence\": entry.get(\"confidence\", \"?\"),\n"
    "                    \"note\": entry.get(\"note\", \"\")\n"
    "                })\n"
)

path = "backend/thermo_pruner.py"
with open(path, "r", encoding="utf-8") as f:
    content = f.read()

if OLD_BLOCK in content:
    content = content.replace(OLD_BLOCK, NEW_BLOCK, 1)
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print("SUCCESS: FBA safety guard inserted into thermo_pruner.py")
else:
    print("ERROR: target block not found. File may have changed.")
    # Show what we have around line 104
    lines = content.splitlines()
    for i, line in enumerate(lines[100:120], 101):
        print(f"  {i}: {repr(line)}")
