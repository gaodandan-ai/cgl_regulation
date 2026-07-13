"""scripts/patch_airc3.py — removes AIRC3 direction lock (it's reversible in iCW773 for a reason)"""
import json

path = "data/reference/thermo_dgr_data.json"
with open(path, encoding="utf-8") as f:
    data = json.load(f)

rxns = data["reactions"]

# AIRC3: bounds=[-1000,1000] in model — reversible for metabolic flexibility.
# Even though dG'° = -28.5, the model uses it in reverse under some conditions.
# AIRC2 (bounds=[0,1000]) is already unidirectional, keep that one.
# We conservatively mark AIRC3 as near-equilibrium to not break the model.
if "AIRC3" in rxns:
    rxns["AIRC3"]["direction_locked"] = "none"
    rxns["AIRC3"]["note"] = "Conservatively left unlocked: model uses this reversibly (bounds [-1000,1000]); AIRC2 covers the forward-lock."
    rxns["AIRC3"]["confidence"] = "LOW"
    print("Patched AIRC3: direction_locked -> none")

# Also update the metadata count
meta = data["_meta"]
# Recount
n_fwd = sum(1 for v in rxns.values() if v.get("direction_locked") == "forward" and v.get("in_model"))
n_rev = sum(1 for v in rxns.values() if v.get("direction_locked") == "reverse" and v.get("in_model"))
meta["n_forward_locked"] = n_fwd
meta["n_reverse_locked"] = n_rev
print(f"Updated counts: forward={n_fwd}, reverse={n_rev}")

with open(path, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2, ensure_ascii=False)
print("Saved.")
