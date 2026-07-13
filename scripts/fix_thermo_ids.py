"""scripts/fix_thermo_ids.py — patches the 2 remaining ID mismatches in thermo_dgr_data.json"""
import json

path = "data/reference/thermo_dgr_data.json"
with open(path, encoding="utf-8") as f:
    data = json.load(f)

rxns = data["reactions"]

# ICHOR -> ICHORS
if "ICHOR" in rxns and "ICHORS" in rxns:
    # ICHOR is already in model as ICHORS; remove the broken one
    del rxns["ICHOR"]
    print("Removed duplicate ICHOR (correct ID is ICHORS)")
elif "ICHOR" in rxns:
    entry = rxns.pop("ICHOR")
    entry["in_model"] = True
    rxns["ICHORS"] = entry
    print("Renamed ICHOR -> ICHORS")

# HSTP -> HISTP (histidinol-phosphatase, different from HSTP transaminase which is HSTPT)
if "HSTP" in rxns and "HISTP" in rxns:
    del rxns["HSTP"]
    print("Removed HSTP (already covered by HISTP)")
elif "HSTP" in rxns:
    entry = rxns.pop("HSTP")
    entry["in_model"] = True
    rxns["HISTP"] = entry
    print("Renamed HSTP -> HISTP")

with open(path, "w", encoding="utf-8") as f:
    json.dump(data, f, indent=2, ensure_ascii=False)
print("Saved. ID misses remaining:", [k for k,v in rxns.items() if v.get("in_model")==False and v.get("dgr_prime_0") is not None])
