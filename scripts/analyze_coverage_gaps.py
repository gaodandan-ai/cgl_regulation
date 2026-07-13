"""
scripts/analyze_coverage_gaps.py
=================================
Shows all model reactions without thermo data, grouped by subsystem/pathway.
Helps identify which metabolic areas need expansion.
"""
import cobra, json, warnings, re
warnings.filterwarnings('ignore')

model = cobra.io.read_sbml_model('backend/models/iCW773.xml')

with open('data/reference/thermo_dgr_data.json', encoding='utf-8') as f:
    data = json.load(f)

covered_ids = {k for k,v in data['reactions'].items() if v.get('dgr_prime_0') is not None}

# Group uncovered reactions by rough category
uncovered = []
for rxn in model.reactions:
    if rxn.id not in covered_ids:
        # Categorize by ID patterns and subsystem
        subsys = rxn.subsystem if hasattr(rxn, 'subsystem') else ''
        uncovered.append((rxn.id, subsys, rxn.name[:60] if rxn.name else ''))

# Show brief summary grouped by first letters / known pathway
print(f"Total reactions   : {len(model.reactions)}")
print(f"With thermo data  : {len(covered_ids)}")
print(f"Without thermo    : {len(uncovered)}")
print()

# Show non-exchange, non-transport reactions that look metabolic
transport_pat = re.compile(r'(tex|tpp|pp|abc|pts|t2|t3|t4|EX_|sink|DM_)', re.I)
metabolic = [(rid, ss, name) for rid, ss, name in uncovered if not transport_pat.search(rid)]

print(f"Metabolic (non-transport) without data: {len(metabolic)}")
print()

# Print first 200 metabolic reactions without data
print("=== Metabolic reactions without thermo data (first 200) ===")
for rid, ss, name in sorted(metabolic)[:200]:
    print(f"  {rid:35s} | {name[:55]}")

print()
# Also show subsystem breakdown of covered reactions
print("=== Subsystem coverage of COVERED reactions ===")
from collections import Counter
covered_subsys = Counter()
for rxn in model.reactions:
    if rxn.id in covered_ids:
        ss = rxn.subsystem if rxn.subsystem else 'Unknown'
        covered_subsys[ss] += 1
for ss, count in covered_subsys.most_common(20):
    print(f"  {count:3d}  {ss}")
