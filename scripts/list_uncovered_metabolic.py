"""
scripts/list_uncovered_metabolic.py
=====================================
Lists all uncovered metabolic reaction IDs in iCW773 (non-transport, non-exchange)
to help identify which ones have known thermodynamic data in literature.
"""
import cobra, json, re, warnings
warnings.filterwarnings('ignore')

model = cobra.io.read_sbml_model('backend/models/iCW773.xml')

with open('data/reference/thermo_dgr_data.json', encoding='utf-8') as f:
    data = json.load(f)

covered = {k for k,v in data['reactions'].items() if v.get('dgr_prime_0') is not None}

# Patterns indicating transport / exchange / biomass
skip_pat = re.compile(r'(tex|tpp|abcpp|ptspp|t2pp|t3pp|t4pp|t2rpp|t3ipp|t3ipp|t7pp|EX_|DM_|sink|BIOMASS|biomass|Stonex|t2_|tonex)', re.I)

metabolic_uncovered = []
for rxn in model.reactions:
    if rxn.id in covered:
        continue
    if skip_pat.search(rxn.id):
        continue
    metabolic_uncovered.append((rxn.id, rxn.name))

print(f"Metabolic reactions without thermo data: {len(metabolic_uncovered)}")
print()
for rid, name in sorted(metabolic_uncovered):
    print(f"{rid:35s}  {name[:65]}")
