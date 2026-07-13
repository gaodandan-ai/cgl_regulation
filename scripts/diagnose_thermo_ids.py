import cobra, json, warnings
warnings.filterwarnings('ignore')

model = cobra.io.read_sbml_model('backend/models/iCW773.xml')
rxn_ids = {r.id for r in model.reactions}

with open('data/reference/thermo_dgr_data.json', encoding='utf-8') as f:
    data = json.load(f)

curated = {k: v for k, v in data['reactions'].items() if v.get('dgr_prime_0') is not None}
matched   = [k for k in curated if k in rxn_ids]
unmatched = [k for k in curated if k not in rxn_ids]

print(f'Model reactions: {len(rxn_ids)}')
print(f'Curated entries: {len(curated)}')
print(f'Matched IDs    : {len(matched)}')
print(f'Unmatched IDs  : {len(unmatched)}')

print('\n=== Matched + direction ===')
for k in sorted(matched):
    d = curated[k]
    print(f"  {k:30s} -> {d['direction_locked']:8s} dG0={d['dgr_prime_0']}")

print('\n=== UNMATCHED (with candidate model IDs) ===')
rxn_list = sorted(rxn_ids)
for k in sorted(unmatched):
    base = k.lower().replace('_reverse','').replace('_num1','').replace('_num2','')
    candidates = [r for r in rxn_list if base in r.lower() or r.lower().startswith(base[:5])][:4]
    print(f"  {k:40s} | candidates: {candidates}")

# Also show a sample of model reaction IDs to understand naming conventions
print('\n=== Sample model reaction IDs (200) ===')
for r in sorted(rxn_ids)[:200]:
    print(' ', r)
