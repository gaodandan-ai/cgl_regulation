"""scripts/diagnose_locks.py — finds which direction lock kills model growth"""
import cobra, json, warnings
warnings.filterwarnings('ignore')

model = cobra.io.read_sbml_model('backend/models/iCW773.xml')
print('WT growth:', round(model.slim_optimize(), 4))

with open('data/reference/thermo_dgr_data.json', encoding='utf-8') as f:
    data = json.load(f)

locked = [(rxn_id, info) for rxn_id, info in data['reactions'].items()
          if info.get('direction_locked') in ('forward','reverse') and info.get('in_model')]

print(f'Total locked reactions: {len(locked)}\n')

for rxn_id, info in locked:
    m = model.copy()
    rxn = m.reactions.get_by_id(rxn_id)
    direction = info['direction_locked']
    old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
    dgr = info['dgr_prime_0']
    conf = info['confidence']

    if direction == 'forward' and rxn.lower_bound < 0:
        rxn.lower_bound = 0
    elif direction == 'reverse' and rxn.upper_bound > 0:
        rxn.upper_bound = 0

    gr = m.slim_optimize()
    if gr is None or gr < 1e-4:
        status = 'KILLS GROWTH'
    else:
        status = 'OK'

    print(f"  {status:14s} {rxn_id:20s} {direction:8s} dG={dgr:7.1f} bounds=[{old_lb},{old_ub}] {conf}")
