import time, urllib.request, json
time.sleep(10)
r = urllib.request.urlopen('http://localhost:8000/api/thermo/pruning-report')
d = json.loads(r.read())
print('n_pruned :', d['n_pruned'])
print('forward  :', d['n_forward_locked'])
print('reverse  :', d['n_reverse_locked'])
print('near_eq  :', d['n_near_equilibrium'])
print('coverage :', d['data_coverage_pct'])
print('\nAll pruned reactions:')
for x in d.get('top_pruned', []):
    print(f"  {x['reaction_id']:22s} {x['direction']:8s} dG0={x['dgr_prime_0']}")
