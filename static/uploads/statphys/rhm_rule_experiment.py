"""Two-layer rule-world experiment; Python standard library only.

Run from the repository root: python3 static/uploads/statphys/rhm_rule_experiment.py
A world has fixed independently sampled root and shared leaf rules; new data are
independent conditional on the world.  The Wald test is approximate and its
null-world rejection rate is reported explicitly.
"""
import json
import math
import random
from pathlib import Path

SEED = 20261030
V, M1, M2, WORLDS, REPEATS = 4, 4, 4, 1500, 12
SIZES = (100, 400, 1200)  # total samples, half per class
rng = random.Random(SEED)

def draw_world():
    roots = [rng.sample(range(V * V), M1) for _ in range(2)]
    rules = [rng.sample(range(V * V), M2) for _ in range(V)]
    probabilities = []
    for root in roots:
        indices = [r // V for r in root]
        # p(leftmost leaf=0), p(leftmost pair=(0,0))
        leaf = sum(sum(pair // V == 0 for pair in rules[i]) / M2 for i in indices) / M1
        pair = sum(sum(pair == 0 for pair in rules[i]) / M2 for i in indices) / M1
        probabilities.append((leaf, pair))
    return probabilities

def observed_counts(p_leaf, p_pair, n):
    # Joint sampling: event (0,0) implies first leaf=0.
    leaf = pair = 0
    for _ in range(n):
        u = rng.random()
        leaf += u < p_leaf
        pair += u < p_pair
    return leaf, pair

def reject(a, b, n):
    p = (a + b) / (2 * n)
    if p == 0 or p == 1:
        return False
    z = abs(a - b) / n / math.sqrt(2 * p * (1 - p) / n)
    return z > 1.95996398454

worlds = [draw_world() for _ in range(WORLDS)]
metrics = {}
for kind, idx in [('leaf', 0), ('pair', 1)]:
    delta = [a[idx] - b[idx] for a, b in worlds]
    ordered = sorted(abs(d) for d in delta)
    metrics[kind] = {'mean_square': sum(d*d for d in delta)/WORLDS,
        'zero_fraction': sum(abs(d) < 1e-12 for d in delta)/WORLDS,
        'median_abs': ordered[WORLDS//2], 'p90_abs': ordered[int(.9*WORLDS)],
        'power': {}, 'null_rejection': {}, 'bins': {}}

# Reuse each generated dataset for the two tests.
for D in SIZES:
    n=D//2
    stats = {kind: {'success': [], 'null': [], 'bins': [[],[],[],[]]} for kind in metrics}
    for world in worlds:
        successes = {kind: 0 for kind in metrics}
        for _ in range(REPEATS):
            a = observed_counts(*world[0], n)
            b = observed_counts(*world[1], n)
            for kind, idx in [('leaf',0),('pair',1)]:
                successes[kind] += reject(a[idx],b[idx],n)
        for kind, idx in [('leaf',0),('pair',1)]:
            p=successes[kind]/REPEATS
            delta=abs(world[0][idx]-world[1][idx])
            bucket=0 if delta < 1e-12 else 1 if delta < .1 else 2 if delta < .2 else 3
            stats[kind]['success'].append(p)
            stats[kind]['bins'][bucket].append(p)
            if bucket == 0: stats[kind]['null'].append(p)
    for kind in metrics:
        def avg(xs): return sum(xs)/len(xs) if xs else None
        xs=stats[kind]['success']; mean=avg(xs)
        se=math.sqrt(sum((x-mean)**2 for x in xs)/(len(xs)-1)/len(xs))
        metrics[kind]['power'][str(D)]={'mean':mean,'world_se':se}
        metrics[kind]['null_rejection'][str(D)]=avg(stats[kind]['null'])
        metrics[kind]['bins'][str(D)]=[{'worlds':len(xs),'rejection':avg(xs)} for xs in stats[kind]['bins']]

out=Path('static/uploads/statphys')
out.mkdir(parents=True,exist_ok=True)
(out/'rhm_rule_experiment.json').write_text(json.dumps({'seed':SEED,'v':V,'m1':M1,'m2':M2,'worlds':WORLDS,'repeats':REPEATS,'metrics':metrics},indent=2)+'\n')
# Stacked panels keep labels readable on narrow screens.
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="560" height="670" viewBox="0 0 560 670" role="img" aria-label="RHM signal distribution and detection rate">',
'<rect width="560" height="670" fill="white"/><g font-family="Arial,sans-serif" fill="#1e293b">',
'<text x="35" y="32" font-size="18" font-weight="bold">Fixed rule worlds: signal and detection</text>',
'<text x="35" y="55" font-size="12">v=4, m1=m2=4; 1500 worlds; 12 datasets each</text>',
'<text x="48" y="90" font-size="15" font-weight="bold">A. Absolute true signal |delta|</text>',
'<text x="48" y="365" font-size="15" font-weight="bold">B. Approximate two-sided test</text>']
colors={'leaf':'#2563eb','pair':'#d97706'}
for j,(kind,idx) in enumerate([('leaf',0),('pair',1)]):
    ds=[abs(a[idx]-b[idx]) for a,b in worlds]
    buckets=[sum(lo <= x < hi for x in ds)/WORLDS for lo,hi in [(0,.0000001),(.0000001,.1),(.1,.2),(.2,1.01)]]
    for k,val in enumerate(buckets):
        x=80+k*105+j*30; h=val*200
        svg.append(f'<rect x="{x}" y="{310-h:.2f}" width="25" height="{h:.2f}" fill="{colors[kind]}"/>')
    for k,D in enumerate(SIZES):
        val=metrics[kind]['power'][str(D)]['mean']; x=100+k*130+j*34; h=val*200
        svg.append(f'<rect x="{x}" y="{590-h:.2f}" width="28" height="{h:.2f}" fill="{colors[kind]}"/>')
for x,lab in [(85,'0'),(178,'(0,.1)'),(274,'[.1,.2)'),(385,'>=.2')]:
    svg.append(f'<text x="{x}" y="330" font-size="12">{lab}</text>')
for x,lab in [(110,'100'),(240,'400'),(370,'1200')]: svg.append(f'<text x="{x}" y="610" font-size="12">D={lab}</text>')
for y in (310,590):
    svg.append(f'<line x1="48" x2="500" y1="{y}" y2="{y}" stroke="#475569"/>')
    for dy,val in [(0,'0'),(-100,'.5'),(-200,'1')]: svg.append(f'<text x="25" y="{y+dy+4}" font-size="11">{val}</text>')
svg += ['<rect x="100" y="637" width="12" height="12" fill="#2563eb"/><text x="118" y="648" font-size="12">one leaf</text>',
'<rect x="245" y="637" width="12" height="12" fill="#d97706"/><text x="263" y="648" font-size="12">ordered pair (0,0)</text>',
'</g></svg>']
(out/'rhm_rule_experiment.svg').write_text(''.join(svg)+'\n')
print(json.dumps(metrics,indent=2))
