#!/usr/bin/env python3
"""Enumerate all 81 labeled two-node signed graphs and generate their atlas.

Classes are structural filters under the assumptions stated in main.md. They
are not the results of a parameter scan or assertions of realized dynamics.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures'; OUT.mkdir(exist_ok=True)
PAIRS=[(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]
SIGN={0:'0',1:'+',-1:'-'}
COLORS=['#d9e5ec','#d1c4e9','#bbdefb','#a5d6a7','#ffe0b2','#ffcc80','#b2dfdb','#ef9a9a']

def family(selfs,cross):
    has_positive=1 in selfs
    if cross==(0,0):return 2 if has_positive else 1
    if 0 in cross:return 4 if has_positive else 3
    if cross[0]*cross[1]>0:return 6 if has_positive else 5
    return 8 if has_positive else 7

def swap(t):
    sx,sy,xy,yx=t
    return sy,sx,yx,xy

records=[]; grid=np.zeros((9,9),int)
for row,selfs in enumerate(PAIRS):
    for col,cross in enumerate(PAIRS):
        graph=selfs+cross;group=family(selfs,cross);grid[row,col]=group
        records.append(dict(id=9*row+col+1,self_x=selfs[0],self_y=selfs[1],x_to_y=cross[0],y_to_x=cross[1],family=f'R{group}',positive_cycle=(1 in selfs or cross[0]*cross[1]>0),hopf_sign_candidate=(group==8)))
counts=Counter(r['family'] for r in records)
assert len(records)==81 and len({(r['self_x'],r['self_y'],r['x_to_y'],r['y_to_x']) for r in records})==81
assert [counts[f'R{i}'] for i in range(1,9)]==[4,5,16,20,8,10,8,10]
assert len({min(t,swap(t)) for t in product([0,1,-1],repeat=4)})==45
assert sum(r['positive_cycle'] for r in records)==53
assert sum(r['hopf_sign_candidate'] for r in records)==10
(ROOT/'code/topology_inventory.json').write_text(json.dumps({'edge_order':['self_x','self_y','x_to_y','y_to_x'],'counts':dict(sorted(counts.items())),'labeled_count':81,'unlabeled_count':45,'records':records},indent=2)+'\n')

plt.rcParams.update({'font.size':10,'savefig.bbox':'tight'})
fig,ax=plt.subplots(figsize=(6.5,5.6))
ax.imshow(grid,cmap=ListedColormap(COLORS),norm=BoundaryNorm(np.arange(.5,9.5),8))
for row in range(9):
    for col in range(9):ax.text(col,row,f'R{grid[row,col]}',ha='center',va='center',fontsize=10)
labels=[''.join(SIGN[v] for v in pair) for pair in PAIRS]
ax.set_xticks(range(9),labels);ax.set_yticks(range(9),labels)
ax.set_xticks(np.arange(-.5,9),minor=True);ax.set_yticks(np.arange(-.5,9),minor=True)
ax.grid(which='minor',color='white',lw=2);ax.tick_params(which='minor',bottom=False,left=False)
ax.set(xlabel='cross edges (X -> Y, Y -> X)',ylabel='self edges (X -> X, Y -> Y)',title='All 81 signed two-node topologies: structural families')
fig.tight_layout();fig.savefig(OUT/'topology_classes.png',dpi=180);plt.close(fig)

# Three PNG atlas sheets; no TeX installation is needed.
from matplotlib.patches import FancyArrowPatch, PathPatch
from matplotlib.path import Path as MplPath

def edge(ax, start, end, sign, rad=0):
    ax.add_patch(FancyArrowPatch(start, end, connectionstyle=f'arc3,rad={rad}',
        arrowstyle='-|>' if sign > 0 else '-[', mutation_scale=12,
        linewidth=1.4, color='#1746b0' if sign > 0 else '#c62828'))

def draw_graph(ax, row, col):
    sx, sy = PAIRS[row]; xy, yx = PAIRS[col]
    group = family((sx, sy), (xy, yx)); ident = 9 * row + col + 1
    for x, label, sign in [(0, 'X', sx), (2, 'Y', sy)]:
        ax.text(x, 0, label, ha='center', va='center', fontsize=12,
            bbox=dict(boxstyle='round,pad=.35', facecolor='#f0effa', edgecolor='#333'))
        if sign:
            # Draw the loop as a curved arrow whose endpoints meet the node's top edge.
            verts = [(x-.16, .24), (x-.8, 1.25), (x+.8, 1.25), (x+.16, .24)]
            path = MplPath(verts, [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
            ax.add_patch(FancyArrowPatch(path=path, arrowstyle='-|>' if sign>0 else '-[',
                mutation_scale=10, linewidth=1.3, color='#1746b0' if sign>0 else '#c62828'))
    if xy: edge(ax, (.32, .05), (1.68, .05), xy, -.25)
    if yx: edge(ax, (1.68, -.05), (.32, -.05), yx, -.25)
    ax.text(1, -.65, f'No. {ident} / R{group}', ha='center', fontsize=9)
    ax.set(xlim=(-.65, 2.65), ylim=(-1, 1.1)); ax.axis('off')

for sheet in range(3):
    fig, axes = plt.subplots(9, 3, figsize=(9, 17))
    for col in range(9):
        for offset, row in enumerate(range(sheet*3, sheet*3+3)):
            ax = axes[col, offset]; draw_graph(ax, row, col)
            if col == 0:
                ax.set_title('self: (' + ', '.join(SIGN[v] for v in PAIRS[row]) + ')', fontsize=11, pad=12)
            if offset == 0:
                ax.text(-1, 0, '(' + ', '.join(SIGN[v] for v in PAIRS[col]) + ')', ha='right', va='center', fontsize=10)
    fig.suptitle(f'Topology atlas {sheet+1}: all signed two-node circuits', fontsize=14, y=.995)
    fig.text(.02, .975, 'Rows: cross edges (X -> Y, Y -> X); columns: self edges (X -> X, Y -> Y).', fontsize=10)
    fig.subplots_adjust(left=.15, right=.99, top=.95, bottom=.01, hspace=.2, wspace=.05)
    fig.savefig(OUT/f'topology_atlas_{sheet+1}.png', dpi=180)
    plt.close(fig)
print('Exhaustive topology counts:',dict(sorted(counts.items())))
print('Unique labeled graphs: 81; node-swap classes: 45; positive-cycle graphs: 53; Hopf sign candidates: 10')
