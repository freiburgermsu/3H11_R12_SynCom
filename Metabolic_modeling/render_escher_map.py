"""Render a cleaned SynCom Escher map JSON to a publication vector PDF.

Draws the community-interaction map (net organism reactions from
make_escher_interactions.py, cleaned by clean_escher_maps.py) with
matplotlib: member-colored reaction edges (3H11 blue, R12 orange -
colorblind-safe categorical pair), neutral metabolite nodes, chemical-name
labels in ink, arrowheads toward products and away from reactants, and no
flux numbers. Vector output on a pure-white background, sized for a
single-column figure.

Run:  ~/Documents/py_venv/bin/python render_escher_map.py \
          [data/escher_syncom_membermap_k1285_cleaned0.json] [out.pdf]
"""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path

from escher_edit.build_map import cross_feeding_segments
from escher_edit.svg_editor import tint_color   # the package's fade

SRC = sys.argv[1] if len(sys.argv) > 1 else './data/escher_syncom_membermap_k1285_cleaned0.json'
OUT = sys.argv[2] if len(sys.argv) > 2 else SRC.replace('_cleaned0.json', '.pdf')

MEMBER_COLOR = {'3H11': '#2a78d6', 'R12': '#eb6834'}   # categorical slots 1-2
INK, MUTED = '#222222', '#8a8985'
FIG_W = 3.35          # inches, single column
FONT_MET, FONT_MEMBER = 6.5, 9.0
LW = 0.9

# escher-edit highlight styling (svg_editor.apply_color_highlights):
# the nitrogen-transformation species keep full member color with filled
# nodes and bold labels; every other exchange is tinted toward white
HIGHLIGHT = {'Nitrate', 'Nitrite', 'Nitrous oxide', 'N2'}
TINT_FACTOR = 0.5     # EscherStyle default
MUTED_LABEL = '#' + tint_color(INK, TINT_FACTOR)
DPI = 300             # raster companion resolution

escher_map = json.load(open(SRC))
meta, body = escher_map
nodes = body['nodes']

# escher-edit reads cross-feeding off the coefficient signs: a node consumed
# by one member and produced by another. Its ids come back s-prefixed for the
# rendered SVG (Escher's JSON has no per-segment style), so strip that to
# match the JSON segment keys.
CROSS_FED = {s.lstrip('s') for s in cross_feeding_segments(escher_map)}

xs = [n['x'] for n in nodes.values()]
ys = [n['y'] for n in nodes.values()]
pad = 260.0
x0, x1 = min(xs) - pad, max(xs) + pad
y0, y1 = min(ys) - pad * 0.55, max(ys) + pad * 0.55

fig_h = FIG_W * (y1 - y0) / (x1 - x0)
fig, ax = plt.subplots(figsize=(FIG_W, fig_h))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

def seg_path(a, b, b1, b2, reverse=False):
    pts = [(a['x'], a['y'])]
    codes = [Path.MOVETO]
    if b1 and b2:
        pts += [(b1['x'], b1['y']), (b2['x'], b2['y']), (b['x'], b['y'])]
        codes += [Path.CURVE4] * 3
    else:
        pts.append((b['x'], b['y']))
        codes.append(Path.LINETO)
    if reverse:
        pts = pts[::-1]
    return Path(pts, codes)

for r in body['reactions'].values():
    color = MEMBER_COLOR.get(r['bigg_id'], INK)
    coef = {m['bigg_id']: m['coefficient'] for m in r['metabolites']}
    for seg_id, seg in r['segments'].items():
        a, b = nodes[seg['from_node_id']], nodes[seg['to_node_id']]
        met = b if b['node_type'] == 'metabolite' else (a if a['node_type'] == 'metabolite' else None)
        if met is None:                       # marker-to-marker backbone
            ax.plot([a['x'], b['x']], [a['y'], b['y']], color=color, lw=LW * 1.4,
                    solid_capstyle='round', zorder=2)
            continue
        producing = coef.get(met['bigg_id'], 0) > 0
        highlighted = met['bigg_id'] in HIGHLIGHT
        seg_color = color if highlighted else '#' + tint_color(color, TINT_FACTOR)
        # arrowhead into the metabolite for products, into the cluster for
        # reactants; shrink at the metabolite end so heads sit off the node
        reverse = (met is b) != producing
        path = seg_path(a, b, seg.get('b1'), seg.get('b2'), reverse=reverse)
        arrow = FancyArrowPatch(path=path, arrowstyle='-|>',
                                mutation_scale=6.5 if highlighted else 6,
                                lw=LW * 1.2 if highlighted else LW,
                                color=seg_color, shrinkA=0,
                                shrinkB=5 if producing else 3,
                                zorder=2.5 if highlighted else 2, fill=True,
                                linestyle=(0, (3, 2)) if seg_id in CROSS_FED else 'solid')
        ax.add_patch(arrow)

for n in nodes.values():
    if n['node_type'] == 'metabolite':
        highlighted = n['bigg_id'] in HIGHLIGHT
        ax.scatter(n['x'], n['y'], s=17 if highlighted else 14,
                   facecolor=INK if highlighted else 'white',
                   edgecolor=INK if highlighted else MUTED,
                   linewidth=0.7, zorder=3.5 if highlighted else 3)
        ha = 'right' if n['x'] < 0 else 'left'
        dx = -14 if n['x'] < 0 else 14
        ax.annotate(n['bigg_id'], (n['x'] + dx, n['y']),
                    color=INK if highlighted else MUTED_LABEL,
                    fontweight='bold' if highlighted else 'normal',
                    fontsize=FONT_MET, ha=ha, va='center', zorder=4)

# member labels sit on the far side of each cluster from the map's centre,
# so neither is overprinted by its own converging edges
mid_ys = {}
for rid, r in body['reactions'].items():
    for seg in r['segments'].values():
        for key in ('from_node_id', 'to_node_id'):
            if nodes[seg[key]]['node_type'] == 'midmarker':
                mid_ys[rid] = nodes[seg[key]]['y']
centre_y = sum(mid_ys.values()) / len(mid_ys)
for rid, r in body['reactions'].items():
    cy = mid_ys.get(rid, r.get('label_y', 0))
    dy = -60 if cy <= centre_y else 60
    ax.annotate(r['bigg_id'], (0, cy + dy),
                color=MEMBER_COLOR.get(r['bigg_id'], INK),
                fontsize=FONT_MEMBER, fontweight='bold',
                ha='center', va='center', zorder=5)

ax.set_xlim(x0, x1)
ax.set_ylim(y1, y0)              # Escher y grows downward
ax.set_aspect('equal')
ax.axis('off')
fig.savefig(OUT, bbox_inches='tight', pad_inches=0.02, facecolor='white', dpi=DPI)
png = OUT.rsplit('.', 1)[0] + '.png'
fig.savefig(png, bbox_inches='tight', pad_inches=0.02, facecolor='white', dpi=DPI)
print(f'wrote {OUT} and {png} ({DPI} dpi, {len(CROSS_FED)} cross-fed edges dashed)')
