"""Render a cleaned SynCom Escher map JSON to a publication vector PDF.

Draws the community-interaction map (net organism reactions from
make_escher_interactions.py, cleaned by clean_escher_maps.py) with
matplotlib: member-colored reaction edges (3H11 blue, R12 orange -
colorblind-safe categorical pair), a member-colored organism box node per
member (escher-edit's svg_editor.mark_asv_rectangles), neutral metabolite
nodes, chemical-name labels in ink, arrowheads toward products and away from
reactants, and no flux numbers. Vector output on a pure-white background, sized for a
single-column figure.

Run:  ~/Documents/py_venv/bin/python render_escher_map.py \
          [data/escher_syncom_membermap_k1285_cleaned0.json] [out.pdf]
"""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
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
# escher-edit's organism box node (svg_editor.mark_asv_rectangles): a filled
# rectangle over each member's marker backbone, so the member reads as a node
# rather than as a bare convergence of edges. The package's 40x15 px is tuned
# to the ABX viewBox; sized here to this map's units, spanning the +/-20
# multimarkers that the edges attach to.
BOX_W, BOX_H = 100.0, 50.0
BOX_GAP = 12.0        # map units held clear around the box for arrowheads

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

def seg_points(a, b, b1, b2, reverse=False):
    """Control points of one segment: a cubic Bezier when Escher stored
    handles, otherwise a straight line. Reversing the list traverses the same
    curve backwards, which is how an arrowhead is aimed at the other end."""
    pts = [(a['x'], a['y'])]
    if b1 and b2:
        pts += [(b1['x'], b1['y']), (b2['x'], b2['y'])]
    pts.append((b['x'], b['y']))
    return pts[::-1] if reverse else pts


def make_path(pts):
    codes = [Path.MOVETO] + ([Path.CURVE4] * 3 if len(pts) == 4 else [Path.LINETO])
    return Path(pts, codes)


def point_at(pts, t):
    if len(pts) == 2:
        (ax, ay), (bx, by) = pts
        return ax + (bx - ax) * t, ay + (by - ay) * t
    (p0, p1, p2, p3), u = pts, 1 - t
    return (u ** 3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t ** 3 * p3[0],
            u ** 3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t ** 3 * p3[1])


def head_at(pts, t):
    """de Casteljau: the leading part of the segment, up to t."""
    if len(pts) == 2:
        return [pts[0], point_at(pts, t)]
    p0, p1, p2, p3 = pts

    def lerp(u, v):
        return u[0] + (v[0] - u[0]) * t, u[1] + (v[1] - u[1]) * t
    q0, q1, q2 = lerp(p0, p1), lerp(p1, p2), lerp(p2, p3)
    r0, r1 = lerp(q0, q1), lerp(q1, q2)
    return [p0, q0, r0, lerp(r0, r1)]


def clip_to_box(pts, cy):
    """Trim a segment that ends on a marker inside the member's box back to
    the box edge, so the consumption arrowhead is not drawn underneath it.

    FancyArrowPatch ignores shrinkA/shrinkB when handed an explicit ``path``
    (it only honours them on the posA/posB route), so the trim has to happen
    on the geometry rather than at draw time."""
    hw, hh = BOX_W / 2 + BOX_GAP, BOX_H / 2 + BOX_GAP

    def inside(p):
        return abs(p[0]) <= hw and abs(p[1] - cy) <= hh

    if inside(pts[0]) or not inside(pts[-1]):
        return pts
    lo, hi = 0.0, 1.0
    for _ in range(40):                       # first crossing into the box
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if inside(point_at(pts, mid)) else (mid, hi)
    return head_at(pts, hi)

# midmarker y per reaction: the centre of that member's box node, and the
# anchor its label is offset from
mid_ys = {}
for rid, r in body['reactions'].items():
    for seg in r['segments'].values():
        for key in ('from_node_id', 'to_node_id'):
            if nodes[seg[key]]['node_type'] == 'midmarker':
                mid_ys[rid] = nodes[seg[key]]['y']

for rid, r in body['reactions'].items():
    color = MEMBER_COLOR.get(r['bigg_id'], INK)
    cy = mid_ys.get(rid, 0.0)
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
        # arrowhead into the metabolite for products, into the box node for
        # reactants, whose tail is trimmed so the head clears the box
        reverse = (met is b) != producing
        pts = seg_points(a, b, seg.get('b1'), seg.get('b2'), reverse=reverse)
        if not producing:
            pts = clip_to_box(pts, cy)
        arrow = FancyArrowPatch(path=make_path(pts), arrowstyle='-|>',
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

for rid, r in body['reactions'].items():
    ax.add_patch(Rectangle((-BOX_W / 2, mid_ys.get(rid, 0.0) - BOX_H / 2),
                           BOX_W, BOX_H,
                           facecolor=MEMBER_COLOR.get(r['bigg_id'], INK),
                           edgecolor='white', linewidth=0.6, zorder=3.2))

# member labels sit on the far side of each box from the map's centre,
# so neither is overprinted by its own converging edges
centre_y = sum(mid_ys.values()) / len(mid_ys)
for rid, r in body['reactions'].items():
    cy = mid_ys.get(rid, r.get('label_y', 0))
    dy = -(BOX_H / 2 + 55) if cy <= centre_y else BOX_H / 2 + 55
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
