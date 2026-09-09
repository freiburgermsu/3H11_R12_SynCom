"""Regenerate the kinetic-sweep summary CSV from a sweep result JSON.

Every flux column names the organism performing it: deltaG of each member's
net cellular reaction (per gDW biomass and per gDW/h), the N/S species of
R12's normalized net reaction, and each member's ATP production and
consumption (equal at steady state, both reported).

Consecutive rungs that report the same answer are collapsed into a single row
(kinetic_coeff written as "1500-2500"), since above the coefficient at which the
kinetic constraint stops binding every further rung reproduces the same optimum
and adds only length. The collapse is measured, not hard-coded: a rung joins the
block only if every column matches the block's first row to within COLLAPSE_TOL
of that column family's largest magnitude anywhere in the sweep, so any rung at
which something genuinely changes keeps its own row. Whatever does vary inside a
block is printed and recorded in the row, never silently dropped. Pass
--no-collapse to write every rung.

Run:  ~/Documents/py_venv/bin/python make_sweep_summary.py \
          [data/kinetic_sweep_net_reactions_800-1400.json] [--no-collapse]
"""
import csv
import json
import sys

args = [a for a in sys.argv[1:] if not a.startswith('--')]
COLLAPSE = '--no-collapse' not in sys.argv[1:]
SRC = args[0] if args else './data/kinetic_sweep_net_reactions_800-1400.json'
OUT = SRC.replace('kinetic_sweep_net_reactions', 'kinetic_sweep_summary').replace('.json', '.csv')
FIT = SRC.replace('kinetic_sweep_net_reactions', 'kinetic_sweep_fit')

# Tolerance is 1% of the largest magnitude the family reaches anywhere in the sweep,
# rather than a relative test against the neighbouring value: a trace column swinging
# over its own tiny range (R12 secretes 0.39 mmol H2S/gDW at one rung and none at the
# next, against nitrogen fluxes of 27-132) is degenerate-optimum reshuffling, while 1%
# of the family scale still catches any change a reader would act on.
COLLAPSE_TOL = 0.01
FAMILIES = ('community_biomass', 'fit_NRMSE', 'kcal_per_gDW_biomass', 'kcal_per_gDW_h',
            'mmol_per_gDW_biomass', 'ATP_')


def family_of(column):
    for fam in FAMILIES:
        if fam in column:
            return fam
    return column


def collapse(rows):
    """Fold maximal runs of rungs that report the same answer into one row."""
    scale = {}
    for row in rows:
        for col, v in row.items():
            if col != 'kinetic_coeff' and isinstance(v, (int, float)):
                fam = family_of(col)
                scale[fam] = max(scale.get(fam, 0.0), abs(v))
    blocks, spreads = [], []
    for row in rows:
        if blocks:
            head, varied = blocks[-1][0], {}
            for col, v in row.items():
                if col == 'kinetic_coeff' or not isinstance(v, (int, float)):
                    continue
                ref = head.get(col)
                if not isinstance(ref, (int, float)):
                    continue
                delta = abs(v - ref)
                if delta > COLLAPSE_TOL * scale[family_of(col)]:
                    varied = None
                    break
                if delta:
                    varied[col] = (ref, v)
            if varied is not None:
                blocks[-1].append(row)
                spreads[-1].update(varied)
                continue
        blocks.append([row])
        spreads.append({})
    out = []
    for block, varied in zip(blocks, spreads):
        row = dict(block[0])
        if len(block) > 1:
            row['kinetic_coeff'] = f"{block[0]['kinetic_coeff']}-{block[-1]['kinetic_coeff']}"
            note = '; '.join(f'{c} {a}->{b}' for c, (a, b) in sorted(varied.items()))
            row['collapsed_rungs'] = ' '.join(str(r['kinetic_coeff']) for r in block)
            row['within_block_variation'] = note or 'none'
            print(f"  collapsed {row['kinetic_coeff']} ({len(block)} rungs): "
                  f"{note or 'every column identical'}")
        else:
            row['collapsed_rungs'] = ''
            row['within_block_variation'] = ''
        out.append(row)
    return out

runs = json.load(open(SRC))['runs']
try:
    fits = json.load(open(FIT))['fits']  # from evaluate_syncom_fit.py
except FileNotFoundError:
    fits = {}
rows = []
for run in runs:
    a = run['members'].get('3H11', {})
    r = run['members'].get('R12', {})
    n = r.get('normalized_coefficients_mmol_per_gDW_biomass')
    an = a.get('normalized_coefficients_mmol_per_gDW_biomass')
    # columns clustered by organism: community (incl. fit) | 3H11 | R12;
    # nitrogen species in pathway order, "in" positive for uptake and "out"
    # positive for secretion (a negative value reverses the direction)
    row = {
        'kinetic_coeff': run['kinetic_coeff'],
        'community_biomass_1_per_h': round(run['community_biomass'], 5) if run['community_biomass'] else 0.0,
        'fit_NRMSE_vs_SynCom_timecourse':
            fits.get(str(run['kinetic_coeff']), {}).get('mean', ''),
        'dG_3H11_kcal_per_gDW_biomass': a.get('deltaG_kcal_per_gDW_biomass', ''),
        'dG_3H11_kcal_per_gDW_h': a.get('deltaG_kcal_per_gDW_h', ''),
        '3H11_NO3_in_mmol_per_gDW_biomass': round(-an.get('Nitrate', 0) + 0.0, 2) if an else '',
        '3H11_NO2_out_mmol_per_gDW_biomass': round(an.get('Nitrite', 0) + 0.0, 2) if an else '',
        '3H11_ATP_production_mmol_per_gDW_h': a.get('atp_production_mmol_gDW_h', ''),
        '3H11_ATP_consumption_mmol_per_gDW_h': a.get('atp_consumption_mmol_gDW_h', ''),
        'dG_R12_kcal_per_gDW_biomass': r.get('deltaG_kcal_per_gDW_biomass', ''),
        'dG_R12_kcal_per_gDW_h': r.get('deltaG_kcal_per_gDW_h', ''),
        'R12_NO3_in_mmol_per_gDW_biomass': round(-n.get('Nitrate', 0) + 0.0, 2) if n else '',
        'R12_NO2_in_mmol_per_gDW_biomass': round(-n.get('Nitrite', 0) + 0.0, 2) if n else '',
        'R12_N2O_out_mmol_per_gDW_biomass': round(n.get('Nitrous oxide', 0), 2) if n else '',
        'R12_N2_out_mmol_per_gDW_biomass': round(n.get('N2', 0), 2) if n else '',
        'R12_H2S_out_mmol_per_gDW_biomass': round(n.get('H2S', 0), 2) if n else '',
        'R12_ATP_production_mmol_per_gDW_h': r.get('atp_production_mmol_gDW_h', ''),
        'R12_ATP_consumption_mmol_per_gDW_h': r.get('atp_consumption_mmol_gDW_h', ''),
    }
    rows.append(row)

if COLLAPSE:
    rungs = len(rows)
    rows = collapse(rows)
    if len(rows) < rungs:
        print(f'{rungs} rungs -> {len(rows)} rows (--no-collapse to keep every rung)')

with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print(f'wrote {OUT} ({len(rows)} rows)')
