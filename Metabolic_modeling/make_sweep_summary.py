"""Regenerate the kinetic-sweep summary CSV from a sweep result JSON.

Every flux column names the organism performing it: each member's growth, the
deltaG of its net cellular reaction (per gDW biomass and per gDW/h), the nitrogen
species of its normalized net reaction, and its ATP production and consumption
(equal at steady state, both reported). Both fit scores from
evaluate_syncom_fit.py sit beside community growth: the NRMSE and the directional
error. A sweep run with the 'off' rung carries it as the first row -- the same
model with every kinetic row removed, the comparison each budget is read against.

Consecutive rungs are collapsed into one row (kinetic_coeff written as
"1400-2500") once the nitrogen chain stops changing. A rung joins the block before
it only while every nitrogen flux, community and member growth, and both fit scores
match the block's first row to within COLLAPSE_TOL of that column family's largest
magnitude anywhere in the sweep. Those are what the sweep is read for. The budget
can go on reorganising a member's energetics inside a block -- 3H11's kinetic row
stays tight to 1500 while its dissipation walks -- so deltaG, ATP turnover and which
rows bind are not gated but recorded as the block's range in
within_block_variation, never silently dropped. The 'off' row is a comparison, not a
rung, and never folds. Pass --no-collapse to write every rung.

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
# rather than a relative test against the neighbouring value, so a trace value swinging
# over its own tiny range cannot hold a block open while 1% of the family scale still
# catches any change a reader would act on (R12's nitrate uptake falling 4.86 -> 0
# keeps 1350 and 1400 apart).
COLLAPSE_TOL = 0.01
FAMILIES = ('community_biomass', 'growth_1_per_h', 'fit_NRMSE', 'fit_directional',
            'kcal_per_gDW_biomass', 'kcal_per_gDW_h', 'mmol_per_gDW_biomass', 'ATP_')
# the families that define a block; everything else is recorded as its range
GATED = ('community_biomass', 'growth_1_per_h', 'fit_NRMSE', 'fit_directional',
         'mmol_per_gDW_biomass')


def family_of(column):
    for fam in FAMILIES:
        if fam in column:
            return fam
    return column


def collapse(rows):
    """Fold maximal runs of rungs over which the nitrogen chain has stopped changing."""
    scale = {}
    for row in rows:
        for col, v in row.items():
            if col != 'kinetic_coeff' and isinstance(v, (int, float)):
                fam = family_of(col)
                scale[fam] = max(scale.get(fam, 0.0), abs(v))
    blocks, spreads = [], []
    for row in rows:
        if blocks and 'off' not in (row['kinetic_coeff'], blocks[-1][0]['kinetic_coeff']):
            head, varied = blocks[-1][0], {}
            for col, v in row.items():
                if col == 'kinetic_coeff':
                    continue
                ref = head.get(col)
                gated = family_of(col) in GATED
                if isinstance(v, (int, float)) and isinstance(ref, (int, float)):
                    if gated and abs(v - ref) > COLLAPSE_TOL * scale[family_of(col)]:
                        varied = None
                        break
                elif gated and v != ref:
                    # a blank fit never compares equal to a real one
                    varied = None
                    break
                if v != ref:
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
    fit_file = json.load(open(FIT))  # from evaluate_syncom_fit.py
except FileNotFoundError:
    fit_file = {}
fits = fit_file.get('fits', {})
directional = fit_file.get('directional', {})
rows = []
for run in runs:
    k = run['kinetic_coeff']
    a = run['members'].get('3H11', {})
    r = run['members'].get('R12', {})
    n = r.get('normalized_coefficients_mmol_per_gDW_biomass')
    an = a.get('normalized_coefficients_mmol_per_gDW_biomass')
    # columns clustered by organism: community (incl. both fits) | 3H11 | R12;
    # nitrogen species in pathway order, "in" positive for uptake and "out"
    # positive for secretion (a negative value reverses the direction).
    # Which members' Eq. (commkin) rows are tight is recorded for every rung; with the
    # rows removed ('off') there is no budget to be tight against.
    binding = [] if k == 'off' else [
        m for m, v in run['members'].items()
        if v.get('flux_per_biomass', 0) >= (1 - 1e-6) * k]
    row = {
        'kinetic_coeff': k,
        'kinetic_row_binding': '+'.join(sorted(binding)) or 'none',
        'community_biomass_1_per_h': round(run['community_biomass'], 5) if run['community_biomass'] else 0.0,
        'fit_NRMSE_vs_SynCom_timecourse': fits.get(str(k), {}).get('mean', ''),
        'fit_directional_error_vs_SynCom_timecourse': directional.get(str(k), {}).get('mean', ''),
        '3H11_growth_1_per_h': a.get('growth_1_per_h', ''),
        'dG_3H11_kcal_per_gDW_biomass': a.get('deltaG_kcal_per_gDW_biomass', ''),
        'dG_3H11_kcal_per_gDW_h': a.get('deltaG_kcal_per_gDW_h', ''),
        '3H11_NO3_in_mmol_per_gDW_biomass': round(-an.get('Nitrate', 0) + 0.0, 2) if an else '',
        '3H11_NO2_out_mmol_per_gDW_biomass': round(an.get('Nitrite', 0) + 0.0, 2) if an else '',
        '3H11_ATP_production_mmol_per_gDW_h': a.get('atp_production_mmol_gDW_h', ''),
        '3H11_ATP_consumption_mmol_per_gDW_h': a.get('atp_consumption_mmol_gDW_h', ''),
        'R12_growth_1_per_h': r.get('growth_1_per_h', ''),
        'dG_R12_kcal_per_gDW_biomass': r.get('deltaG_kcal_per_gDW_biomass', ''),
        'dG_R12_kcal_per_gDW_h': r.get('deltaG_kcal_per_gDW_h', ''),
        'R12_NO3_in_mmol_per_gDW_biomass': round(-n.get('Nitrate', 0) + 0.0, 2) if n else '',
        'R12_NO2_in_mmol_per_gDW_biomass': round(-n.get('Nitrite', 0) + 0.0, 2) if n else '',
        'R12_N2O_out_mmol_per_gDW_biomass': round(n.get('Nitrous oxide', 0), 2) if n else '',
        'R12_ATP_production_mmol_per_gDW_h': r.get('atp_production_mmol_gDW_h', ''),
        'R12_ATP_consumption_mmol_per_gDW_h': r.get('atp_consumption_mmol_gDW_h', ''),
    }
    rows.append(row)

if COLLAPSE:
    rungs = len(rows)
    rows = collapse(rows)
    if len(rows) < rungs:
        print(f'{rungs} rungs -> {len(rows)} rows (--no-collapse to keep every rung)')
    # the comparison row stands alone: never a block's head, never folded into one
    assert all(not r['collapsed_rungs'] for r in rows if r['kinetic_coeff'] == 'off')
    assert all('off' not in r['collapsed_rungs'].split() for r in rows)

with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print(f'wrote {OUT} ({len(rows)} rows)')
