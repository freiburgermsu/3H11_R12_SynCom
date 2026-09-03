"""Regenerate the kinetic-sweep summary CSV from a sweep result JSON.

Every flux column names the organism performing it: deltaG of each member's
net cellular reaction (per gDW biomass and per gDW/h), the N/S species of
R12's normalized net reaction, and each member's ATP production and
consumption (equal at steady state, both reported).

Run:  ~/Documents/py_venv/bin/python make_sweep_summary.py \
          [data/kinetic_sweep_net_reactions_800-1400.json]
"""
import csv
import json
import sys

SRC = sys.argv[1] if len(sys.argv) > 1 else './data/kinetic_sweep_net_reactions_800-1400.json'
OUT = SRC.replace('kinetic_sweep_net_reactions', 'kinetic_sweep_summary').replace('.json', '.csv')
FIT = SRC.replace('kinetic_sweep_net_reactions', 'kinetic_sweep_fit')

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
    # columns clustered by organism: community | 3H11 | R12 | community fit;
    # nitrogen species in pathway order, "in" positive for uptake and "out"
    # positive for secretion (a negative value reverses the direction)
    row = {
        'kinetic_coeff': run['kinetic_coeff'],
        'community_biomass_1_per_h': round(run['community_biomass'], 5) if run['community_biomass'] else 0.0,
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
        'fit_NRMSE_vs_SynCom_timecourse':
            fits.get(str(run['kinetic_coeff']), {}).get('mean', ''),
    }
    rows.append(row)

with open(OUT, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys())
    w.writeheader()
    w.writerows(rows)
print(f'wrote {OUT} ({len(rows)} rows)')
