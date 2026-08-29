"""Simulate the parsed MSCommunity model on the community (GSP) medium.

Loads data/comm_model_mscommunity_parsed.json, applies the community medium
from analysis.ipynb (acetate 20, nitrate 12, minerals/vitamins), maximizes
community biomass (bio1, 40% 3H11 / 60% R12) with pFBA, and reports growth,
ATP fluxes, the denitrification chain per member, and each member's
uptake/secretion of the shared e0 pool.

Run:  ~/Documents/py_venv/bin/python simulate_community.py
"""
import warnings

warnings.filterwarnings('ignore')

import cobra

MODEL = './data/comm_model_mscommunity_parsed.json'

GSP_MEDIUM = {
    'EX_cpd00001_e0': 1000.0, 'EX_cpd00013_e0': 1000.0,
    'EX_cpd00209_e0': 12.0, 'EX_cpd00029_e0': 20.0,
    'EX_cpd00218_e0': 100.0, 'EX_cpd00220_e0': 100.0,
    'EX_cpd00644_e0': 0.0002281, 'EX_cpd00305_e0': 100.0,
    'EX_cpd00393_e0': 100.0, 'EX_cpd03424_e0': 100.0,
    'EX_cpd00443_e0': 100.0, 'EX_cpd00263_e0': 100.0,
    'EX_cpd00048_e0': 100.0, 'EX_cpd00009_e0': 100.0,
    'EX_cpd00242_e0': 29.759425, 'EX_cpd00205_e0': 1.3415688,
    'EX_cpd00063_e0': 100.0, 'EX_cpd00971_e0': 34.9324073,
    'EX_cpd00099_e0': 100.0, 'EX_cpd00254_e0': 100.0,
    'EX_cpd00030_e0': 100.0, 'EX_cpd00058_e0': 100.0,
    'EX_cpd00034_e0': 100.0, 'EX_cpd10515_e0': 100.0,
    'EX_cpd00149_e0': 100.0, 'EX_cpd00244_e0': 100.0,
    'EX_cpd11574_e0': 100.0, 'EX_cpd15574_e0': 100.0,
    'EX_cpd00067_e0': 100.0,
}

MEMBER = {frozenset({'e0', 'c1'}): '3H11', frozenset({'e0', 'c2'}): 'R12'}


def member_exchange(model, sol, cpd_e0):
    """Net production (+) / consumption (-) of an e0 metabolite per member."""
    met = model.metabolites.get_by_id(cpd_e0)
    totals = {'3H11': 0.0, 'R12': 0.0}
    for r in met.reactions:
        v = sol.fluxes[r.id]
        who = MEMBER.get(frozenset(r.compartments))
        if who and v:
            totals[who] += v * r.metabolites[met]
    return totals


if __name__ == '__main__':
    model = cobra.io.load_json_model(MODEL)
    model.solver = 'glpk'
    for r in model.reactions:
        if r.id.startswith('EX_'):
            r.bounds = (0, 1000)
    for ex_id, v in GSP_MEDIUM.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -v
    model.objective = 'bio1'
    sol = cobra.flux_analysis.pfba(model)
    print(f'pFBA status: {sol.status}')
    print(f'\ncommunity biomass (bio1) : {sol.fluxes["bio1"]:.6f} 1/h')
    print(f'  3H11 biomass (bio2)    : {sol.fluxes["bio2"]:.6f} (40% of bio1 = {0.4 * sol.fluxes["bio1"]:.6f})')
    print(f'  R12  biomass (bio3)    : {sol.fluxes["bio3"]:.6f} (60% of bio1 = {0.6 * sol.fluxes["bio1"]:.6f})')
    print(f'  ATP synthase 3H11/R12  : {sol.fluxes["rxn08173_c1"]:.3f} / {sol.fluxes["rxn08173_c2"]:.3f}')

    print('\ndenitrification chain (flux, mmol/gDW/h):')
    chain = [('Nar (NO3->NO2)', 'dnr00001'), ('Nir (NO2->NO)', 'dnr00002'),
             ('Nor (NO->N2O)', 'dnr00003'), ('Nos (N2O->N2)', 'dnr00004'),
             ('cyt bc1/Azurin', 'dnr00005')]
    print(f'  {"step":<16} {"3H11":>10} {"R12":>10}')
    for name, base in chain:
        v1 = sol.fluxes.get(f'{base}_c1', float('nan'))
        v2 = sol.fluxes.get(f'{base}_c2', float('nan'))
        f1 = f'{v1:.3f}' if f'{base}_c1' in sol.fluxes.index else 'absent'
        f2 = f'{v2:.3f}' if f'{base}_c2' in sol.fluxes.index else 'absent'
        print(f'  {name:<16} {f1:>10} {f2:>10}')

    print('\nshared-pool exchange per member (+ = secreted to e0, - = taken up):')
    key_cpds = [('acetate', 'cpd00029_e0'), ('nitrate', 'cpd00209_e0'),
                ('nitrite', 'cpd00075_e0'), ('NO', 'cpd00418_e0'),
                ('N2O', 'cpd00659_e0'), ('N2', 'cpd00528_e0'),
                ('leucine', 'cpd00107_e0'), ('H+', 'cpd00067_e0'),
                ('CO2', 'cpd00011_e0'), ('NH3', 'cpd00013_e0')]
    print(f'  {"compound":<10} {"3H11":>10} {"R12":>10} {"medium EX":>10}')
    for name, cid in key_cpds:
        if cid not in model.metabolites:
            continue
        t = member_exchange(model, sol, cid)
        ex = sol.fluxes.get(f'EX_{cid}', 0.0)
        print(f'  {name:<10} {t["3H11"]:>10.3f} {t["R12"]:>10.3f} {ex:>10.3f}')

    print('\nall nonzero medium exchanges (EX_*, mmol/gDW/h, - = uptake):')
    for r in sorted(model.reactions, key=lambda x: x.id):
        if r.id.startswith('EX_') and abs(sol.fluxes[r.id]) > 1e-9:
            met = list(r.metabolites)[0]
            print(f'  {r.id:<20} {sol.fluxes[r.id]:>10.4f}  {met.name}')
