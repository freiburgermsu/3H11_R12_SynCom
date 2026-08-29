"""Net cell reaction and its deltaG for each SynCom member.

Re-runs the pFBA simulation of the parsed MSCommunity model on the community
(GSP) medium (same setup as simulate_community.py), then for each member
(3H11 = c1, R12 = c2):

1. sums every flux crossing the member's membrane (reactions spanning
   {e0, cX}) per extracellular compound -> the member's net exchange with the
   shared pool, plus its biomass production = the net reaction of that cell;
2. normalizes all coefficients by the member's biomass flux (coefficients are
   then mmol per gDW of new biomass, biomass coefficient = 1);
3. computes deltaG of the net reaction from the ModelSEED formation energies
   in data/compound_formation_energies.json (sum of nu_i * deltaGf_i);
   compounds without a formation energy (biomass, trace metals) are excluded
   and reported.

Results are printed and saved to data/net_cell_reactions.json.

Run:  ~/Documents/py_venv/bin/python net_cell_reactions.py
"""
import json
import re
import warnings
from datetime import date

warnings.filterwarnings('ignore')

import cobra

from simulate_community import GSP_MEDIUM, MODEL

FORMATION = './data/compound_formation_energies.json'
OUT = './data/net_cell_reactions.json'
MEMBERS = {'c1': '3H11', 'c2': 'R12'}
EPS = 1e-9


def net_exchange(model, sol, comp):
    """Net production (+) / consumption (-) of each e0 compound by member comp."""
    totals = {}
    for met in model.metabolites:
        if met.compartment != 'e0':
            continue
        acc = 0.0
        for r in met.reactions:
            if comp in r.compartments and 'e0' in r.compartments:
                acc += sol.fluxes[r.id] * r.metabolites[met]
        if abs(acc) > EPS:
            totals[met.id] = acc
    return totals


def element_balance(net, model):
    """Element/charge imbalance of the net exchange (products - substrates)."""
    balance = {}
    charge = 0.0
    for cid, v in net.items():
        met = model.metabolites.get_by_id(cid)
        charge += (met.charge or 0) * v
        for el, n in (met.elements or {}).items():
            balance[el] = balance.get(el, 0.0) + n * v
    balance = {el: round(n, 6) for el, n in balance.items() if abs(n) > 1e-6}
    return balance, round(charge, 6)


def equation_string(coeffs, names):
    subs = [f'{-v:.4g} {names[c]}' for c, v in sorted(coeffs.items(), key=lambda x: x[1]) if v < 0]
    prods = [f'{v:.4g} {names[c]}' for c, v in sorted(coeffs.items(), key=lambda x: -x[1]) if v > 0]
    return ' + '.join(subs) + ' --> ' + ' + '.join(prods)


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
    print(f'pFBA status: {sol.status}, community biomass {sol.fluxes["bio1"]:.6f}\n')

    formation = json.load(open(FORMATION))['compounds']
    dgf = {}   # base cpd id -> deltaGf (kcal/mol), first instance wins
    for mid, e in formation.items():
        if e.get('seed_id') and e.get('deltag') is not None:
            dgf.setdefault(e['seed_id'], e['deltag'])

    results = {}
    biomass_rxn = {'c1': 'bio2', 'c2': 'bio3'}
    for comp, name in MEMBERS.items():
        growth = sol.fluxes[biomass_rxn[comp]]
        net = net_exchange(model, sol, comp)
        balance, charge = element_balance(net, model)
        names = {c: model.metabolites.get_by_id(c).name.replace(' [e0]', '') for c in net}
        names['biomass'] = f'biomass_{name}'

        norm = {c: v / growth for c, v in net.items()}
        norm['biomass'] = 1.0

        dg_known = 0.0
        excluded = []
        for cid, v in norm.items():
            if cid == 'biomass':
                excluded.append('biomass (no deltaGf)')
                continue
            base = re.match(r'^(cpd\d+)', cid).group(1)
            if base in dgf:
                dg_known += v * dgf[base] / 1000.0  # mmol -> mol
            else:
                excluded.append(f'{names[cid]} ({base}, nu={v:.3g}, no deltaGf)')

        results[name] = {
            'growth_1_per_h': growth,
            'net_exchange_mmol_gDW_h': {names[c]: round(v, 6) for c, v in net.items()},
            'normalized_coefficients_mmol_per_gDW_biomass': {names[c]: round(v, 4) for c, v in norm.items()},
            'net_reaction_normalized': equation_string(norm, names),
            'deltaG_kcal_per_gDW_biomass': round(dg_known, 3),
            'deltaG_kcal_per_gDW_h': round(dg_known * growth, 4),
            'excluded_from_deltaG': excluded,
            'element_imbalance_products_minus_substrates': balance,
            'charge_imbalance': charge,
        }
        print(f'=== {name} (growth {growth:.6f} 1/h) ===')
        print('net cell reaction (per 1 gDW biomass):')
        print(' ', results[name]['net_reaction_normalized'])
        print(f'  deltaG = {dg_known:.2f} kcal/gDW biomass '
              f'({dg_known * growth:.3f} kcal/gDW/h at the pFBA growth rate)')
        print(f'  excluded from deltaG: {excluded}')
        print(f'  element imbalance (≈ biomass drain): {balance}, charge: {charge}\n')

    out = {
        '_metadata': {
            'model': MODEL,
            'medium': 'GSP community medium (analysis.ipynb)',
            'method': 'pFBA max bio1; member net exchange = sum of fluxes of reactions '
                      'spanning {e0, cX} per e0 compound; normalized by member biomass flux; '
                      'deltaG = sum nu_i * deltaGf_i (ModelSEED, kcal/mol, pH 7)',
            'formation_energies': FORMATION,
            'date': date.today().isoformat(),
        },
        'members': results,
    }
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f'saved {OUT}')
