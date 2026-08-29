"""Kinetic-coefficient sweep of the parsed MSCommunity model on the GSP medium.

Applies the MSCommunity/CommKineticPkg community kinetic constraint
(sum of |flux| over each member's reactions <= K * member biomass flux;
implemented directly with optlang since the shipped CommKineticPkg targets an
older MSCommunity API) for a spectrum of 10 kinetic coefficients K, re-runs
pFBA for each, and records each member's net cellular reaction (membrane-
crossing fluxes normalized per gDW biomass) and its deltaG from the ModelSEED
formation energies (data/compound_formation_energies.json).

Baseline (unconstrained) pFBA ratios sum|v|/mu are ~1336 (3H11) and ~1243
(R12), so K spans 100 (strongly limiting) to 2000 (non-binding).

Results: printed summary + data/kinetic_sweep_net_reactions.json (default K
ladder) or data/kinetic_sweep_net_reactions_<min>-<max>.json (CLI K values).

Run:  ~/Documents/py_venv/bin/python kinetic_sweep.py [K1 K2 ...]
"""
import json
import re
import sys
import warnings
from datetime import date

warnings.filterwarnings('ignore')

import cobra
from optlang.symbolics import Zero

from simulate_community import GSP_MEDIUM, MODEL
from net_cell_reactions import net_exchange, equation_string, FORMATION

if len(sys.argv) > 1:
    K_VALUES = [int(x) for x in sys.argv[1:]]
    OUT = f'./data/kinetic_sweep_net_reactions_{K_VALUES[0]}-{K_VALUES[-1]}.json'
else:
    K_VALUES = [100, 150, 250, 400, 600, 800, 1000, 1200, 1500, 2000]
    OUT = './data/kinetic_sweep_net_reactions.json'
MEMBERS = {'c1': ('3H11', 'bio2'), 'c2': ('R12', 'bio3')}
EPS = 1e-9


def member_reactions(model, comp, bio_id):
    return [r for r in model.reactions
            if comp in r.compartments and not r.id.startswith(('EX_', 'SK_'))
            and r.id != bio_id]


def add_kinetic_constraints(model, k):
    constraints = []
    for comp, (_, bio_id) in MEMBERS.items():
        cons = model.problem.Constraint(Zero, ub=0, name=f'commkin_{comp}')
        model.add_cons_vars(cons)
        model.solver.update()
        bio = model.reactions.get_by_id(bio_id)
        coefs = {}
        for r in member_reactions(model, comp, bio_id):
            coefs[r.forward_variable] = 1
            coefs[r.reverse_variable] = 1
        coefs[bio.forward_variable] = -k
        coefs[bio.reverse_variable] = k
        cons.set_linear_coefficients(coefs)
        constraints.append(cons)
    return constraints


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

    formation = json.load(open(FORMATION))['compounds']
    dgf = {}
    for mid, e in formation.items():
        if e.get('seed_id') and e.get('deltag') is not None:
            dgf.setdefault(e['seed_id'], e['deltag'])

    runs = []
    print(f'{"K":>6} {"bio1":>9} {"mu_3H11":>9} {"mu_R12":>9} '
          f'{"dG_3H11":>9} {"dG_R12":>9}  (dG in kcal/gDW biomass)')
    for k in K_VALUES:
        cons = add_kinetic_constraints(model, k)
        try:
            sol = cobra.flux_analysis.pfba(model)
            status = sol.status
        except Exception as exc:
            sol, status = None, f'failed ({exc})'
        run = {'kinetic_coeff': k, 'status': status,
               'community_biomass': sol.fluxes['bio1'] if sol else None, 'members': {}}
        row = {}
        if sol is not None:
            for comp, (name, bio_id) in MEMBERS.items():
                mu = sol.fluxes[bio_id]
                sumflux = sum(abs(sol.fluxes[r.id]) for r in member_reactions(model, comp, bio_id))
                atp = model.metabolites.get_by_id(f'cpd00002_{comp}')
                atp_prod = atp_cons = 0.0
                for r in atp.reactions:
                    t = sol.fluxes[r.id] * r.metabolites[atp]
                    if t > 0:
                        atp_prod += t
                    else:
                        atp_cons -= t
                entry = {'growth_1_per_h': round(mu, 6),
                         'sum_flux_mmol_gDW_h': round(sumflux, 3),
                         'flux_per_biomass': round(sumflux / mu, 1) if mu > EPS else None,
                         'atp_production_mmol_gDW_h': round(atp_prod, 4),
                         'atp_consumption_mmol_gDW_h': round(atp_cons, 4)}
                if mu > EPS:
                    net = net_exchange(model, sol, comp)
                    names = {c: model.metabolites.get_by_id(c).name.replace(' [e0]', '') for c in net}
                    names['biomass'] = f'biomass_{name}'
                    norm = {c: v / mu for c, v in net.items()}
                    norm['biomass'] = 1.0
                    dg = 0.0
                    excluded = []
                    for cid, v in norm.items():
                        if cid == 'biomass':
                            excluded.append('biomass')
                            continue
                        base = re.match(r'^(cpd\d+)', cid).group(1)
                        if base in dgf:
                            dg += v * dgf[base] / 1000.0
                        else:
                            excluded.append(base)
                    entry.update({
                        'net_exchange_mmol_gDW_h': {names[c]: round(v, 6) for c, v in net.items()},
                        'normalized_coefficients_mmol_per_gDW_biomass':
                            {names[c]: round(v, 4) for c, v in norm.items()},
                        'net_reaction_normalized': equation_string(norm, names),
                        'deltaG_kcal_per_gDW_biomass': round(dg, 3),
                        'deltaG_kcal_per_gDW_h': round(dg * mu, 4),
                        'excluded_from_deltaG': excluded,
                    })
                    row[name] = (mu, dg)
                else:
                    entry['note'] = 'no growth at this kinetic coefficient'
                    row[name] = (0.0, float('nan'))
                run['members'][name] = entry
        runs.append(run)
        model.remove_cons_vars(cons)
        b = run['community_biomass']
        print(f'{k:>6} {b if b is not None else float("nan"):>9.5f} '
              f'{row.get("3H11", (float("nan"),))[0]:>9.5f} {row.get("R12", (float("nan"),))[0]:>9.5f} '
              f'{row.get("3H11", (0, float("nan")))[1]:>9.3f} {row.get("R12", (0, float("nan")))[1]:>9.3f}')

    out = {
        '_metadata': {
            'model': MODEL,
            'medium': 'GSP community medium (analysis.ipynb)',
            'method': 'pFBA max bio1 with community kinetic constraint '
                      'sum|v_member| <= K * mu_member (CommKineticPkg formulation, '
                      'implemented via optlang); net reaction and deltaG computed as in '
                      'net_cell_reactions.py',
            'kinetic_coefficients': K_VALUES,
            'baseline_unconstrained_flux_per_biomass': {'3H11': 1335.7, 'R12': 1243.0},
            'formation_energies': FORMATION,
            'date': date.today().isoformat(),
        },
        'runs': runs,
    }
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f'\nsaved {OUT}')
