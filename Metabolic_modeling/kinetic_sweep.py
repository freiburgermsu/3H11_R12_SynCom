"""Kinetic-coefficient sweep of the parsed MSCommunity model on the GSP medium.

Applies the MSCommunity community kinetic constraint (sum of |flux| over each
member's reactions <= K * member biomass flux) through the standalone
MSCommunity package (mscommunity.mscommsim.MSCommunity.add_commkinetics) for a
spectrum of kinetic coefficients K, re-runs pFBA for each, and records each
member's net cellular reaction (membrane-crossing fluxes normalized per gDW
biomass) and its deltaG from the ModelSEED formation energies
(data/compound_formation_energies.json).

Pass `off` among the K values to add a comparison rung: the same model with every
member's kinetic row removed. It runs first, its flux-per-biomass ratios are written
to the metadata as the unconstrained baseline (they used to be a hard-coded constant
from before the member drains were closed), and the output name is still taken from
the numeric ladder alone.

Results: printed summary + data/kinetic_sweep_net_reactions.json (default K
ladder) or data/kinetic_sweep_net_reactions_<min>-<max>.json (CLI K values).

Member biomass drains are closed (close_member_drains=True) and R12's nosZ is knocked
out per Carr et al. 2025; the abundance is the measured proteomic mass fraction. See the
comments at ABUNDANCE and in wrap_community.

With --fitted, the sweep additionally imposes the data-derived constraints:
the community acetate consumption is FIXED (equality, both bounds) at the
rate fitted from the measured SynCom time course
(fit_acetate_uptake.py -> data/fitted_acetate_uptake.json) - every solution
must consume exactly the prescribed acetate - and each member carries its
non-growth maintenance ATP floor - the
calc_max_ATPM values of analysis.ipynb (9.6757 for 3H11, 23.8012 for R12,
mmol/gDW/h in each monoculture) scaled by the member's 0.4/0.6 abundance,
since community-model fluxes are expressed per gDW of total community
biomass. The unscaled monoculture floors sum to 33.5, above the community's
simultaneous maintenance capacity of 25.7 on this medium, and are infeasible.

Run:  ~/Documents/py_venv/bin/python kinetic_sweep.py [--fitted] [off] [K1 K2 ...]
"""
import json
import re
import sys
import warnings
from datetime import date

warnings.filterwarnings('ignore')

import cobra
from mscommunity.commkineticpkg import member_kinetic_reactions
from mscommunity.mscommsim import MSCommunity

from simulate_community import GSP_MEDIUM, MODEL
from net_cell_reactions import (assert_biomass_is_retained, net_exchange,
                                equation_string, FORMATION)

ATPM_MONOCULTURE = {'3H11': 9.6757, 'R12': 23.8012}  # calc_max_ATPM, analysis.ipynb
# Measured steady-state composition, not an assumption: the unique-peptide mass
# fraction of ../Proteomics/data/SynCom-Nitrate-composition-data.csv at 10 mM nitrate
# is 65.8 +/- 1.2 % 3H11 (n = 8), and the whole 1-40 mM series spans 62.3-69.8 %.
# ../Proteomics/data/SynCom-only-composition-data.csv shows the same value reached from
# inocula spanning 5-98 % 3H11 (mean 63.6 %), and Carr et al. 2025 (ISME J 19:wraf093)
# report the convergence to "~65% 3H11 and ~35% R12" as a result of the paper.
# This was 0.4/0.6 -- near enough to inverted -- which is what forced the model to
# discard 3H11 biomass through the drain to satisfy its kinetic row.
ABUNDANCE = {'3H11': 0.658, 'R12': 0.342}

args = sys.argv[1:]
FITTED = '--fitted' in args
args = [a for a in args if a != '--fitted']
# taken out before K_VALUES is built: wrap_community sizes the initial rows from the
# numeric ladder, and the output name must not move when the comparison rung is added
KINETICS_OFF = 'off' in args
args = [a for a in args if a != 'off']
tag = '_fitted' if FITTED else ''
if args:
    K_VALUES = [int(x) for x in args]
    OUT = f'./data/kinetic_sweep_net_reactions{tag}_{K_VALUES[0]}-{K_VALUES[-1]}.json'
else:
    K_VALUES = [100, 150, 250, 400, 600, 800, 1000, 1200, 1500, 2000]
    OUT = f'./data/kinetic_sweep_net_reactions{tag}.json'
LADDER = (['off'] if KINETICS_OFF else []) + K_VALUES
MEMBERS = {'c1': ('3H11', 'bio2'), 'c2': ('R12', 'bio3')}
EPS = 1e-9


def member_reactions(msc, name):
    """Reactions counted in the reported per-member flux sums.

    Delegates to the same helper the CommKinetics constraint is built from, so the
    reported sum|v| and the quantity the constraint bounds cannot drift apart. The
    local version used to exclude SK_/DM_ reactions that the constraint included.
    """
    return list(member_kinetic_reactions(msc, msc.members.get_by_id(name)))


def wrap_community(model):
    """Wrap the loaded community model in the MSCommunity package class,
    preserving the model's 40/60 (3H11/R12) biomass coupling."""
    abundances = {
        '3H11': {'abundance': ABUNDANCE['3H11'],
                 'biomass_compound': model.metabolites.get_by_id('cpd11416_c1')},
        'R12': {'abundance': ABUNDANCE['R12'],
                'biomass_compound': model.metabolites.get_by_id('cpd11416_c2')},
    }
    # construction installs the kinetic rows; start non-binding, each sweep
    # rung replaces them via add_commkinetics.
    #
    # close_member_drains=True is required for this sweep to mean anything: with the
    # member biomass sinks open, a member can synthesise biomass and discard it, so its
    # biomass flux exceeds abundance * bio1 and inflates the right-hand side of its own
    # kinetic row. 83% of 3H11's biomass left through the drain at K = 1285.
    return MSCommunity(model=model, abundances=abundances,
                       kinetic_coeff=max(K_VALUES + [2000]), ID='SynCom',
                       close_member_drains=True)


if __name__ == '__main__':
    model = cobra.io.load_json_model(MODEL)
    model.solver = 'glpk'
    msc = wrap_community(model)
    model = msc.util.model  # the package works on its own copy
    model.solver = 'glpk'
    for r in model.reactions:
        if r.id.startswith('EX_'):
            r.bounds = (0, 1000)
    for ex_id, v in GSP_MEDIUM.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -v
    fitted_info = {}
    if FITTED:
        fitted_info = json.load(open('./data/fitted_acetate_uptake.json'))
        q_ac = fitted_info['q_acetate_mmol_per_gDW_h']
        q_no3 = fitted_info['q_nitrate_mmol_per_gDW_h']
        # calc_max_ATPM (analysis.py) bounded biomass with growth_OD * OD_coeff while
        # bounding uptake with the raw mM drawdown, so the floors carry exactly the same
        # per-litre/per-vessel confusion the uptake fit did and must be rescaled with it.
        # Left unscaled they exceed the ATP the rescaled acetate can supply and the model
        # is infeasible at every kinetic coefficient.
        biomass_scale = fitted_info['biomass_scale_factor']
        atpm = {m: v / biomass_scale for m, v in ATPM_MONOCULTURE.items()}
        model.reactions.EX_cpd00029_e0.bounds = (-q_ac, -q_ac)  # consumption forced
        # Nitrate is capped at its data-fitted rate, which after the gDW/L conversion
        # sits below the GSP medium's own cap of 12 and is therefore the binding
        # constraint on N reduction throughout the sweep (3H11 draws 132.28 mmol
        # NO3 per gDW biomass at mu = 0.01091 /h, i.e. exactly q_no3). Note the
        # nitrate fit is poor (R^2 = 0.107) because the measured NO3 collapses
        # between 42 and 53 h rather than declining with the biomass integral.
        model.reactions.EX_cpd00209_e0.lower_bound = -q_no3
        model.reactions.ATPM_c1.lower_bound = ABUNDANCE['3H11'] * atpm['3H11']
        model.reactions.ATPM_c2.lower_bound = ABUNDANCE['R12'] * atpm['R12']
        print(f"fitted constraints: acetate consumption fixed at {q_ac}, "
              f"nitrate uptake <= {q_no3}, "
              f"ATPM_c1 >= {model.reactions.ATPM_c1.lower_bound:.3f}, "
              f"ATPM_c2 >= {model.reactions.ATPM_c2.lower_bound:.3f} mmol/gDW/h")
    # Carr et al. 2025 (ISME J 19:wraf093) report R12 as "an incomplete denitrifier with
    # a non-functional nosZ, requiring a partner", with 3H11 the primary N2O reducer.
    # Leaving dnr00004_c2 open let R12 carry all of the model's N2 production.
    model.reactions.dnr00004_c2.bounds = (0, 0)

    model.objective = 'bio1'

    formation = json.load(open(FORMATION))['compounds']
    dgf = {}
    for mid, e in formation.items():
        if e.get('seed_id') and e.get('deltag') is not None:
            dgf.setdefault(e['seed_id'], e['deltag'])

    runs = []
    print(f'{"K":>6} {"bio1":>9} {"mu_3H11":>9} {"mu_R12":>9} '
          f'{"dG_3H11":>9} {"dG_R12":>9}  (dG in kcal/gDW biomass)')
    for k in LADDER:
        if k == 'off':
            for member in msc.members:
                cons_id = f'{member.id}_commKin'
                if cons_id in model.constraints:
                    model.remove_cons_vars(model.constraints[cons_id])
        else:
            msc.add_commkinetics(k)  # replaces each member's _commKin row
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
                sumflux = sum(abs(sol.fluxes[r.id]) for r in member_reactions(msc, name))
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
                    # every coefficient below is divided by `mu`, so `mu` must be
                    # biomass that actually reaches the community, not gross synthesis
                    assert_biomass_is_retained(sol, bio_id, ABUNDANCE[name])
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
        b = run['community_biomass']
        print(f'{k:>6} {b if b is not None else float("nan"):>9.5f} '
              f'{row.get("3H11", (float("nan"),))[0]:>9.5f} {row.get("R12", (float("nan"),))[0]:>9.5f} '
              f'{row.get("3H11", (0, float("nan")))[1]:>9.3f} {row.get("R12", (0, float("nan")))[1]:>9.3f}')

    out = {
        '_metadata': {
            'model': MODEL,
            'medium': 'GSP community medium (analysis.ipynb)',
            'method': 'pFBA max bio1 with community kinetic constraint '
                      'sum|v_member| <= K * mu_member, installed per rung through '
                      'mscommunity.mscommsim.MSCommunity.add_commkinetics; net reaction '
                      'and deltaG computed as in net_cell_reactions.py',
            'kinetic_coefficients': LADDER,
            'baseline_unconstrained_flux_per_biomass': next(
                ({m: e['flux_per_biomass'] for m, e in r['members'].items()}
                 for r in runs if r['kinetic_coeff'] == 'off'), None),
            'fitted_constraints': ({'acetate_consumption_fixed_mmol_per_gDW_h': fitted_info['q_acetate_mmol_per_gDW_h'],
                                    'acetate_fit_r_squared': fitted_info['acetate_r_squared'],
                                    'nitrate_uptake_cap_mmol_per_gDW_h': fitted_info['q_nitrate_mmol_per_gDW_h'],
                                    'nitrate_fit_r_squared': fitted_info['nitrate_r_squared'],
                                    'ATPM_floors_mmol_per_gDW_h': {
                                        m: round(ABUNDANCE[m] * ATPM_MONOCULTURE[m]
                                                 / fitted_info['biomass_scale_factor'], 4)
                                        for m in ABUNDANCE},
                                    'biomass_scale_factor': fitted_info['biomass_scale_factor'],
                                    'gdw_per_L_per_OD600': fitted_info['gdw_per_L_per_OD600'],
                                    'ATPM_provenance': 'calc_max_ATPM monoculture values '
                                        '(analysis.ipynb), scaled by the measured abundances '
                                        'and divided by the OD -> gDW/L biomass scale factor'}
                                   if FITTED else None),
            'formation_energies': FORMATION,
            'date': date.today().isoformat(),
        },
        'runs': runs,
    }
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f'\nsaved {OUT}')
