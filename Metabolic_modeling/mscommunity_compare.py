"""Build two MSCommunity models of the 3H11/R12 SynCom and compare them.

Version 1 ("parsed"):   members are parsed back out of the downloaded
                        community model (data/comm_model.json), with GPRs
                        grafted from the downloaded isolate models
                        (data/acido_model.json, data/rhoda_model.json),
                        since the custom community builder dropped genes.
Version 2 ("isolates"): members are the repo's gap-filled isolate models
                        (data/model_acido_gf.json, data/model_rhoda_gf.json).

Both are merged with modelseedpy MSCommunity.build_from_species_models
(3H11 -> index 1 / c1, R12 -> index 2 / c2, shared e0) at 40/60 abundance.

Run:  ~/Documents/py_venv/bin/python mscommunity_compare.py
"""
import re
import warnings

warnings.filterwarnings('ignore')

import cobra
from cobra.core import Model, Metabolite, Reaction

from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.fbahelper import FBAHelper

# modelseedpy.community.__init__ imports modules missing from this checkout
# (commkineticpkg, mscompatibility, ...), so load mscommunity.py directly.
import importlib.util
import modelseedpy

_msc_path = modelseedpy.__path__[0] + '/community/mscommunity.py'
_spec = importlib.util.spec_from_file_location('mscommunity_standalone', _msc_path)
_msc_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_msc_mod)
MSCommunity = _msc_mod.MSCommunity

# FBAHelper.add_atp_hydrolysis no longer exists in modelseedpy 0.4.3 but
# CommunityModelSpecies still calls it; delegate to the MSModelUtil method.
if not hasattr(FBAHelper, 'add_atp_hydrolysis'):
    FBAHelper.add_atp_hydrolysis = staticmethod(
        lambda model, compartment: MSModelUtil.get(model).add_atp_hydrolysis(compartment))

# mscommunity.build_from_species_models expects parse_id's compartment index
# as a string (len(output[2])), but this checkout returns an int.
_orig_parse_id = MSModelUtil.parse_id
def _parse_id_compat(cobra_obj):
    out = _orig_parse_id(cobra_obj)
    if out is not None and len(out) == 3:
        return (out[0], out[1], str(out[2]))
    return out
MSModelUtil.parse_id = staticmethod(_parse_id_compat)

DATA = './data'
NAMES = ['3H11', 'R12']
ABUNDANCES = {'3H11': 0.4, 'R12': 0.6}

# community medium from analysis.ipynb (acetate + nitrate + minerals)
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


def parse_member_from_comm(model_comm, token, member_id, gpr_source):
    """Extract one member (token 'A' or 'R') from the community model back
    into a standalone c0/e0 isolate model, grafting GPRs from gpr_source."""
    member = Model(member_id)
    member.compartments = {'c0': 'Cytosol', 'e0': 'Extracellular'}

    def map_met(m):
        if m.id.endswith('_e0'):
            new_id, comp = m.id, 'e0'
        elif m.id.endswith(f'_c{token}'):
            new_id, comp = m.id[:-len(f'_c{token}')] + '_c0', 'c0'
        else:
            return None
        if new_id in member.metabolites:
            return member.metabolites.get_by_id(new_id)
        name = re.sub(r'\s*\[[A-Za-z0-9]+\]$', '', m.name or new_id)
        met = Metabolite(new_id, m.formula, f'{name} [{comp}]', m.charge, comp)
        member.add_metabolites([met])
        return met

    n_rxn = 0
    for r in model_comm.reactions:
        if r.id == f'bio1_{token}':
            new_id = 'bio1'
        elif r.id.endswith(f'_c{token}'):
            new_id = r.id[:-len(f'_c{token}')] + '_c0'
        elif r.id.endswith(f'_e{token}'):
            new_id = r.id[:-len(f'_e{token}')] + '_e0'
        else:
            continue  # other member, community bio1, or EX_
        name = re.sub(r'\s*\[[A-Za-z0-9]+\]$', '', r.name or new_id)
        r_copy = Reaction(new_id, name, r.subsystem, r.lower_bound, r.upper_bound)
        member.add_reactions([r_copy])
        r_copy.add_metabolites({map_met(m): v for m, v in r.metabolites.items()})
        n_rxn += 1

    # exchanges for every extracellular metabolite (comm EX_ bounds carry the
    # saved analysis state, so regenerate them open instead of copying)
    ex = []
    for m in member.metabolites:
        if m.compartment == 'e0':
            r_ex = Reaction(f'EX_{m.id}', f'Exchange for {m.name}', 'exchange', -1000, 1000)
            ex.append((r_ex, m))
    member.add_reactions([r for r, _ in ex])
    for r_ex, m in ex:
        r_ex.add_metabolites({m: -1})

    # reset bounds that were left over from the saved simulation state
    member.reactions.bio1.bounds = (0, 1000)
    if 'ATPM_c0' in member.reactions:
        member.reactions.ATPM_c0.bounds = (0, 1000)

    # graft GPRs from the isolate model
    grafted, no_rule, missing = 0, 0, []
    for r in member.reactions:
        if r.id.startswith('EX_') or r.id == 'bio1':
            continue
        if r.id in gpr_source.reactions:
            rule = gpr_source.reactions.get_by_id(r.id).gene_reaction_rule
            if rule:
                r.gene_reaction_rule = rule
                grafted += 1
            else:
                no_rule += 1
        else:
            missing.append(r.id)
    print(f'[{member_id}] parsed {n_rxn} rxns + {len(ex)} exchanges | GPR grafted: {grafted}, '
          f'in source w/o rule: {no_rule}, absent from source: {len(missing)}')
    print(f'[{member_id}] no GPR source match: {sorted(missing)}')
    return member


def build_community(members, mdlid):
    msc = MSCommunity.build_from_species_models(
        members, mdlid=mdlid, name=mdlid, names=NAMES, abundances=dict(ABUNDANCES))
    msc.model.id = mdlid
    return msc


def apply_medium(model):
    for r in model.reactions:
        if r.id.startswith('EX_'):
            r.lower_bound = 0
            r.upper_bound = 1000
    for ex_id, v in GSP_MEDIUM.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -v


def species_rxns(model, idx):
    """Non-boundary reactions belonging to species index idx, keyed by base id."""
    out = {}
    for r in model.reactions:
        if r.id.startswith(('EX_', 'SK_', 'DM_')):
            continue
        m = re.match(r'^(.*_[a-z])(\d+)$', r.id) or re.match(r'^(.*\.)(\d+)$', r.id)
        if m and int(m.group(2)) == idx:
            out[m.group(1)] = r
    return out


def shared_rxns(model):
    """Boundary + purely-extracellular reactions (shared e0 pool)."""
    return {r.id: r for r in model.reactions
            if r.id.startswith(('EX_', 'SK_', 'DM_')) or r.id.endswith('_e0')}


def stoich_dict(r):
    return {m.id: round(v, 8) for m, v in r.metabolites.items()}


def compare(msc1, msc2, label1, label2):
    m1, m2 = msc1.model, msc2.model
    print('\n' + '=' * 70)
    print(f'STRUCTURE           {label1:>25} {label2:>25}')
    print(f'  reactions         {len(m1.reactions):>25} {len(m2.reactions):>25}')
    print(f'  metabolites       {len(m1.metabolites):>25} {len(m2.metabolites):>25}')
    print(f'  genes             {len(m1.genes):>25} {len(m2.genes):>25}')
    gpr1 = sum(1 for r in m1.reactions if r.gene_reaction_rule)
    gpr2 = sum(1 for r in m2.reactions if r.gene_reaction_rule)
    print(f'  rxns with GPR     {gpr1:>25} {gpr2:>25}')

    for idx, name in [(1, '3H11'), (2, 'R12')]:
        s1, s2 = species_rxns(m1, idx), species_rxns(m2, idx)
        only1, only2 = sorted(set(s1) - set(s2)), sorted(set(s2) - set(s1))
        both = set(s1) & set(s2)
        print(f'\n--- {name} (c{idx}) --- {label1}: {len(s1)} rxns | {label2}: {len(s2)} rxns')
        print(f'  only in {label1} ({len(only1)}):')
        for b in only1:
            print(f'    {b}{idx}  {s1[b].name}')
        print(f'  only in {label2} ({len(only2)}):')
        for b in only2:
            print(f'    {b}{idx}  {s2[b].name}')
        bd = [b for b in both if s1[b].bounds != s2[b].bounds]
        st = [b for b in both if stoich_dict(s1[b]) != stoich_dict(s2[b])]
        gp = [b for b in both
              if {g.id for g in s1[b].genes} != {g.id for g in s2[b].genes}]
        order_only = sum(1 for b in both
                         if {g.id for g in s1[b].genes} == {g.id for g in s2[b].genes}
                         and s1[b].gene_reaction_rule != s2[b].gene_reaction_rule)
        print(f'  shared {len(both)}: bounds differ {len(bd)}, stoichiometry differs {len(st)}, '
              f'gene set differs {len(gp)} (rule-order-only diffs ignored: {order_only})')
        for b in bd:
            print(f'    [bounds] {b}{idx}: {s1[b].bounds} vs {s2[b].bounds}')
        for b in st:
            print(f'    [stoich] {b}{idx}:')
            print(f'       {label1}: {s1[b].build_reaction_string()}')
            print(f'       {label2}: {s2[b].build_reaction_string()}')
        for b in gp:
            g1 = {g.id for g in s1[b].genes}
            g2 = {g.id for g in s2[b].genes}
            print(f'    [genes]  {b}{idx}: only {label1} {sorted(g1 - g2)} | only {label2} {sorted(g2 - g1)}')

    e1, e2 = shared_rxns(m1), shared_rxns(m2)
    only1, only2 = sorted(set(e1) - set(e2)), sorted(set(e2) - set(e1))
    print(f'\n--- shared e0 / boundary --- {label1}: {len(e1)} | {label2}: {len(e2)}')
    print(f'  only in {label1} ({len(only1)}): {only1}')
    print(f'  only in {label2} ({len(only2)}): {only2}')


def simulate(msc, label):
    model = msc.model
    model.solver = 'glpk'
    apply_medium(model)
    model.objective = 'bio1'
    try:
        sol = cobra.flux_analysis.pfba(model)
    except Exception:
        sol = model.optimize()
    print(f'\n### {label}: pFBA on GSP medium -> {sol.status}, community biomass = {sol.fluxes["bio1"]:.6f}')
    if sol.status != 'optimal':
        return
    watch = {
        'bio2': 'biomass 3H11', 'bio3': 'biomass R12',
        'dnr00001_c1': 'Nar 3H11', 'dnr00001_c2': 'Nar R12',
        'dnr00002_c1': 'Nir 3H11', 'dnr00002_c2': 'Nir R12',
        'dnr00003_c1': 'Nor 3H11', 'dnr00003_c2': 'Nor R12',
        'dnr00004_c1': 'Nos 3H11', 'dnr00004_c2': 'Nos R12',
        'LeuE_c1': 'LeuE 3H11', 'rxn05161_c2': 'Leu ABC R12',
        'EX_cpd00029_e0': 'acetate ex', 'EX_cpd00209_e0': 'nitrate ex',
        'EX_cpd00075_e0': 'nitrite ex', 'EX_cpd00418_e0': 'NO ex',
        'EX_cpd00659_e0': 'N2O ex', 'EX_cpd00528_e0': 'N2 ex',
        'EX_cpd00107_e0': 'leucine ex',
    }
    for rid, desc in watch.items():
        if rid in model.reactions:
            v = sol.fluxes[rid]
            if abs(v) > 1e-9:
                print(f'  {desc:<14} {rid:<16} {v: .4f}')
        else:
            print(f'  {desc:<14} {rid:<16} (absent)')
    # cooperativity check: cut leucine export from 3H11
    if 'LeuE_c1' in model.reactions:
        with model:
            model.reactions.LeuE_c1.knock_out()
            v = model.slim_optimize(error_value=float('nan'))
        print(f'  LeuE_c1 knockout -> max community biomass = {v:.6f}')


if __name__ == '__main__':
    print('== loading models ==')
    model_comm = cobra.io.load_json_model(f'{DATA}/comm_model.json')
    gpr_acido = cobra.io.load_json_model(f'{DATA}/acido_model.json')
    gpr_rhoda = cobra.io.load_json_model(f'{DATA}/rhoda_model.json')

    print('\n== version 1: members parsed from comm_model.json + grafted GPRs ==')
    member_a = parse_member_from_comm(model_comm, 'A', '3H11', gpr_acido)
    member_r = parse_member_from_comm(model_comm, 'R', 'R12', gpr_rhoda)
    msc_parsed = build_community([member_a, member_r], 'SynCom_from_comm')

    print('\n== version 2: members from repo isolate models ==')
    iso_a = cobra.io.load_json_model(f'{DATA}/model_acido_gf.json')
    iso_r = cobra.io.load_json_model(f'{DATA}/model_rhoda_gf.json')
    msc_isolates = build_community([iso_a, iso_r], 'SynCom_from_isolates')

    cobra.io.save_json_model(msc_parsed.model, f'{DATA}/comm_model_mscommunity_parsed.json')
    cobra.io.save_json_model(msc_isolates.model, f'{DATA}/comm_model_mscommunity_isolates.json')
    print('\nsaved data/comm_model_mscommunity_parsed.json and data/comm_model_mscommunity_isolates.json')

    compare(msc_parsed, msc_isolates, 'parsed-from-comm', 'from-isolates')
    simulate(msc_parsed, 'parsed-from-comm')
    simulate(msc_isolates, 'from-isolates')

    # attribute the growth gap between the two versions: knock out reaction
    # groups only present in the from-isolates members, and apply the comm
    # model's directionality curations, then re-optimize
    print('\n### growth-gap attribution (from-isolates, max community biomass)')
    m2 = msc_isolates.model
    groups = {
        'rxn14427 (Nar cyt-c +2H+)': ['rxn14427_c1', 'rxn14427_c2'],
        'rxn11937 (NAD N2O-forming)': ['rxn11937_c1', 'rxn11937_c2'],
        'SDH (rxn10126 + rxn00288)': ['rxn10126_c1', 'rxn10126_c2', 'rxn00288_c1', 'rxn00288_c2'],
        'rxn08792 (lactate oxidation)': ['rxn08792_c1'],
    }
    comm_bounds = [('rxn00083', (-1000, 1000)), ('rxn00154', (0, 1000)),
                   ('rxn00548', (0, 1000)), ('rxn00611', (-1000, 0))]
    for name, rids in groups.items():
        with m2:
            for rid in rids:
                if rid in m2.reactions:
                    m2.reactions.get_by_id(rid).knock_out()
            print(f'  KO {name}: {m2.slim_optimize(error_value=float("nan")):.6f}')
    with m2:
        for base, b in comm_bounds:
            for idx in (1, 2):
                rid = f'{base}_c{idx}'
                if rid in m2.reactions:
                    m2.reactions.get_by_id(rid).bounds = b
        print(f'  comm directionality bounds on rxn00083/154/548/611: '
              f'{m2.slim_optimize(error_value=float("nan")):.6f}')
