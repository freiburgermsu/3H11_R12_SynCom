"""Build a community-interaction Escher map from the fitted kinetic sweep.

Applies the editEscher method (escher_edit.build_map, installed from
~/Documents/editEscher): each member of each condition becomes one net
organism reaction whose reactants are the compounds the member consumes and
whose products are the compounds it excretes, drawn as a self-contained
cluster.

This script closes the input gap: it converts the per-member net exchange
fluxes of the fitted sweep (data/kinetic_sweep_net_reactions_fitted_*.json,
raw mmol/gDW/h, consumption negative / secretion positive - already the
builder's convention) into the wide ASVMetaboliteInteractions-style matrix
the builder parses (rows = members; columns = <cpd>_<condition> with
ModelSEED compound ids recovered from the canonical community model), writes
the parallel names matrix for display labels, and generates the map. The
degenerate post-saturation rungs (kappa > 2500) are excluded.

Outputs (data/):
  escher_syncom_interactions.csv        wide flux matrix, cpd-id headers
  escher_syncom_interactions_names.csv  same matrix, display-name headers
  escher_syncom_membermap.json          the Escher map

Run:  ~/Documents/py_venv/bin/python make_escher_interactions.py \
          [data/kinetic_sweep_net_reactions_fitted_250-3000.json]
"""
import csv
import json
import sys

from escher_edit.build_map import build_map_from_interactions

SRC = sys.argv[1] if len(sys.argv) > 1 else './data/kinetic_sweep_net_reactions_fitted_250-3000.json'
MODEL = './data/comm_model_mscommunity_parsed.json'
CSV_OUT = './data/escher_syncom_interactions.csv'
NAMES_OUT = './data/escher_syncom_interactions_names.csv'
MAP_OUT = './data/escher_syncom_membermap.json'
MIN_FLUX = 0.01          # drops the ~4e-4 trace-metal exchanges
MAX_KAPPA = 2500         # rungs above saturation are degenerate duplicates
MEMBERS = ['3H11', 'R12']

# display name (as recorded in the sweep JSON) -> ModelSEED compound id,
# recovered from the canonical model's extracellular metabolites
name_to_cpd = {}
for met in json.load(open(MODEL))['metabolites']:
    if met['id'].endswith('_e0'):
        name_to_cpd[met['name'].replace(' [e0]', '')] = met['id'][:-3]

runs = json.load(open(SRC))['runs']
conditions, fluxes = [], {}   # fluxes[(member, cpd, condition)] = value
compounds, cpd_name = [], {}
for run in runs:
    if run['kinetic_coeff'] > MAX_KAPPA:
        continue
    members = {m: e.get('net_exchange_mmol_gDW_h') for m, e in run['members'].items()}
    if not any(members.values()):
        continue                       # infeasible rung
    cond = f"k{run['kinetic_coeff']}"
    conditions.append(cond)
    for member, net in members.items():
        for name, value in (net or {}).items():
            cpd = name_to_cpd[name]
            if cpd not in cpd_name:
                compounds.append(cpd)
                cpd_name[cpd] = name
            fluxes[(member, cpd, cond)] = value

id_header = [''] + [f'{cpd}_{cond}' for cond in conditions for cpd in compounds]
name_header = [''] + [cpd_name[cpd] for _ in conditions for cpd in compounds]
for path, header in [(CSV_OUT, id_header), (NAMES_OUT, name_header)]:
    with open(path, 'w', newline='') as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for member in MEMBERS:
            writer.writerow([member] + [
                fluxes.get((member, cpd, cond), 0)
                for cond in conditions for cpd in compounds])
print(f'wrote {CSV_OUT} and {NAMES_OUT}: {len(MEMBERS)} members x '
      f'{len(conditions)} conditions x {len(compounds)} compounds')

build_map_from_interactions(
    CSV_OUT,
    output_path=MAP_OUT,
    names_csv=NAMES_OUT,
    min_abs_flux=MIN_FLUX,
)

# each kappa community fully isolated: one map per condition, so no nodes,
# layout, or Escher same-species identity are shared between simulations
paths = build_map_from_interactions(
    CSV_OUT,
    output_path=MAP_OUT,
    names_csv=NAMES_OUT,
    min_abs_flux=MIN_FLUX,
    separate_maps=True,
)
print('per-condition maps:', ', '.join(p.name for p in paths))
