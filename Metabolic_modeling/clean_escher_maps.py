"""Clean the SynCom Escher maps: chemical names on the labels, no flux numbers.

Applies the editEscher cleaning stage (escher_edit.clean_json.cleanEscherJSON)
to every map make_escher_interactions.py generates:

* MSID -> chemical name: Escher renders a node's bigg_id as the on-map label,
  so the maps display raw ModelSEED ids (cpd00075, ...). The abbreviation map
  handed to cleanEscherJSON is built from the local ModelSEED Database
  (modelseedpy.biochem.from_local on ../../ModelSEEDDatabase), following
  clean_json.build_name_abbrev_table's pattern, so every metabolite label
  becomes the database's chemical name (cpd00075 -> Nitrite).
* Reaction labels are reduced to the bare member name (3H11k1285 -> 3H11)
  through cleanEscherJSON's rxn_abbrev_map; the kappa condition stays in the
  block caption / file name.
* Quantitative flux labels are removed: the builder writes each compound's
  net flux as its stoichiometric coefficient, and Escher draws every non-unit
  coefficient as a number beside the metabolite (the "stoichiometry-labels"
  group that svg_editor.remove_label_groups strips from exported SVGs). At
  the JSON stage the equivalent is normalizing coefficients to +/-1 -
  direction is kept, the numbers disappear.

cleanEscherJSON writes *_cleaned0.json siblings; the raw maps keep the
quantitative fluxes for analysis.

Run:  ~/Documents/py_venv/bin/python clean_escher_maps.py
"""
import glob
import json
import warnings

warnings.filterwarnings('ignore')

from escher_edit import cleanEscherJSON
from modelseedpy.biochem import from_local

MSDB_PATH = '/home/freiburger/Documents/ModelSEEDDatabase'
MAPS = sorted(glob.glob('./data/escher_syncom_membermap*.json'))
MAPS = [m for m in MAPS if '_cleaned' not in m]

msdb = from_local(MSDB_PATH)

# collect every compound and reaction id used across the maps
cpds, rxns = set(), set()
for path in MAPS:
    body = json.load(open(path))[1]
    cpds |= {n['bigg_id'] for n in body['nodes'].values()
             if n.get('node_type') == 'metabolite'}
    rxns |= set(r['bigg_id'] for r in body['reactions'].values())

abbrev_map = {cpd: msdb.compounds.get_by_id(cpd).name for cpd in sorted(cpds)}
rxn_abbrev_map = {r: r.split('k')[0] for r in rxns}   # 3H11k1285 -> 3H11
print(f'{len(abbrev_map)} compounds named from {MSDB_PATH}, '
      f'{len(rxn_abbrev_map)} reaction labels reduced')

for path in MAPS:
    out = cleanEscherJSON(path, abbrev_map=abbrev_map,
                          rxn_abbrev_map=rxn_abbrev_map, rewrite_bigg_id=True)
    # remove the quantitative flux labels: Escher only renders stoichiometry
    # numbers for non-unit coefficients
    cleaned = json.load(open(out))
    for r in cleaned[1]['reactions'].values():
        for met in r.get('metabolites', []):
            met['coefficient'] = 1 if met['coefficient'] > 0 else -1
    json.dump(cleaned, open(out, 'w'), indent=1)
    print(f'cleaned {out}')
