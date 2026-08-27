"""Attach ModelSEED reaction thermodynamics to the parsed MSCommunity model.

Pulls the reaction table served by https://modelseed.org/biochem/reactions
(the site's Solr backend) and joins it against every reaction in
data/comm_model_mscommunity_parsed.json, keyed by the model reaction id.
Values: deltag / deltagerr (kcal/mol, pH 7), status, reversibility,
is_transport. The sentinel 10000000 (unknown ΔG) is stored as null.

Run:  ~/Documents/py_venv/bin/python parse_thermodynamics.py [cached_solr.json]
"""
import json
import re
import sys
import urllib.request
from datetime import date

MODEL = './data/comm_model_mscommunity_parsed.json'
OUT = './data/reaction_thermodynamics.json'
SOLR_URL = ('https://modelseed.org/solr/reactions/select?wt=json&q=*:*&rows=100000'
            '&fl=id,name,deltag,deltagerr,status,reversibility,direction,is_transport')
UNKNOWN = 1e6  # ModelSEED marks unknown deltag/deltagerr as 10000000

# model reactions that are curated additions with an exact ModelSEED equivalent
EXPLICIT_SEED = {'ATPM': 'rxn00062'}  # ATP hydrolysis


def load_seed_reactions(cached=None):
    if cached:
        with open(cached) as fh:
            data = json.load(fh)
    else:
        with urllib.request.urlopen(SOLR_URL, timeout=300) as fh:
            data = json.load(fh)
    docs = data['response']['docs']
    print(f'ModelSEED reactions fetched: {len(docs)}')
    return {d['id']: d for d in docs}


def clean(value):
    if value is None or abs(value) >= UNKNOWN:
        return None
    return value


if __name__ == '__main__':
    seed = load_seed_reactions(sys.argv[1] if len(sys.argv) > 1 else None)
    with open(MODEL) as fh:
        model = json.load(fh)

    entries = {}
    n_hit = n_custom = n_boundary = n_miss = 0
    for r in model['reactions']:
        rid = r['id']
        if rid.startswith(('EX_', 'SK_', 'DM_')) or rid.startswith('bio'):
            entries[rid] = {'seed_id': None, 'note': 'boundary/biomass reaction',
                            'deltag': None, 'deltagerr': None, 'status': None,
                            'reversibility': None, 'is_transport': None}
            n_boundary += 1
            continue
        m = re.match(r'^(rxn\d+)', rid)
        base = m.group(1) if m else EXPLICIT_SEED.get(rid.split('_')[0])
        doc = seed.get(base) if base else None
        if doc:
            entries[rid] = {
                'seed_id': base,
                'seed_name': doc.get('name'),
                'deltag': clean(doc.get('deltag')),
                'deltagerr': clean(doc.get('deltagerr')),
                'status': doc.get('status'),
                'reversibility': doc.get('reversibility'),
                'is_transport': doc.get('is_transport'),
            }
            n_hit += 1
        else:
            entries[rid] = {'seed_id': base, 'note': 'no ModelSEED biochemistry entry',
                            'deltag': None, 'deltagerr': None, 'status': None,
                            'reversibility': None, 'is_transport': None}
            if m:
                n_miss += 1
            else:
                entries[rid]['note'] = 'custom curated reaction (not in ModelSEED)'
                n_custom += 1

    with_dg = sum(1 for e in entries.values() if e['deltag'] is not None)
    out = {
        '_metadata': {
            'model': MODEL,
            'source': 'https://modelseed.org/biochem/reactions (Solr backend)',
            'retrieved': date.today().isoformat(),
            'units': 'kcal/mol (ModelSEED group-contribution estimates, pH 7)',
            'notes': 'deltag/deltagerr of null = unknown in ModelSEED (sentinel 10000000); '
                     'ATPM_* mapped to rxn00062 (ATP hydrolysis)',
            'reactions_total': len(entries),
            'matched_to_modelseed': n_hit,
            'with_deltag': with_dg,
            'boundary_or_biomass': n_boundary,
            'custom_reactions': n_custom,
            'seed_id_without_entry': n_miss,
        },
        'reactions': entries,
    }
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out['_metadata'], indent=2))
    unmatched = sorted({e['seed_id'] or rid for rid, e in entries.items()
                        if e['deltag'] is None and not rid.startswith(('EX_', 'SK_', 'bio'))})
    print(f'reactions without deltag ({len(unmatched)} unique):')
    print(unmatched)
