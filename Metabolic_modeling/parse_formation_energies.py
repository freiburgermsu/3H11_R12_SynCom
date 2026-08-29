"""Attach ModelSEED compound formation energies to the parsed MSCommunity model.

Reads the compound tables of a local ModelSEEDDatabase checkout (default:
../ModelSEEDDatabase relative to the repo root, i.e. sibling of the repo) and
joins them against every metabolite in data/comm_model_mscommunity_parsed.json,
keyed by the model metabolite id. Values: deltag / deltagerr of formation
(kcal/mol, group-contribution / eQuilibrator estimates at pH 7), formula,
charge, mass. The sentinel 10000000 (unknown) is stored as null.

Run:  ~/Documents/py_venv/bin/python parse_formation_energies.py [db_path]
"""
import glob
import json
import re
import sys
from datetime import date
from pathlib import Path

MODEL = './data/comm_model_mscommunity_parsed.json'
OUT = './data/compound_formation_energies.json'
UNKNOWN = 1e6  # ModelSEED marks unknown deltag/deltagerr as 10000000


def load_seed_compounds(db_path):
    files = sorted(glob.glob(str(Path(db_path) / 'Biochemistry' / 'compound_*.json')))
    if not files:
        raise SystemExit(f'no compound_*.json files under {db_path}/Biochemistry')
    compounds = {}
    for f in files:
        with open(f) as fh:
            for c in json.load(fh):
                compounds[c['id']] = c
    print(f'ModelSEED compounds loaded: {len(compounds)} from {len(files)} files')
    return compounds


def clean(value):
    if value is None or abs(value) >= UNKNOWN:
        return None
    return value


if __name__ == '__main__':
    repo_root = Path(__file__).resolve().parents[1]
    db_path = sys.argv[1] if len(sys.argv) > 1 else repo_root.parent / 'ModelSEEDDatabase'
    seed = load_seed_compounds(db_path)

    with open(MODEL) as fh:
        model = json.load(fh)

    entries = {}
    n_hit = n_miss = n_custom = 0
    for met in model['metabolites']:
        mid = met['id']
        m = re.match(r'^(cpd\d+)', mid)
        base = m.group(1) if m else None
        doc = seed.get(base) if base else None
        if doc:
            entries[mid] = {
                'seed_id': base,
                'seed_name': doc.get('name'),
                'formula': doc.get('formula'),
                'charge': doc.get('charge'),
                'mass': clean(doc.get('mass')),
                'deltag': clean(doc.get('deltag')),
                'deltagerr': clean(doc.get('deltagerr')),
                'is_obsolete': bool(doc.get('is_obsolete')),
            }
            n_hit += 1
        elif base:
            entries[mid] = {'seed_id': base, 'note': 'no ModelSEED biochemistry entry',
                            'deltag': None, 'deltagerr': None}
            n_miss += 1
        else:
            entries[mid] = {'seed_id': None, 'note': 'custom metabolite (not in ModelSEED)',
                            'deltag': None, 'deltagerr': None}
            n_custom += 1

    with_dg = sum(1 for e in entries.values() if e['deltag'] is not None)
    out = {
        '_metadata': {
            'model': MODEL,
            'source': f'local ModelSEEDDatabase checkout ({Path(db_path).resolve()})',
            'retrieved': date.today().isoformat(),
            'units': 'kcal/mol formation energy (ModelSEED estimates, pH 7)',
            'notes': 'deltag/deltagerr of null = unknown in ModelSEED (sentinel 10000000)',
            'metabolites_total': len(entries),
            'matched_to_modelseed': n_hit,
            'with_deltag': with_dg,
            'seed_id_without_entry': n_miss,
            'custom_metabolites': n_custom,
        },
        'compounds': entries,
    }
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out['_metadata'], indent=2))
    no_dg = sorted({e['seed_id'] or mid for mid, e in entries.items() if e['deltag'] is None})
    print(f'compounds without deltag ({len(no_dg)} unique):')
    print(no_dg)
