"""Score each kinetic-sweep solution against the measured SynCom time course.

analysis.ipynb judges the fit of the model's delegation to the data visually:
CommPlots.get_exp_syncom() supplies the measured community time course
(biomass, acetate, NO3, NO2, N2O at 0-119 h) and generate_total_acc_data
integrates simulated fluxes hourly for plot_total_acc to overlay on the
measurements. This script quantifies that same comparison for every kappa of
a sweep: each solution's community-level specific exchange rates
(mmol/gDW/h, summed over both members) are forward-integrated hourly with
growing biomass,

    dX/dt = mu * X ,   dC_i/dt = q_i * X ,

with all fluxes scaled down proportionally in any hour a consumed substrate
would be overdrawn (the flux vector is fixed, so proportional scaling
preserves its stoichiometry). Predictions are compared to the measurements
at the experimental time points as a normalized RMSE per series (RMSE over
the series' measured range) and averaged across the five measured series;
lower is better, and a run with no growth scores the flat "dead community"
baseline.

Results: data/kinetic_sweep_fit_<suffix>.json, which make_sweep_summary.py
merges into the summary CSV as a fit column.

Run:  ~/Documents/py_venv/bin/python evaluate_syncom_fit.py \
          [data/kinetic_sweep_net_reactions_800-1400.json]
"""
import json
import math
import sys
from datetime import date

from plots import CommPlots

SRC = sys.argv[1] if len(sys.argv) > 1 else './data/kinetic_sweep_net_reactions_800-1400.json'
OUT = SRC.replace('kinetic_sweep_net_reactions', 'kinetic_sweep_fit')
SERIES = {'acetate': 'Acetate', 'no3': 'Nitrate', 'no2': 'Nitrite', 'n2o': 'Nitrous oxide'}
DT = 1.0  # h, the hourly stepping of generate_total_acc_data


def community_rates(run):
    """Community-level specific exchange rates, summed over both members."""
    rates = {}
    for member in run['members'].values():
        for name, v in member.get('net_exchange_mmol_gDW_h', {}).items():
            rates[name] = rates.get(name, 0.0) + v
    return rates


def project(mu, rates, exp, horizon):
    """Hourly Euler projection from the measured initial state."""
    x = exp['biomass'][0]
    conc = {key: exp[key][0] for key in SERIES}
    q = {key: rates.get(name, 0.0) for key, name in SERIES.items()}
    t = 0.0
    trajectory = {0.0: (x, dict(conc))}
    while t < horizon:
        factor = 1.0
        for key, qi in q.items():
            drawn = -qi * x * DT
            if drawn > 0 and drawn > conc[key]:
                factor = min(factor, conc[key] / drawn)
        for key, qi in q.items():
            conc[key] = max(conc[key] + qi * x * DT * factor, 0.0)
        x += mu * x * DT * factor
        t += DT
        trajectory[t] = (x, dict(conc))
    return trajectory


def nrmse(run, exp):
    mu = run['community_biomass'] or 0.0
    rates = community_rates(run)
    traj = project(mu, rates, exp, max(exp['i']))
    errors = {}
    for key in list(SERIES) + ['biomass']:
        lo, hi = min(exp[key]), max(exp[key])
        span = (hi - lo) or 1.0
        sq = 0.0
        for t, measured in zip(exp['i'], exp[key]):
            x, conc = traj[float(t)]
            predicted = x if key == 'biomass' else conc[key]
            sq += ((predicted - measured) / span) ** 2
        errors[key] = round(math.sqrt(sq / len(exp['i'])), 4)
    errors['mean'] = round(sum(errors.values()) / len(errors), 4)
    return errors


if __name__ == '__main__':
    exp = CommPlots.get_exp_syncom()
    runs = json.load(open(SRC))['runs']
    fits = {}
    print(f'{"K":>6} {"NRMSE":>7}   (per series: biomass, acetate, no3, no2, n2o)')
    for run in runs:
        e = nrmse(run, exp)
        fits[str(run['kinetic_coeff'])] = e
        print(f"{run['kinetic_coeff']:>6} {e['mean']:>7.3f}   "
              f"({e['biomass']:.3f}, {e['acetate']:.3f}, {e['no3']:.3f}, {e['no2']:.3f}, {e['n2o']:.3f})")
    with open(OUT, 'w') as fh:
        json.dump({'_metadata': {
            'source': SRC,
            'experiment': 'CommPlots.get_exp_syncom() (plots.py), the measured SynCom '
                          'time course of analysis.ipynb, 0-119 h',
            'method': 'hourly Euler projection of each solution\'s community specific '
                      'rates with growing biomass and proportional substrate-exhaustion '
                      'scaling; NRMSE per measured series (normalized by its range), '
                      'averaged over biomass, acetate, NO3, NO2, N2O; lower is better',
            'date': date.today().isoformat(),
        }, 'fits': fits}, fh, indent=1)
    print(f'saved {OUT}')
