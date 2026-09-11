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

The same projection is also scored for direction alone. Each series is split
into its sampling intervals and every interval classed as a rise, a fall or
flat; the directional error is the fraction of intervals whose class the
projection misses, averaged over the same five series. DEAD_BAND sets "flat"
for the measurements only -- it is there to absorb measurement noise, which the
projection does not have, so predictions are classed by exact sign. Any band
from 1.7 to 3.3 % of a series' range classes the measured series identically:
the largest noise step is 1.6 % (N2O after its collapse) and the smallest real
change 3.4 % (acetate's first interval).

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
DEAD_BAND = 0.025  # fraction of a measured series' range within which a change is noise
EXACT = 1e-9       # the projection is noise-free: only a numerically zero change is flat


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


def predictions(run, exp):
    """Each series projected to the experimental sampling times."""
    mu = run['community_biomass'] or 0.0
    traj = project(mu, community_rates(run), exp, max(exp['i']))
    return {key: [traj[float(t)][0] if key == 'biomass' else traj[float(t)][1][key]
                  for t in exp['i']]
            for key in list(SERIES) + ['biomass']}


def nrmse(pred, exp):
    errors = {}
    for key, predicted in pred.items():
        lo, hi = min(exp[key]), max(exp[key])
        span = (hi - lo) or 1.0
        sq = 0.0
        for p, measured in zip(predicted, exp[key]):
            sq += ((p - measured) / span) ** 2
        errors[key] = round(math.sqrt(sq / len(exp['i'])), 4)
    errors['mean'] = round(sum(errors.values()) / len(errors), 4)
    return errors


def direction(delta, band):
    return 0 if abs(delta) <= band else (1 if delta > 0 else -1)


def directional_error(pred, exp):
    """Fraction of sampling intervals whose rise/fall/flat class the projection misses."""
    errors = {}
    for key, predicted in pred.items():
        measured = exp[key]
        span = (max(measured) - min(measured)) or 1.0
        wrong = sum(direction(measured[i + 1] - measured[i], DEAD_BAND * span)
                    != direction(predicted[i + 1] - predicted[i], EXACT * span)
                    for i in range(len(measured) - 1))
        errors[key] = round(wrong / (len(measured) - 1), 4)
    errors['mean'] = round(sum(errors.values()) / len(errors), 4)
    return errors


if __name__ == '__main__':
    exp = CommPlots.get_exp_syncom()
    # The projection integrates dC/dt = q*X with q in mmol/gDW/h, so X must be gDW/L --
    # the same conversion fit_acetate_uptake.py applies to derive q. Scoring the raw
    # OD_coeff series against a gDW/L trajectory makes the biomass NRMSE meaningless.
    scale = json.load(open('./data/fitted_acetate_uptake.json'))['biomass_scale_factor']
    exp['biomass'] = [b * scale for b in exp['biomass']]
    runs = json.load(open(SRC))['runs']
    fits, directional = {}, {}
    print(f'{"K":>6} {"NRMSE":>7} {"dir":>6}   (NRMSE per series: biomass, acetate, no3, no2, n2o)')
    for run in runs:
        pred = predictions(run, exp)
        e, d = nrmse(pred, exp), directional_error(pred, exp)
        fits[str(run['kinetic_coeff'])] = e
        directional[str(run['kinetic_coeff'])] = d
        print(f"{run['kinetic_coeff']:>6} {e['mean']:>7.3f} {d['mean']:>6.3f}   "
              f"({e['biomass']:.3f}, {e['acetate']:.3f}, {e['no3']:.3f}, {e['no2']:.3f}, {e['n2o']:.3f})")
    # a community that never grows: every series holds its initial measurement
    dead = predictions({'community_biomass': 0.0, 'members': {}}, exp)
    baseline = {'NRMSE': nrmse(dead, exp)['mean'],
                'directional_error': directional_error(dead, exp)['mean']}
    print(f"{'dead':>6} {baseline['NRMSE']:>7.3f} {baseline['directional_error']:>6.3f}   "
          f"(flat no-growth baseline)")
    with open(OUT, 'w') as fh:
        json.dump({'_metadata': {
            'source': SRC,
            'experiment': 'CommPlots.get_exp_syncom() (plots.py), the measured SynCom '
                          'time course of analysis.ipynb, 0-119 h',
            'method': 'hourly Euler projection of each solution\'s community specific '
                      'rates with growing biomass and proportional substrate-exhaustion '
                      'scaling; NRMSE per measured series (normalized by its range), '
                      'averaged over biomass, acetate, NO3, NO2, N2O; lower is better. '
                      'Directional error over the same projection: the fraction of '
                      'sampling intervals whose rise/fall/flat class is missed, measured '
                      'changes within dead_band_fraction_of_range counted flat and '
                      'predictions classed by exact sign, averaged over the same series; '
                      'lower is better',
            'dead_band_fraction_of_range': DEAD_BAND,
            'no_growth_baseline': baseline,
            'date': date.today().isoformat(),
        }, 'fits': fits, 'directional': directional}, fh, indent=1)
    print(f'saved {OUT}')
