"""Fit the community's total (specific) acetate consumption from the data.

Uses the measured SynCom time course of analysis.ipynb
(CommPlots.get_exp_syncom): with the measured biomass curve X(t) integrated
by trapezoid, a fixed specific uptake rate q predicts

    acetate(t) = acetate(0) - q * integral_0^t X dtau ,

and q is fitted by least squares over the measured acetate series. The fit
is the data-derived community acetate uptake bound for the constrained
kinetic sweep (kinetic_sweep.py --fitted).

Run:  ~/Documents/py_venv/bin/python fit_acetate_uptake.py
"""
import json
from datetime import date

from plots import CommPlots

OUT = './data/fitted_acetate_uptake.json'

def fit_series(t, x, series):
    """Least-squares q for series(t) = series(0) - q * integral X dtau."""
    integral = [0.0]
    for k in range(1, len(t)):
        integral.append(integral[-1] + (x[k] + x[k - 1]) / 2 * (t[k] - t[k - 1]))
    drawdown = [series[0] - a for a in series]
    q = (sum(d * i for d, i in zip(drawdown, integral))
         / sum(i * i for i in integral))
    predicted = [series[0] - q * i for i in integral]
    ss_res = sum((p - a) ** 2 for p, a in zip(predicted, series))
    mean = sum(series) / len(series)
    r2 = 1 - ss_res / sum((a - mean) ** 2 for a in series)
    return q, r2, predicted


if __name__ == '__main__':
    exp = CommPlots.get_exp_syncom()
    t, x = exp['i'], exp['biomass']
    out = {'method': 'least-squares fit of C(t) = C(0) - q * trapezoid-integral of the '
                     'measured biomass curve, over the measured SynCom series '
                     '(CommPlots.get_exp_syncom, 0-119 h)',
           'date': date.today().isoformat()}
    for key, label in [('acetate', 'acetate'), ('no3', 'nitrate')]:
        q, r2, predicted = fit_series(t, x, exp[key])
        print(f'fitted specific {label} uptake q = {q:.2f} mmol/gDW/h (R^2 = {r2:.3f})')
        for tt, a, p in zip(t, exp[key], predicted):
            print(f'  t={tt:>4} h   measured {a:6.2f}   fitted {p:6.2f} mM')
        out[f'q_{label}_mmol_per_gDW_h'] = round(q, 3)
        out[f'{label}_r_squared'] = round(r2, 4)
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f'saved {OUT}')
