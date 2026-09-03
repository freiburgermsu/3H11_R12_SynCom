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

if __name__ == '__main__':
    exp = CommPlots.get_exp_syncom()
    t, x, ac = exp['i'], exp['biomass'], exp['acetate']
    # cumulative integral of measured biomass (trapezoid), gDW h / L
    integral = [0.0]
    for k in range(1, len(t)):
        integral.append(integral[-1] + (x[k] + x[k - 1]) / 2 * (t[k] - t[k - 1]))
    drawdown = [ac[0] - a for a in ac]
    q = (sum(d * i for d, i in zip(drawdown, integral))
         / sum(i * i for i in integral))
    predicted = [ac[0] - q * i for i in integral]
    ss_res = sum((p - a) ** 2 for p, a in zip(predicted, ac))
    mean = sum(ac) / len(ac)
    ss_tot = sum((a - mean) ** 2 for a in ac)
    r2 = 1 - ss_res / ss_tot
    print(f'fitted specific acetate uptake q = {q:.2f} mmol/gDW/h (R^2 = {r2:.3f})')
    for tt, a, p in zip(t, ac, predicted):
        print(f'  t={tt:>4} h   measured {a:6.2f}   fitted {p:6.2f} mM')
    with open(OUT, 'w') as fh:
        json.dump({'q_acetate_mmol_per_gDW_h': round(q, 3), 'r_squared': round(r2, 4),
                   'method': 'least-squares fit of acetate(t) = acetate(0) - q * '
                             'trapezoid-integral of the measured biomass curve, over the '
                             'measured SynCom acetate series (CommPlots.get_exp_syncom, 0-119 h)',
                   'date': date.today().isoformat()}, fh, indent=1)
    print(f'saved {OUT}')
