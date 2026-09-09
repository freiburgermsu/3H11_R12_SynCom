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

# The measured `biomass` series of CommPlots.get_exp_syncom is OD600 * OD_coeff, with
# OD_coeff = 0.006 (3H11) / 0.008 (R12) from analysis.ipynb cell 4. Those values are not
# gDW/L: the conventional OD600 -> dry weight conversion is 0.30-0.50 gDW/L, 40-80x
# larger. Carr et al. 2025 report 10-ml Balch tubes, so 0.006 g per OD unit in a 10-ml
# culture is 0.6 gDW/L per OD -- the right order, i.e. the coefficient reads as grams in
# the vessel rather than grams per litre, though 0.6 still sits above every published
# calibration.
#
# The series has to be gDW/L to be dimensionally consistent with the mM concentrations it
# is divided into; otherwise q = d[C]/integral(X dt) inherits a factor of 1/V and comes
# out ~50x too large (44.5 rather than ~0.8 mmol/gDW/h), which is what forced the FBA
# model to grow 35x faster than the measured curve.
#
# 0.35 is the conventional mid-range value and is what the data supports: it puts the
# observed yield at 0.376 gDW per g acetate against the model's own 0.362 (a 4% match).
# Alternatives: 0.30 -> 0.322 g/g, 0.50 -> 0.537 g/g, 0.60 (pure 1/V) -> 0.644 g/g.
PUBLISHED_OD_COEFF = 0.006      # g per OD600 unit, analysis.ipynb cell 4
GDW_PER_L_PER_OD600 = 0.35      # conventional OD600 -> dry weight conversion
BIOMASS_SCALE = GDW_PER_L_PER_OD600 / PUBLISHED_OD_COEFF

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
    # convert the measured series from the published OD_coeff units to gDW/L
    t, x = exp['i'], [b * BIOMASS_SCALE for b in exp['biomass']]
    out = {'method': 'least-squares fit of C(t) = C(0) - q * trapezoid-integral of the '
                     'measured biomass curve, over the measured SynCom series '
                     '(CommPlots.get_exp_syncom, 0-119 h), with the biomass series '
                     'converted from OD600 * 0.006 g to gDW/L at '
                     f'{GDW_PER_L_PER_OD600} gDW/L per OD600',
           'gdw_per_L_per_OD600': GDW_PER_L_PER_OD600,
           'published_OD_coeff': PUBLISHED_OD_COEFF,
           'biomass_scale_factor': BIOMASS_SCALE,
           'date': date.today().isoformat()}
    for key, label in [('acetate', 'acetate'), ('no3', 'nitrate')]:
        q, r2, predicted = fit_series(t, x, exp[key])
        print(f'fitted specific {label} uptake q = {q:.2f} mmol/gDW/h (R^2 = {r2:.3f})')
        for tt, a, p in zip(t, exp[key], predicted):
            print(f'  t={tt:>4} h   measured {a:6.2f}   fitted {p:6.2f} mM')
        # 6 dp, not 3: after the OD -> gDW/L conversion q is O(1) rather than O(50), so
        # 3 dp is a ~5e-4 relative truncation -- enough to move the feasibility threshold
        out[f'q_{label}_mmol_per_gDW_h'] = round(q, 6)
        out[f'{label}_r_squared'] = round(r2, 4)
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print(f'saved {OUT}')
