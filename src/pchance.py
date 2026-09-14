"""Null sample -> chance probability, shared by every scoring script.

`p_chance` is always the same object in this project: the fraction of null
positions whose best association_TS is at least as high as the observed one,

    p(t) = P(null best_TS >= t)

read off the EMPIRICAL survival function of the null sample and interpolated
onto a 4000-point grid — the same grid size the official
`random_gaia_ts_prob.txt` mapping uses, so our tables are line-comparable with
it. Outside the sampled range the interpolator returns 1 below and 0 above:
a p of exactly 0 does not mean "impossible", it means "no null draw reached
this TS", i.e. p < 1/len(null). Report the null size alongside it.

Callers: scripts/wise_pchance/score_wise.py, scripts/racs_pchance/score_racs.py,
scripts/null_resample/{rescore_floors,make_galactic_mappings}.py — which each
carried a byte-identical private copy of this until 2026-09.
"""
import numpy as np
from scipy.interpolate import interp1d

GRID_PTS = 4000


def survival_table(ts, n_grid=GRID_PTS):
    """(grid, p) — the empirical survival function of a null TS sample.

    Non-finite entries (empty background cones) are dropped, not counted as
    zeros: a position with no catalog source in its cone has no best match,
    so it is outside the population the p-value is about.
    """
    ts = np.asarray(ts, float)
    ts = ts[np.isfinite(ts)]
    if ts.size == 0:
        raise ValueError('null sample has no finite entries')
    grid = np.linspace(ts.min(), ts.max(), n_grid)
    p = 1.0 - np.searchsorted(np.sort(ts), grid, side='right') / ts.size
    return grid, p


def survival_interp(ts, n_grid=GRID_PTS):
    """Cubic interpolator over `survival_table`, filling (1, 0) outside.

    Matches get_gaia_prob's own interpolation, so a TS scored here and a TS
    scored by the official code against the same mapping give the same p.
    """
    grid, p = survival_table(ts, n_grid=n_grid)
    return interp1d(grid, p, kind='cubic', fill_value=(1, 0),
                    bounds_error=False)
