"""best association_TS for many (sigma_ra, sigma_dec) pairs in one pass.

Same math as scripts/null_resample/null_ts.py (which replicates
spt3g.sources.gaia.get_gaia_association), but sigma is supplied directly
instead of being derived from an event TS, and every sigma pair is scored
against the same cone so the expensive 0.12-deg neighbour query is done once.

The spherical geometry (ARCSEC/DEG, unit vectors, chord length, Vincenty
separation) is imported from null_ts rather than re-implemented: null_ts is
the copy that was verified against spt3g's get_gaia_association to 2e-11, so
it stays the single definition and this module cannot silently drift from it.
"""
import os
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "null_resample"))
from null_ts import (ARCSEC, DEG, _chord, _unit_vectors,  # noqa: F401
                     vincenty_sep)


class MultiSigmaEngine:
    def __init__(self, ra_deg, dec_deg, gmag):
        self.ra = np.deg2rad(np.asarray(ra_deg, float))
        self.dec = np.deg2rad(np.asarray(dec_deg, float))
        self.gmag = np.asarray(gmag, float)
        self.tree = cKDTree(_unit_vectors(self.ra, self.dec))

    def best_ts(self, ev_ra_deg, ev_dec_deg, sigma_pairs_arcsec,
                bg_radius=0.12 * DEG, cand_radius=60 * ARCSEC,
                guarantee_radius=40 * ARCSEC, batch=1024, workers=-1,
                progress=None, exact=False):
        """(n_events, n_sigma) array of best association_TS.

        exact=True scans the whole background cone for every event (used for
        the handful of real sources); otherwise the same candidate-pruning
        guarantee as null_ts.py applies -- valid a fortiori here because these
        sigmas are smaller than the ones it was derived for, which only steepens
        the chi2 penalty on distant stars.
        """
        sig = np.asarray(sigma_pairs_arcsec, float) * ARCSEC
        s_ra = sig[:, 0][:, None]
        s_dec = sig[:, 1][:, None]
        ev_ra = np.deg2rad(np.asarray(ev_ra_deg, float))
        ev_dec = np.deg2rad(np.asarray(ev_dec_deg, float))
        n, k = len(ev_ra), len(sig)
        out = np.full((n, k), np.nan)
        ev_vec = _unit_vectors(ev_ra, ev_dec)

        for lo in range(0, n, batch):
            hi = min(lo + batch, n)
            cones = self.tree.query_ball_point(ev_vec[lo:hi], _chord(bg_radius),
                                               workers=workers)
            if not exact:
                cands = self.tree.query_ball_point(
                    ev_vec[lo:hi], _chord(cand_radius), workers=workers)
            for j in range(hi - lo):
                i = lo + j
                cone_idx = np.asarray(cones[j], dtype=np.intp)
                if cone_idx.size == 0:
                    continue
                if exact:
                    cand_idx = cone_idx
                else:
                    cand_idx = np.asarray(cands[j], dtype=np.intp)
                    use_full = cand_idx.size == 0
                    if not use_full:
                        d = vincenty_sep(ev_ra[i], ev_dec[i],
                                         self.ra[cand_idx], self.dec[cand_idx])
                        use_full = d.min() > guarantee_radius
                    if use_full:
                        cand_idx = cone_idx
                sorted_mags = np.sort(self.gmag[cone_idx])
                rank = sorted_mags.searchsorted(self.gmag[cand_idx], "left") + 1
                ra_dist = vincenty_sep(ev_ra[i], ev_dec[i],
                                       self.ra[cand_idx], ev_dec[i])
                dec_dist = np.abs(self.dec[cand_idx] - ev_dec[i])
                base = -np.log(rank) + np.log(4 * np.pi) + np.log(cone_idx.size)
                chi2 = (ra_dist[None, :] / s_ra) ** 2 + \
                       (dec_dist[None, :] / s_dec) ** 2
                out[i] = (-0.5 * chi2 + base[None, :]).max(axis=1)
            if progress is not None:
                progress(hi, n)
        return out
