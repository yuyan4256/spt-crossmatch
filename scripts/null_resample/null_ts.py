"""Vectorized replication of spt3g.sources.gaia.get_gaia_association best-match
association_TS, for large batches of (random) event positions.

Exact math (G3Units.rad == 1, so everything below is plain radians):
  bg cone   : great-circle dist <= 1 deg (the function's dec/RA pre-masks are
              supersets of this final cut for |dec| <= 72 deg, so the final
              star set is identical)
  rank      : num_brighter = searchsorted(sorted cone mags, mag, 'left') + 1
  num_total : stars in cone
  sigma_ra  : sqrt(floor_ra^2  + (50 arcsec / sqrt(event_ts))^2)
  sigma_dec : sqrt(floor_dec^2 + (49 arcsec / sqrt(event_ts))^2)
  ra_dist   : separation((ev_ra, ev_dec), (star_ra, ev_dec))  [Vincenty]
  dec_dist  : |star_dec - ev_dec|
  TS        : -chi2/2 - ln(num_brighter) + ln(4*pi) + ln(num_total)
  best TS   : max over cone stars

Speed trick: only stars within CAND_RADIUS of the event can win the max
whenever at least one star lies within GUARANTEE_RADIUS (worst-case sigma
bound; see pilot notes). Positions with no star inside GUARANTEE_RADIUS fall
back to scanning the full 1-deg cone.
"""
import numpy as np
from scipy.spatial import cKDTree

ARCSEC = np.deg2rad(1.0 / 3600.0)
DEG = np.deg2rad(1.0)
BG_RADIUS = 1.0 * DEG
CAND_RADIUS = 0.15 * DEG
GUARANTEE_RADIUS = 0.12 * DEG
# floors hard-coded in the Feb-2023 get_gaia_association (== current defaults)
FLOOR_RA = 3.73 * ARCSEC
FLOOR_DEC = 4.35 * ARCSEC


def _unit_vectors(ra_rad, dec_rad):
    c = np.cos(dec_rad)
    return np.column_stack((c * np.cos(ra_rad), c * np.sin(ra_rad),
                            np.sin(dec_rad)))


def _chord(angle_rad):
    return 2.0 * np.sin(angle_rad / 2.0)


def vincenty_sep(ra1, dec1, ra2, dec2):
    """Great-circle separation, atan2 (Vincenty) form == astropy separation."""
    sd1, cd1 = np.sin(dec1), np.cos(dec1)
    sd2, cd2 = np.sin(dec2), np.cos(dec2)
    dra = ra2 - ra1
    sdr, cdr = np.sin(dra), np.cos(dra)
    num = np.hypot(cd2 * sdr, cd1 * sd2 - sd1 * cd2 * cdr)
    den = sd1 * sd2 + cd1 * cd2 * cdr
    return np.arctan2(num, den)


class NullTSEngine:
    def __init__(self, cat_ra_deg, cat_dec_deg, cat_gmag):
        self.ra = np.deg2rad(np.asarray(cat_ra_deg, dtype=np.float64))
        self.dec = np.deg2rad(np.asarray(cat_dec_deg, dtype=np.float64))
        self.gmag = np.asarray(cat_gmag, dtype=np.float64)
        self.tree = cKDTree(_unit_vectors(self.ra, self.dec))

    def best_ts(self, ev_ra_deg, ev_dec_deg, ev_ts, batch=2048, workers=-1,
                progress=None, floor_ra=FLOOR_RA, floor_dec=FLOOR_DEC,
                bg_radius=BG_RADIUS, cand_radius=None, guarantee_radius=None):
        """Best (max) association_TS per event; NaN if the bg cone is empty.

        floor_ra/floor_dec in radians (defaults: winter-field 3.73"/4.35";
        pass 3.05"/4.42" [Wan+2025] for the Galactic field). bg_radius sets
        the background cone; the scoring run for the 73 sources fed 0.12-deg
        catalogs to get_gaia_prob, so the Galactic null must use 0.12 too.
        cand/guarantee radii auto-scale if not given.
        """
        if cand_radius is None:
            cand_radius = min(CAND_RADIUS, bg_radius)
        if guarantee_radius is None:
            guarantee_radius = 0.8 * cand_radius
        ev_ra = np.deg2rad(np.asarray(ev_ra_deg, dtype=np.float64))
        ev_dec = np.deg2rad(np.asarray(ev_dec_deg, dtype=np.float64))
        ev_ts = np.asarray(ev_ts, dtype=np.float64)
        n = len(ev_ra)
        sigma_ra = np.sqrt(floor_ra ** 2 + (50 * ARCSEC / np.sqrt(ev_ts)) ** 2)
        sigma_dec = np.sqrt(floor_dec ** 2 + (49 * ARCSEC / np.sqrt(ev_ts)) ** 2)

        out = np.full(n, np.nan)
        ev_vec = _unit_vectors(ev_ra, ev_dec)
        for lo in range(0, n, batch):
            hi = min(lo + batch, n)
            vec = ev_vec[lo:hi]
            cones = self.tree.query_ball_point(
                vec, _chord(bg_radius), workers=workers)
            cands = self.tree.query_ball_point(
                vec, _chord(cand_radius), workers=workers)
            for k in range(hi - lo):
                i = lo + k
                cone_idx = np.asarray(cones[k], dtype=np.intp)
                if cone_idx.size == 0:
                    continue
                cand_idx = np.asarray(cands[k], dtype=np.intp)
                use_full = cand_idx.size == 0
                if not use_full:
                    d_cand = vincenty_sep(ev_ra[i], ev_dec[i],
                                          self.ra[cand_idx], self.dec[cand_idx])
                    if d_cand.min() > guarantee_radius:
                        use_full = True
                if use_full:
                    cand_idx = cone_idx
                sorted_mags = np.sort(self.gmag[cone_idx])
                rank = sorted_mags.searchsorted(
                    self.gmag[cand_idx], side="left") + 1
                ra_dist = vincenty_sep(ev_ra[i], ev_dec[i],
                                       self.ra[cand_idx], ev_dec[i])
                dec_dist = np.abs(self.dec[cand_idx] - ev_dec[i])
                chi2 = (ra_dist / sigma_ra[i]) ** 2 + (dec_dist / sigma_dec[i]) ** 2
                ts = (-0.5 * chi2 - np.log(rank)
                      + np.log(4 * np.pi) + np.log(cone_idx.size))
                out[i] = ts.max()
            if progress is not None and (hi // batch) % 16 == 0:
                progress(hi, n)
        return out
