"""SPT-3G positional uncertainty.

Wan et al. 2025 (arXiv:2509.08962) Sec 4.4, Eqs 4-6, calibrated ON the
SPT-3G Galactic plane field:
    sigma_RA^2  = 3.05^2 + (50/sqrt(TS))^2 + sigma_rel_RA^2
    sigma_Dec^2 = 4.42^2 + (49/sqrt(TS))^2 + sigma_rel_Dec^2
with TS = SNR^2, so sqrt(TS) = |SNR|.

sigma_sys,abs (3.05/4.42) is the AT20G-vs-average-map offset scatter for this
field. sigma_sys,rel is not included yet: Wan's definition (per-obs maps during
a flare vs the average map) does not transfer to yearly-coadd centroids; the
yearly-map registration analog is pending. The hook is `SIGMA_REL_RA/DEC`,
keyed by obs_max year — empty means 0.
"""
import numpy as np

SYS_ABS_RA, SYS_ABS_DEC = 3.05, 4.42     # arcsec, AT20G vs average map
BEAM_RA, BEAM_DEC = 50.0, 49.0           # arcsec, theta_beam/sqrt(TS)

SIGMA_REL_RA = {}    # e.g. {2023: 0.8, 2024: 0.7}
SIGMA_REL_DEC = {}


def sigma_pos_wan2025(snr, year=None):
    """Radial positional uncertainty [arcsec]: sqrt(sigma_RA^2 + sigma_Dec^2)."""
    s = max(abs(float(snr)), 1e-3)
    rel_ra = SIGMA_REL_RA.get(year, 0.0)
    rel_dec = SIGMA_REL_DEC.get(year, 0.0)
    sig_ra_sq = SYS_ABS_RA**2 + (BEAM_RA / s)**2 + rel_ra**2
    sig_dec_sq = SYS_ABS_DEC**2 + (BEAM_DEC / s)**2 + rel_dec**2
    return float(np.sqrt(sig_ra_sq + sig_dec_sq))
