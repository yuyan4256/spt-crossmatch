"""The (sigma_ra, sigma_dec) variants scored for the RACS-mid AGN re-run.

sigma is built exactly like the existing Gaia pipeline builds it for SPT
(published survey floor in quadrature with the per-source statistical term):

  floor : Duchesne+2024 (RACS-mid, PASA 41, e003) Table 5, primary catalogue
          vs ICRF3 -> d(alpha)cos(dec) = -0.14 +/- 0.65", d(dec)_corr =
          -0.03 +/- 0.97"
  fit   : catalog columns e_RAJ2000 and e_DEdegc (the latter already folds in
          the declination-correction model, hence variant B/C below)

Variants
  A  fit (+) floor        primary; conservative (dec floor partly double-counts
                          the model term already inside e_DEdegc)
  B  floor only           published ICRF3 scatter, identical for every source

  C  fit only             RETIRED 2026-08-23. Catalog formal uncertainties with
                          no systematic term: unphysical (it ignores the +/-1-2"
                          beam-to-beam astrometric scatter documented in
                          Duchesne+2024 Sec 4.10, giving sigma_RA as small as
                          0.002"), and the resulting sigma anisotropy of ~300
                          breaks the candidate-pruning the null relies on
                          (verify_pruning.py). Kept in sigmas() for reference
                          only; not in VARIANTS, so nothing scores it.
"""
import numpy as np

FLOOR_RA = 0.65   # arcsec
FLOOR_DEC = 0.97  # arcsec
VARIANTS = ("A_fit_plus_floor", "B_floor_only")
RETIRED = ("C_fit_only",)  # see module docstring


def sigmas(df, variant):
    """(sigma_ra, sigma_dec) in arcsec for every row of racs_agn_positions."""
    fit_ra = df["e_ra_fit_arcsec"].to_numpy(float)
    fit_dec = df["e_dec_fit_arcsec"].to_numpy(float)
    if variant == "A_fit_plus_floor":
        return np.hypot(fit_ra, FLOOR_RA), np.hypot(fit_dec, FLOOR_DEC)
    if variant == "B_floor_only":
        return np.full(len(df), FLOOR_RA), np.full(len(df), FLOOR_DEC)
    if variant == "C_fit_only":
        return fit_ra, fit_dec
    raise ValueError(variant)


def unique_configs(df, variants=VARIANTS, decimals=3):
    """Deduped list of (sigma_ra, sigma_dec) arcsec pairs to run the null on,
    plus a {variant: index-into-configs per source} lookup."""
    pairs, index = [], {}
    seen = {}
    for v in variants:
        sr, sd = sigmas(df, v)
        idx = []
        for a, b in zip(np.round(sr, decimals), np.round(sd, decimals)):
            key = (a, b)
            if key not in seen:
                seen[key] = len(pairs)
                pairs.append(key)
            idx.append(seen[key])
        index[v] = np.array(idx)
    return np.array(pairs, dtype=float), index
