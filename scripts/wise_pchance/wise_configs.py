"""Sigma grid for the AllWISE association_TS / p_chance run.

The positional sigma is a property of the SPT event, not of the counterpart
catalog, so it is built exactly as the Gaia pipeline builds it -- published
pointing floor in quadrature with the per-source statistical term:

  sigma_ra  = sqrt(floor_ra^2  + (50" / sqrt(TS))^2)
  sigma_dec = sqrt(floor_dec^2 + (49" / sqrt(TS))^2)       TS = snr_max^2

Floors are the ones already in the repo, unchanged:
  GAL     3.05" / 4.42"  Wan+2025 Sec 4.4, THIS field -- primary, the pair the
                         production Galactic Gaia null and
                         outputs/gaia_pchance_galactic 3.05'' 4.42''.csv use
  WINTER  3.73" / 4.35"  the floors hard-coded in the Feb-2023
                         get_gaia_association (null_ts.FLOOR_RA/FLOOR_DEC),
                         carried along as a cross-check

TS_GRID is the same 13-point event-TS grid as the Galactic mapping family
(make_galactic_mappings.TS_GRID), so a source picks its mapping by the same
nearest-in-log-TS rule and WISE / Gaia p-values are directly comparable.
"""
import numpy as np

BEAM_RA = 50.0   # arcsec, get_gaia_association
BEAM_DEC = 49.0  # arcsec
FLOORS = {"gal": (3.05, 4.42), "winter": (3.73, 4.35)}
PRIMARY = "gal"
TS_GRID = [50, 70, 100, 150, 220, 320, 500, 800, 1300, 2500, 6000, 15000, 42000]


def sigma(ts, floor):
    """(sigma_ra, sigma_dec) in arcsec for event TS and a (ra, dec) floor."""
    ts = np.asarray(ts, float)
    return (np.sqrt(floor[0] ** 2 + (BEAM_RA / np.sqrt(ts)) ** 2),
            np.sqrt(floor[1] ** 2 + (BEAM_DEC / np.sqrt(ts)) ** 2))


def config_table(floors=FLOORS, ts_grid=TS_GRID):
    """(cfg, index): cfg is the (n, 2) arcsec sigma grid the null is run on;
    index[(floor_name, ts)] is the column of cfg holding that combination."""
    cfg, index = [], {}
    for name, fl in floors.items():
        for ts in ts_grid:
            sr, sd = sigma(ts, fl)
            index[(name, ts)] = len(cfg)
            cfg.append((float(sr), float(sd)))
    return np.array(cfg, float), index
