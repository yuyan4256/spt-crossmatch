#!/usr/bin/env python
"""Galactic-field null: random positions from the coadd footprint + best
association_TS against the Galactic Gaia catalog.

Differences vs the winter null (all mirroring how the 73 real sources were
scored in gaia_official_10arcsec.py):
  - positions: uniform over VALID pixels of the ZEA coadd (equal-area proj,
    so uniform-in-pixel == uniform-on-sky), not a RA/Dec rectangle
  - pointing floors: 3.05"/4.42" (Wan+2025 Sec 4.4, this field)
  - bg cone: 0.12 deg (the scoring run's BG_CONE_DEG; association_TS depends
    on background only through the brightness quantile, but null and scoring
    must use the SAME cone to calibrate away any residual)
  - event TS: fixed per run (--ts); a family of runs over the sources' TS
    range (snr_max^2 spans ~50..42000) gives a TS-dependent mapping like the
    official g3 file

Usage: run_null_galactic.py --n 400000 --seed 3001 --ts 250
Saves data/gaia_null/ts_gal_ts{ts}_n{n}.npz
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from null_ts import NullTSEngine, ARCSEC, DEG

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "gaia_null")
COADD = os.path.join(REPO, "data", "coadd", "galaxy_3yr_pc_gc_v2_150.fits")

FLOOR_RA_GAL = 3.05 * ARCSEC
FLOOR_DEC_GAL = 4.42 * ARCSEC
BG_GAL = 0.12 * DEG
CAND_GAL = (60.0 / 3600.0) * DEG
GUAR_GAL = (40.0 / 3600.0) * DEG


def sample_footprint(n, rng):
    """Uniform random sky positions over the valid coadd footprint."""
    from astropy.io import fits
    from astropy.wcs import WCS

    f = fits.open(COADD)
    d0, d1 = f[0].data, f[1].data
    valid = (np.isfinite(d0) & (d0 != 0) & np.isfinite(d1) & (d1 != 0))
    ys, xs = np.nonzero(valid)
    pick = rng.integers(0, len(xs), n)
    x = xs[pick] + rng.uniform(-0.5, 0.5, n)
    y = ys[pick] + rng.uniform(-0.5, 0.5, n)
    w = WCS(f[1].header)
    ra, dec = w.pixel_to_world_values(x, y)
    f.close()
    return np.asarray(ra, dtype=np.float64), np.asarray(dec, dtype=np.float64)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=400000)
    ap.add_argument("--seed", type=int, default=3001)
    ap.add_argument("--ts", type=float, required=True)
    ap.add_argument("--bg-deg", type=float, default=0.12)
    ap.add_argument("--tag", default="gal")
    args = ap.parse_args()

    d = np.load(os.path.join(DATADIR, "gaia_galactic_g20.npz"))
    print(f"catalog: {len(d['ra'])} stars", flush=True)
    eng = NullTSEngine(d["ra"], d["dec"], d["g_mag"])
    print("tree built", flush=True)

    rng = np.random.default_rng(args.seed)
    ev_ra, ev_dec = sample_footprint(args.n, rng)

    t0 = time.time()
    bg = args.bg_deg * DEG
    ts = eng.best_ts(
        ev_ra, ev_dec, np.full(args.n, args.ts),
        floor_ra=FLOOR_RA_GAL, floor_dec=FLOOR_DEC_GAL,
        bg_radius=bg, cand_radius=min(CAND_GAL, bg),
        guarantee_radius=min(GUAR_GAL, 0.8 * bg),
        progress=lambda i, n_: print(f"  {i}/{n_} ({time.time()-t0:.0f}s)",
                                     flush=True))
    dt = time.time() - t0
    print(f"{args.n} events in {dt:.0f}s ({1000*dt/args.n:.2f} ms/event)",
          flush=True)
    print(f"TS: min {np.nanmin(ts):.2f} med {np.nanmedian(ts):.2f} "
          f"max {np.nanmax(ts):.2f} | NaN {np.isnan(ts).sum()}", flush=True)

    lab = f"{args.tag}_ts{args.ts:g}_n{args.n}"
    if args.bg_deg != 0.12:
        lab += f"_bg{args.bg_deg:g}"
    out = os.path.join(DATADIR, f"ts_{lab}.npz")
    np.savez(out, ev_ra=ev_ra, ev_dec=ev_dec, best_ts=ts, seed=args.seed,
             event_ts=args.ts, bg_deg=args.bg_deg)
    print("saved", out, flush=True)


if __name__ == "__main__":
    main()
