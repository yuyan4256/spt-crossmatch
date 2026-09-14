#!/usr/bin/env python
"""AllWISE null for the SPT-3G Galactic field.

Identical in every respect to scripts/null_resample/run_null_galactic.py --
random positions drawn uniformly over the valid coadd footprint, best
association_TS in a 0.12-deg background cone, sigma = repo pointing floor (+)
50"/49" over sqrt(TS) -- except the counterpart catalog is AllWISE ranked by
w1mpro instead of Gaia DR3 ranked by Gmag.

All 26 sigma pairs (13 event-TS grid points x 2 floors, see wise_configs.py)
are scored in one pass, so the expensive cone query is done once per event.

Usage: run_null_wise.py --n 400000 --seed 9001
Saves data/wise_null/wise_null_n{n}_seed{seed}.npz
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "scripts", "racs_pchance"))
sys.path.insert(0, os.path.join(REPO, "scripts", "null_resample"))

from null_multi_sigma import MultiSigmaEngine, ARCSEC, DEG
from run_null_galactic import sample_footprint
from wise_configs import config_table

DATADIR = os.path.join(REPO, "data", "wise_null")
BG_DEG = 0.12
CAND = 60.0 * ARCSEC
GUAR = 40.0 * ARCSEC


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=400000)
    ap.add_argument("--seed", type=int, default=9001)
    ap.add_argument("--bg-deg", type=float, default=BG_DEG)
    args = ap.parse_args()

    cfg, index = config_table()
    print(f"{len(cfg)} sigma pairs: ra {cfg[:,0].min():.2f}-{cfg[:,0].max():.2f}\", "
          f"dec {cfg[:,1].min():.2f}-{cfg[:,1].max():.2f}\"", flush=True)

    d = np.load(os.path.join(DATADIR, "allwise_field.npz"))
    print(f"catalog: {d['ra'].size} AllWISE sources", flush=True)
    t0 = time.time()
    eng = MultiSigmaEngine(d["ra"], d["dec"], d["w1mpro"])
    print(f"tree built in {time.time()-t0:.0f}s", flush=True)

    rng = np.random.default_rng(args.seed)
    ev_ra, ev_dec = sample_footprint(args.n, rng)
    print(f"{args.n} footprint positions drawn", flush=True)

    t0 = time.time()
    last = [0.0]

    def prog(i, n_):
        if time.time() - last[0] > 60:
            last[0] = time.time()
            el = time.time() - t0
            print(f"  {i}/{n_} ({el:.0f}s, eta {el*(n_-i)/max(i,1):.0f}s)",
                  flush=True)

    bg = args.bg_deg * DEG
    ts = eng.best_ts(ev_ra, ev_dec, cfg, bg_radius=bg,
                     cand_radius=min(CAND, bg), guarantee_radius=min(GUAR, 0.8 * bg),
                     progress=prog)
    dt = time.time() - t0
    print(f"{args.n} events x {len(cfg)} sigmas in {dt:.0f}s "
          f"({1000*dt/args.n:.2f} ms/event)", flush=True)
    nan = np.isnan(ts[:, 0]).sum()
    print(f"empty cones: {nan}", flush=True)
    for key in [("gal", 50), ("gal", 500), ("gal", 42000)]:
        col = ts[:, index[key]]
        col = col[np.isfinite(col)]
        print(f"  {key}: med {np.median(col):.2f} p99 {np.percentile(col, 99):.2f} "
              f"max {col.max():.2f}", flush=True)

    out = os.path.join(DATADIR, f"wise_null_n{args.n}_seed{args.seed}.npz")
    np.savez_compressed(out, best_ts=ts.astype(np.float32), cfg=cfg,
                        ev_ra=ev_ra, ev_dec=ev_dec, seed=args.seed,
                        bg_deg=args.bg_deg,
                        keys=np.array([f"{a}|{b}" for a, b in index], dtype=str),
                        key_cols=np.array([index[k] for k in index]))
    print("saved", out, flush=True)
    print("WISE_NULL_DONE", flush=True)


if __name__ == "__main__":
    main()
