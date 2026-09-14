#!/usr/bin/env python
"""Galactic-field null for the RACS-mid-position re-scoring.

Identical to scripts/null_resample/run_null_galactic.py -- random positions
drawn uniformly over the valid coadd footprint, scored against the 162M-star
Galactic Gaia DR3 G<=20 catalog with a 0.12-deg background cone -- except the
positional sigma comes from RACS-mid instead of SPT, so every (sigma_ra,
sigma_dec) pair needed by the three variants is scored in one pass.

Usage: run_null_racs.py --n 250000 --seed 7001
Saves data/gaia_null/racs_null_n{n}_seed{seed}.npz
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "null_resample"))

from null_multi_sigma import MultiSigmaEngine, ARCSEC, DEG
from racs_configs import unique_configs
from run_null_galactic import sample_footprint

REPO = os.path.dirname(os.path.dirname(HERE))
DATADIR = os.path.join(REPO, "data", "gaia_null")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=250000)
    ap.add_argument("--seed", type=int, default=7001)
    ap.add_argument("--exact", action="store_true",
                    help="scan the full 0.12-deg cone for every event; "
                         "required because the circular candidate-radius "
                         "guarantee breaks for strongly anisotropic sigma "
                         "(see verify_pruning.py)")
    args = ap.parse_args()

    src = pd.read_csv(os.path.join(REPO, "outputs", "racs_agn_positions.csv"))
    cfg, index = unique_configs(src)
    print(f"{len(src)} sources -> {len(cfg)} unique sigma pairs "
          f"(ra {cfg[:,0].min():.2f}-{cfg[:,0].max():.2f}\", "
          f"dec {cfg[:,1].min():.2f}-{cfg[:,1].max():.2f}\")", flush=True)

    d = np.load(os.path.join(DATADIR, "gaia_galactic_g20.npz"))
    print(f"catalog: {len(d['ra'])} stars", flush=True)
    t0 = time.time()
    eng = MultiSigmaEngine(d["ra"], d["dec"], d["g_mag"])
    print(f"tree built in {time.time()-t0:.0f}s", flush=True)

    rng = np.random.default_rng(args.seed)
    ev_ra, ev_dec = sample_footprint(args.n, rng)

    t0 = time.time()
    last = [0.0]

    def prog(i, n_):
        if time.time() - last[0] > 60:
            last[0] = time.time()
            el = time.time() - t0
            print(f"  {i}/{n_} ({el:.0f}s, eta {el*(n_-i)/max(i,1):.0f}s)",
                  flush=True)

    ts = eng.best_ts(ev_ra, ev_dec, cfg, bg_radius=0.12 * DEG,
                     cand_radius=60 * ARCSEC, guarantee_radius=40 * ARCSEC,
                     progress=prog, exact=args.exact)
    dt = time.time() - t0
    print(f"{args.n} events x {len(cfg)} sigmas in {dt:.0f}s "
          f"({1000*dt/args.n:.2f} ms/event)", flush=True)

    tag = "_exact" if args.exact else ""
    out = os.path.join(
        DATADIR, f"racs_null_n{args.n}_seed{args.seed}{tag}.npz")
    np.savez_compressed(out, best_ts=ts.astype(np.float32), cfg=cfg,
                        ev_ra=ev_ra, ev_dec=ev_dec, seed=args.seed)
    print("saved", out, flush=True)
    print("RACS_NULL_DONE", flush=True)


if __name__ == "__main__":
    main()
