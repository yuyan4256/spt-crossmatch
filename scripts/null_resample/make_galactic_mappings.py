#!/usr/bin/env python
"""Generate the TS-dependent association->p-value mapping family for the
Galactic field, in the official g3 format (readable by
get_gaia_prob(..., ts_dependent=True, pval_filename=...)).

One process, tree built once, TS values run serially (250k events each).
Also runs a bg-cone invariance check (0.2 vs 0.12 deg at TS=250, 100k) to
quantify the quantile-invariance assumption the 73-source scoring relied on.

Output:
  data/gaia_null/ts_gal_ts{TS}_n250000.npz      per-TS raw samples
  data/gaia_null/Gaia_association_to_pval_mappings_galactic.g3
  outputs/mappings/galactic_mapping_ts{TS}.txt  two-column per-TS tables
Run: <spt3g-env python> make_galactic_mappings.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "src"))
sys.path.insert(0, os.path.expanduser("~/spt3g_software/build"))

import numpy as np
from null_ts import NullTSEngine, DEG
from pchance import survival_table
from run_null_galactic import (
    sample_footprint, FLOOR_RA_GAL, FLOOR_DEC_GAL, CAND_GAL, GUAR_GAL)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "gaia_null")

TS_GRID = [50, 70, 100, 150, 220, 320, 500, 800, 1300, 2500, 6000, 15000, 42000]
N_PER_TS = 250000
N_INVAR = 100000


def run_one(eng, n, seed, event_ts, bg_deg):
    rng = np.random.default_rng(seed)
    ev_ra, ev_dec = sample_footprint(n, rng)
    t0 = time.time()
    bg = bg_deg * DEG
    ts = eng.best_ts(
        ev_ra, ev_dec, np.full(n, float(event_ts)),
        floor_ra=FLOOR_RA_GAL, floor_dec=FLOOR_DEC_GAL,
        bg_radius=bg, cand_radius=min(CAND_GAL, bg),
        guarantee_radius=min(GUAR_GAL, 0.8 * bg), batch=1024)
    print(f"  TS={event_ts} bg={bg_deg}: {n} events in {time.time()-t0:.0f}s, "
          f"med {np.nanmedian(ts):.2f}", flush=True)
    return ev_ra, ev_dec, ts



def main():
    d = np.load(os.path.join(DATADIR, "gaia_galactic_g20.npz"))
    print(f"catalog: {len(d['ra'])} stars", flush=True)
    t0 = time.time()
    eng = NullTSEngine(d["ra"], d["dec"], d["g_mag"])
    print(f"tree built in {time.time()-t0:.0f}s", flush=True)

    results = {}
    for i, tsv in enumerate(TS_GRID):
        out = os.path.join(DATADIR, f"ts_gal_ts{tsv:g}_n{N_PER_TS}.npz")
        if os.path.isfile(out):
            results[tsv] = np.load(out)["best_ts"]
            print(f"  TS={tsv}: cached", flush=True)
            continue
        ev_ra, ev_dec, ts = run_one(eng, N_PER_TS, 3000 + i, tsv, 0.12)
        np.savez(out, ev_ra=ev_ra, ev_dec=ev_dec, best_ts=ts,
                 seed=3000 + i, event_ts=tsv, bg_deg=0.12)
        results[tsv] = ts

    # bg-cone invariance check at TS=250
    inv = os.path.join(DATADIR, "ts_gal_invariance.npz")
    if not os.path.isfile(inv):
        _, _, a = run_one(eng, N_INVAR, 4001, 250, 0.12)
        _, _, b = run_one(eng, N_INVAR, 4001, 250, 0.20)
        np.savez(inv, ts_bg012=a, ts_bg020=b)

    # g3 file in the official ts_dependent format
    from spt3g import core
    fr = core.G3Frame()
    m_assoc = core.G3MapVectorDouble()
    m_pval = core.G3MapVectorDouble()
    for tsv in TS_GRID:
        grid, p = survival_table(results[tsv])
        key = f"{tsv:g}"
        m_assoc[key] = core.G3VectorDouble(grid)
        m_pval[key] = core.G3VectorDouble(p)
        mapdir = os.path.join(REPO, "outputs", "mappings")
        os.makedirs(mapdir, exist_ok=True)
        np.savetxt(
            os.path.join(mapdir, f"galactic_mapping_ts{tsv:g}.txt"),
            np.column_stack([grid, p]), fmt="%.6e",
            header=(f"Galactic-field mapping, event_ts={tsv:g}, "
                    f"n={N_PER_TS}, bg=0.12deg, floors 3.05/4.42 arcsec\n"
                    "association_ts prob"))
    fr["association TS"] = m_assoc
    fr["p-value"] = m_pval
    g3out = os.path.join(
        DATADIR, "Gaia_association_to_pval_mappings_galactic.g3")
    w = core.G3Writer(g3out)
    w(fr)
    w(core.G3Frame(core.G3FrameType.EndProcessing))
    del w
    print("wrote", g3out, flush=True)
    print("ALL_MAPPINGS_DONE", flush=True)


if __name__ == "__main__":
    main()
