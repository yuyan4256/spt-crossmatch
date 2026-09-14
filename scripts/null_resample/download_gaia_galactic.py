#!/usr/bin/env python
"""Download Gaia DR3 G<=20 over the SPT-3G Galactic field (Wan+2025 footprint).

Same source_id healpix-range strategy as the winter download, but adaptive:
bulge densities reach ~1M/deg2, so a fixed level would blow VizieR's MAXREC.
Start at level 5 (3.36 deg2); on overflow split the pixel into its 4 children
and recurse (level 8 = 0.05 deg2 handles anything the bulge offers).

Box: coadd footprint RA [252.1,275.7] Dec [-37.5,-26.7] padded 0.3 deg
(covers the 0.12/0.2-deg scoring cones; extend later if 1-deg cones needed).

3 parallel workers over top-level pixels; per-tile npz, resume-safe.
Run:  <spt3g-env python> download_gaia_galactic.py
"""
import glob
import io
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from healpix_nest import ang2pix_nest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "gaia_null")
OUTDIR = os.path.join(DATADIR, "tiles_galactic")
os.makedirs(OUTDIR, exist_ok=True)

URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
MAX_G = 20.0
MAXREC = 900000
RA_LO, RA_HI = 251.8, 276.0
DEC_LO, DEC_HI = -37.8, -26.4
BASE_LEVEL = 5
MAX_LEVEL = 8


def sid_range(pix, level):
    shift = 35 + 2 * (12 - level)
    return pix << shift, ((pix + 1) << shift) - 1


def tile_path(pix, level):
    return os.path.join(OUTDIR, f"hpx{level}_{pix:06d}.npz")


def fetch(pix, level):
    a, b = sid_range(pix, level)
    q = ('SELECT Source,RA_ICRS,DE_ICRS,Gmag FROM "I/355/gaiadr3" '
         f"WHERE Source BETWEEN {a} AND {b} AND Gmag <= {MAX_G}")
    for attempt in (1, 2, 3, 4, 5):
        try:
            r = requests.post(
                URL,
                data=dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv",
                          MAXREC=str(MAXREC), QUERY=q),
                headers={"Connection": "close"}, timeout=(30, 240))
            r.raise_for_status()
            return pd.read_csv(io.StringIO(r.text))
        except Exception as e:
            print(f"hpx{level}_{pix} attempt {attempt} "
                  f"{type(e).__name__}: {e}", flush=True)
            time.sleep(10 * attempt)
    raise RuntimeError(f"hpx{level}_{pix} failed after retries")


def get_tile(pix, level):
    """Download pixel (recursing into children on MAXREC overflow)."""
    if os.path.isfile(tile_path(pix, level)):
        return 0
    children = [(pix << 2) + i for i in range(4)]
    if level < MAX_LEVEL and any(
            glob.glob(os.path.join(OUTDIR, f"hpx{lv}_{c << (2*(lv-level-1)):06d}*"))
            or glob.glob(os.path.join(OUTDIR, f"hpx{level+1}_{c:06d}.npz"))
            for lv in (level + 1,) for c in children):
        # a previous run already split this pixel
        return sum(get_tile(c, level + 1) for c in children)
    df = fetch(pix, level)
    if len(df) >= MAXREC and level < MAX_LEVEL:
        print(f"hpx{level}_{pix}: {len(df)} rows >= MAXREC, splitting",
              flush=True)
        return sum(get_tile(c, level + 1) for c in children)
    tmp = tile_path(pix, level) + ".tmp.npz"
    np.savez_compressed(
        tmp,
        source_id=df["Source"].to_numpy(np.int64),
        ra=df["RA_ICRS"].to_numpy(np.float64),
        dec=df["DE_ICRS"].to_numpy(np.float64),
        g_mag=df["Gmag"].to_numpy(np.float64))
    os.replace(tmp, tile_path(pix, level))
    print(f"hpx{level}_{pix}: {len(df)} rows saved", flush=True)
    time.sleep(0.2)
    return len(df)


def worker(pix):
    try:
        return get_tile(int(pix), BASE_LEVEL)
    except Exception as e:
        print(f"WORKER FAILED hpx{BASE_LEVEL}_{pix}: {e}", flush=True)
        return -1


def main():
    nside = 2 ** BASE_LEVEL
    ras = np.arange(RA_LO - 1.5, RA_HI + 1.5, 0.05)
    decs = np.arange(DEC_LO - 1.5, DEC_HI + 1.5, 0.05)
    G = np.meshgrid(ras, decs)
    pixels = np.unique(ang2pix_nest(nside, G[0].ravel(), G[1].ravel()))
    print(f"{len(pixels)} level-{BASE_LEVEL} top pixels", flush=True)
    t0 = time.time()
    with Pool(3) as pool:
        rows = pool.map(worker, pixels)
    bad = sum(1 for r in rows if r == -1)
    print(f"done in {time.time()-t0:.0f}s, new rows {sum(r for r in rows if r > 0)}, "
          f"failed subtrees {bad}", flush=True)


if __name__ == "__main__":
    main()
