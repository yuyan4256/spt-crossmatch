#!/usr/bin/env python
"""Merge the per-declination-strip AllWISE tiles into one catalog npz.

Output: data/wise_null/allwise_field.npz (ra, dec, w1mpro + quality columns)
"""
import glob
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "wise_null")
KEYS = ("designation", "ra", "dec", "w1mpro", "w1sigmpro", "w1snr", "w2mpro",
        "cc_flags", "ext_flg", "ph_qual")


def main():
    tiles = sorted(glob.glob(os.path.join(DATADIR, "tiles", "dec*.npz")))
    print(f"{len(tiles)} tiles")
    cols = {k: [] for k in KEYS}
    for t in tiles:
        d = np.load(t, allow_pickle=False)
        for k in KEYS:
            cols[k].append(d[k])
    out = {k: np.concatenate(v) for k, v in cols.items()}
    n = out["ra"].size
    good = np.isfinite(out["ra"]) & np.isfinite(out["dec"]) & np.isfinite(out["w1mpro"])
    if not good.all():
        out = {k: v[good] for k, v in out.items()}
        print(f"dropped {n - good.sum()} rows with non-finite ra/dec/w1mpro")
    print(f"{out['ra'].size} sources | W1 {np.min(out['w1mpro']):.2f} .. "
          f"{np.max(out['w1mpro']):.2f} (median {np.median(out['w1mpro']):.2f})")
    path = os.path.join(DATADIR, "allwise_field.npz")
    np.savez(path, **out)
    print("saved", path)


if __name__ == "__main__":
    main()
