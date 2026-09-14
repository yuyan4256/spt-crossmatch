#!/usr/bin/env python
"""Merge adaptive galactic tiles into one catalog npz, with sanity checks."""
import glob
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "gaia_null")


def main():
    tiles = sorted(glob.glob(os.path.join(DATADIR, "tiles_galactic", "hpx*.npz")))
    print(f"{len(tiles)} tiles")
    sid, ra, dec, g = [], [], [], []
    for t in tiles:
        d = np.load(t)
        sid.append(d["source_id"])
        ra.append(d["ra"])
        dec.append(d["dec"])
        g.append(d["g_mag"])
    sid = np.concatenate(sid)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    g = np.concatenate(g)
    # adaptive tiling must partition source_id space: duplicates = a pixel
    # downloaded both whole and split
    assert len(np.unique(sid)) == len(sid), "duplicate source_ids after merge"
    assert g.max() <= 20.0
    print(f"{len(sid)} rows | RA [{ra.min():.2f},{ra.max():.2f}] "
          f"Dec [{dec.min():.2f},{dec.max():.2f}] G [{g.min():.2f},{g.max():.2f}]")
    out = os.path.join(DATADIR, "gaia_galactic_g20.npz")
    np.savez(out, source_id=sid, ra=ra, dec=dec, g_mag=g)
    print("saved", out)


if __name__ == "__main__":
    main()
