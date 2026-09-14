#!/usr/bin/env python
"""Completeness check of the merged Galactic catalog vs live VizieR counts."""
import io
import os
import sys

import numpy as np
import pandas as pd
import requests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from null_ts import _unit_vectors, _chord  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "gaia_null")
URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

CIRCLES = [  # (ra, dec, radius_deg) — inside the padded galactic box
    (252.5, -37.0, 0.2),
    (256.0, -34.0, 0.2),
    (264.25, -32.25, 0.2),   # field center, b~0
    (268.0, -30.0, 0.15),    # bulge-side, ~1M/deg2
    (270.0, -27.0, 0.15),
    (275.5, -26.8, 0.2),
]


def vizier_count(ra, dec, r):
    q = ('SELECT COUNT(*) AS n FROM "I/355/gaiadr3" WHERE 1=CONTAINS('
         f"POINT('ICRS',RA_ICRS,DE_ICRS), CIRCLE('ICRS',{ra},{dec},{r})) "
         "AND Gmag <= 20")
    resp = requests.post(URL, data=dict(REQUEST="doQuery", LANG="ADQL",
                                        FORMAT="csv", QUERY=q), timeout=300)
    resp.raise_for_status()
    return int(pd.read_csv(io.StringIO(resp.text))["n"].iloc[0])


def main():
    from scipy.spatial import cKDTree

    d = np.load(os.path.join(DATADIR, "gaia_galactic_g20.npz"))
    tree = cKDTree(_unit_vectors(np.deg2rad(d["ra"]), np.deg2rad(d["dec"])))
    ok = True
    for ra, dec, r in CIRCLES:
        n_local = tree.query_ball_point(
            _unit_vectors(np.deg2rad([ra]), np.deg2rad([dec]))[0],
            _chord(np.deg2rad(r)), return_length=True)
        n_remote = vizier_count(ra, dec, r)
        status = "OK" if n_local == n_remote else "MISMATCH"
        ok &= n_local == n_remote
        print(f"({ra:7.2f},{dec:6.2f}) r={r}: local {n_local:>7} "
              f"vizier {n_remote:>7}  {status}", flush=True)
    print("PASS" if ok else "FAIL")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
