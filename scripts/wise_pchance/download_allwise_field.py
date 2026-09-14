#!/usr/bin/env python
"""Download AllWISE over the SPT-3G Galactic field, for the WISE version of the
association_TS / p_chance null.

Same footprint box as the Gaia download (download_gaia_galactic.py):
coadd RA [252.1,275.7] Dec [-37.5,-26.7] padded 0.3 deg, so the 0.12-deg
background cones of any event inside the footprint are fully covered.

Catalog: AllWISE Source Catalog (allwise_p3as_psd) via the IRSA TAP service --
the same catalog the v4 crossmatch table already quotes W1 from, so the WISE
p_chance and the WISE fluxes in outputs/ refer to the same sources.
Ranking magnitude: w1mpro (Vega), the analogue of Gaia Gmag in the official
get_gaia_association math.

Tiled in 0.5-deg declination strips (~250k rows each), per-tile npz, resume-safe.
Run: python download_allwise_field.py
"""
import io
import os
import time

import numpy as np
import pandas as pd
import requests

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUTDIR = os.path.join(REPO, "data", "wise_null", "tiles")
os.makedirs(OUTDIR, exist_ok=True)

URL = "https://irsa.ipac.caltech.edu/TAP/sync"
RA_LO, RA_HI = 251.8, 276.0
DEC_LO, DEC_HI = -37.8, -26.4
STRIP = 0.5
COLS = ("designation,ra,dec,w1mpro,w1sigmpro,w1snr,w2mpro,cc_flags,ext_flg,"
        "ph_qual")


def fetch(d0, d1):
    q = (f"SELECT {COLS} FROM allwise_p3as_psd "
         f"WHERE ra BETWEEN {RA_LO} AND {RA_HI} "
         f"AND dec >= {d0} AND dec < {d1} AND w1mpro IS NOT NULL")
    for attempt in (1, 2, 3, 4, 5):
        try:
            r = requests.post(URL, data=dict(QUERY=q, FORMAT="csv"),
                              timeout=(30, 900))
            r.raise_for_status()
            return pd.read_csv(io.StringIO(r.text))
        except Exception as e:
            print(f"  dec[{d0},{d1}) attempt {attempt} {type(e).__name__}: {e}",
                  flush=True)
            time.sleep(10 * attempt)
    raise RuntimeError(f"dec[{d0},{d1}) failed")


def main():
    edges = np.arange(DEC_LO, DEC_HI, STRIP)
    if edges[-1] < DEC_HI - 1e-9:
        edges = np.append(edges, DEC_HI)   # keep the top partial strip
    t0 = time.time()
    total = 0
    for d0, d1 in zip(edges[:-1], edges[1:]):
        path = os.path.join(OUTDIR, f"dec{d0:+07.2f}.npz")
        if os.path.isfile(path):
            total += int(np.load(path)["ra"].size)
            print(f"dec[{d0:.2f},{d1:.2f}) cached", flush=True)
            continue
        df = fetch(d0, d1)
        tmp = path + ".tmp.npz"
        np.savez_compressed(
            tmp,
            designation=df["designation"].to_numpy(str),
            ra=df["ra"].to_numpy(np.float64),
            dec=df["dec"].to_numpy(np.float64),
            w1mpro=df["w1mpro"].to_numpy(np.float64),
            w1sigmpro=pd.to_numeric(df["w1sigmpro"], errors="coerce")
                .to_numpy(np.float64),
            w1snr=pd.to_numeric(df["w1snr"], errors="coerce")
                .to_numpy(np.float64),
            w2mpro=pd.to_numeric(df["w2mpro"], errors="coerce")
                .to_numpy(np.float64),
            cc_flags=df["cc_flags"].to_numpy(str),
            ext_flg=pd.to_numeric(df["ext_flg"], errors="coerce")
                .to_numpy(np.float64),
            ph_qual=df["ph_qual"].to_numpy(str))
        os.replace(tmp, path)
        total += len(df)
        print(f"dec[{d0:.2f},{d1:.2f}) {len(df)} rows "
              f"({time.time()-t0:.0f}s, {total} total)", flush=True)
    print(f"DOWNLOAD_DONE {total} rows in {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
