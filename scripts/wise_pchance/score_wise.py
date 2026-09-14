#!/usr/bin/env python
"""Score the 73 SPT-3G sources against AllWISE and turn association_TS into a
chance probability with the AllWISE field null.

Mirrors scripts/racs_pchance/score_racs.py and
scripts/null_resample/rescore_73_galactic.py exactly:
  - association_TS of the best AllWISE source, EXACT full 0.12-deg cone
    (no candidate pruning for the real sources)
  - sigma from the source's own TS = snr_max^2 and the repo pointing floors
  - p_chance = P(null best_TS >= observed), read off the empirical survival
    function of the null column with the nearest-in-log event TS, cubic
    interpolation, fill (1, 0) -- i.e. get_gaia_prob(ts_dependent=True)

Writes outputs/wise_pchance.csv (one row per source per floor is folded into
columns) and prints the pass counts at p<0.02 / p<0.05.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "src"))
sys.path.insert(0, os.path.join(REPO, "scripts", "racs_pchance"))

from null_multi_sigma import MultiSigmaEngine, ARCSEC, DEG, vincenty_sep
from pchance import survival_interp
from wise_configs import FLOORS, PRIMARY, TS_GRID, sigma, config_table

DATADIR = os.path.join(REPO, "data", "wise_null")
NULL = os.path.join(DATADIR, "wise_null_n400000_seed9001.npz")
TABLE = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
OUT = os.path.join(REPO, "outputs", "wise_pchance.csv")
BG = 0.12 * DEG
MATCH_RADIUS_ARCSEC = 10.0



def score_one(eng, ra_deg, dec_deg, sig_pairs):
    """Exact full-cone best association_TS for each (sigma_ra, sigma_dec)."""
    sig = np.asarray(sig_pairs, float) * ARCSEC
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    vec = np.array([[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                     np.sin(dec)]])
    idx = np.asarray(eng.tree.query_ball_point(
        vec, 2.0 * np.sin(BG / 2.0), workers=-1)[0], dtype=np.intp)
    if idx.size == 0:
        return None
    mags = eng.gmag[idx]
    rank = np.sort(mags).searchsorted(mags, "left") + 1
    ra_dist = vincenty_sep(ra, dec, eng.ra[idx], dec)
    dec_dist = np.abs(eng.dec[idx] - dec)
    base = -np.log(rank) + np.log(4 * np.pi) + np.log(idx.size)
    chi2 = (ra_dist[None, :] / sig[:, 0][:, None]) ** 2 + \
           (dec_dist[None, :] / sig[:, 1][:, None]) ** 2
    ts = -0.5 * chi2 + base[None, :]
    best = ts.argmax(axis=1)
    sep = np.rad2deg(vincenty_sep(ra, dec, eng.ra[idx], eng.dec[idx])) * 3600.0
    return dict(ts=ts[np.arange(len(sig)), best], star=idx[best],
                n_cone=idx.size, sep=sep[best], w1=mags[best],
                rank=rank[best], n10=int((sep <= MATCH_RADIUS_ARCSEC).sum()))


def main():
    src = pd.read_csv(TABLE)
    ev_ts = src["snr_max"].to_numpy(float) ** 2
    print(f"{len(src)} sources, event TS {ev_ts.min():.0f}..{ev_ts.max():.0f}",
          flush=True)

    nz = np.load(NULL)
    cfg, index = np.asarray(nz["cfg"], float), config_table()[1]
    null_ts = nz["best_ts"].astype(float)
    interps = {k: survival_interp(null_ts[:, j]) for k, j in index.items()}
    sorted_null = {k: np.sort(null_ts[np.isfinite(null_ts[:, j]), j])
                   for k, j in index.items()}
    n_null = null_ts.shape[0]
    print(f"null: {null_ts.shape[0]} events, {null_ts.shape[1]} sigma pairs",
          flush=True)

    d = np.load(os.path.join(DATADIR, "allwise_field.npz"))
    eng = MultiSigmaEngine(d["ra"], d["dec"], d["w1mpro"])
    print(f"catalog: {d['ra'].size} AllWISE sources, tree built", flush=True)

    tsgrid = np.array(TS_GRID, float)
    nearest = tsgrid[np.argmin(np.abs(np.log(ev_ts)[:, None]
                                      - np.log(tsgrid)[None, :]), axis=1)]

    names = list(FLOORS)
    rows = []
    for i, r in src.iterrows():
        pairs = [sigma(ev_ts[i], FLOORS[n]) for n in names]
        res = score_one(eng, r["ra_deg"], r["dec_deg"], pairs)
        row = dict(spt_id=r["id"], ra_deg=r["ra_deg"], dec_deg=r["dec_deg"],
                   snr_max=r["snr_max"], event_ts=ev_ts[i],
                   mapping_ts=nearest[i])
        if res is None:
            row.update(n_cone=0)
            rows.append(row)
            continue
        j = res["star"][names.index(PRIMARY)]
        row.update(n_cone=res["n_cone"], n_wise_10arcsec=res["n10"],
                   wise_designation=str(d["designation"][j]),
                   wise_ra=float(d["ra"][j]), wise_dec=float(d["dec"][j]),
                   w1mpro=float(d["w1mpro"][j]), w2mpro=float(d["w2mpro"][j]),
                   cc_flags=str(d["cc_flags"][j]), ext_flg=float(d["ext_flg"][j]),
                   ph_qual=str(d["ph_qual"][j]),
                   sep_arcsec=float(res["sep"][names.index(PRIMARY)]),
                   w1_rank=int(res["rank"][names.index(PRIMARY)]))
        for k, n in enumerate(names):
            sr, sd = pairs[k]
            row[f"sigma_ra_{n}"] = sr
            row[f"sigma_dec_{n}"] = sd
            row[f"association_TS_{n}"] = float(res["ts"][k])
            row[f"p_wise_{n}"] = float(interps[(n, nearest[i])](res["ts"][k]))
            sn = sorted_null[(n, nearest[i])]
            row[f"n_null_ge_{n}"] = int(sn.size - np.searchsorted(
                sn, res["ts"][k], side="left"))
        rows.append(row)
        if (i + 1) % 10 == 0:
            print(f"  {i+1}/{len(src)}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print("saved", OUT, flush=True)

    for n in names:
        col = df[f"p_wise_{n}"]
        print(f"floor {n} {FLOORS[n]}: p<0.02 {(col < 0.02).sum()}, "
              f"p<0.05 {(col < 0.05).sum()} of {col.notna().sum()}")
    print(f"null sample size {n_null} -> smallest resolvable p = "
          f"{1.0/n_null:.2e}; "
          f"{(df[f'n_null_ge_{PRIMARY}'] == 0).sum()} sources are off the top "
          f"of the null (reported p = 0)")
    print("\nlowest p_chance (primary floor):")
    for _, r in df.nsmallest(10, f"p_wise_{PRIMARY}").iterrows():
        print(f"  {r.spt_id}  TS={r[f'association_TS_{PRIMARY}']:.2f} "
              f"sep={r.sep_arcsec:.2f}\" W1={r.w1mpro:.2f} "
              f"p={r[f'p_wise_{PRIMARY}']:.4f} "
              f"(n_null>=TS: {int(r[f'n_null_ge_{PRIMARY}'])})")


if __name__ == "__main__":
    main()
