#!/usr/bin/env python
"""Score the AGN-marked sources at their RACS-mid positions and turn the
association_TS into a chance probability with the Galactic-plane null.

For every source and every sigma variant (see racs_configs.py):
  - association_TS of the best Gaia DR3 (G<=20) star, exact full 0.12-deg cone
  - p_chance = P(null best_TS >= observed), read off the empirical survival
    function of the matching null sigma pair, cubic-interpolated exactly as
    scripts/null_resample/rescore_73_galactic.py does for the SPT run

Writes outputs/racs_agn_pchance.csv (long) and outputs/racs_agn_pchance_wide.csv
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO_, "src"))
from pchance import survival_interp
from null_multi_sigma import MultiSigmaEngine, ARCSEC, DEG, vincenty_sep
from racs_configs import sigmas, VARIANTS

REPO = os.path.dirname(os.path.dirname(HERE))
DATADIR = os.path.join(REPO, "data", "gaia_null")
NULL = os.path.join(DATADIR, "racs_null_n250000_seed7001.npz")



def score(eng, ra_deg, dec_deg, cfg, bg=0.12 * DEG):
    """Exact full-cone best TS + winning star, for every sigma pair."""
    sig = np.asarray(cfg, float) * ARCSEC
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    vec = np.array([[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                     np.sin(dec)]])
    idx = np.asarray(eng.tree.query_ball_point(
        vec, 2.0 * np.sin(bg / 2.0), workers=-1)[0], dtype=np.intp)
    mags = eng.gmag[idx]
    rank = np.sort(mags).searchsorted(mags, "left") + 1
    ra_dist = vincenty_sep(ra, dec, eng.ra[idx], dec)
    dec_dist = np.abs(eng.dec[idx] - dec)
    base = -np.log(rank) + np.log(4 * np.pi) + np.log(idx.size)
    chi2 = (ra_dist[None, :] / sig[:, 0][:, None]) ** 2 + \
           (dec_dist[None, :] / sig[:, 1][:, None]) ** 2
    ts = -0.5 * chi2 + base[None, :]
    best = ts.argmax(axis=1)
    sep = vincenty_sep(ra, dec, eng.ra[idx], eng.dec[idx])
    return dict(
        ts=ts[np.arange(len(sig)), best], star=idx[best], n_cone=idx.size,
        sep=np.rad2deg(sep[best]) * 3600.0, gmag=mags[best], rank=rank[best])


def main():
    src = pd.read_csv(os.path.join(REPO, "outputs", "racs_agn_positions.csv"))
    nz = np.load(NULL)
    # The null was run over the sigma grid of every variant ever defined, so
    # look each source's sigma up in it rather than requiring the grid to match
    # the current VARIANTS exactly (retiring a variant must not force a re-run).
    cfg = nz["cfg"]
    lookup = {(a, b): j for j, (a, b) in enumerate(map(tuple, np.round(cfg, 3)))}
    index = {}
    for v in VARIANTS:
        sr, sd = sigmas(src, v)
        index[v] = np.array([lookup[(a, b)] for a, b in
                             zip(np.round(sr, 3), np.round(sd, 3))])
    null_ts = nz["best_ts"].astype(float)
    interps = [(survival_interp(c), int(np.isfinite(c).sum()))
               for c in (null_ts[:, j] for j in range(len(cfg)))]
    print(f"null: {null_ts.shape[0]} events, {len(cfg)} sigma pairs", flush=True)

    d = np.load(os.path.join(DATADIR, "gaia_galactic_g20.npz"))
    eng = MultiSigmaEngine(d["ra"], d["dec"], d["g_mag"])
    sid = d["source_id"]
    print("tree built", flush=True)

    rows = []
    for i, s in src.iterrows():
        r = score(eng, s["racs_RAJ2000"], s["racs_DEdegc"], cfg)
        for v in VARIANTS:
            j = index[v][i]
            f, nfin = interps[j]
            rows.append(dict(
                name=s["name"], classification=s["classification"],
                agn_group=s["agn_group"], variant=v,
                sigma_ra_arcsec=cfg[j, 0], sigma_dec_arcsec=cfg[j, 1],
                racs_ra=s["racs_RAJ2000"], racs_dec=s["racs_DEdegc"],
                racs_sep_from_spt_arcsec=s["racs_sep_arcsec"],
                racs_Ftot_mJy=s["racs_Ftot"], racs_SCode=s["racs_SCode"],
                racs_Flag=s["racs_Flag"],
                association_TS=r["ts"][j], gaia_source_id=sid[r["star"][j]],
                gaia_ra=np.rad2deg(eng.ra[r["star"][j]]),
                gaia_dec=np.rad2deg(eng.dec[r["star"][j]]),
                gaia_sep_arcsec=r["sep"][j], gaia_g_mag=r["gmag"][j],
                num_brighter=r["rank"][j], n_cone=r["n_cone"],
                p_chance_racs=float(f(r["ts"][j])),
                null_ts_max=np.nanmax(null_ts[:, j]),
                p_floor=1.0 / nfin))
        nv = len(VARIANTS)
        print(f"  {s['name']}: "
              + "  ".join(f"{v.split('_')[0]} TS={rows[k - nv]['association_TS']:.1f} "
                          f"p={rows[k - nv]['p_chance_racs']:.4f}"
                          for k, v in enumerate(VARIANTS)), flush=True)

    out = pd.DataFrame(rows)
    p1 = os.path.join(REPO, "outputs", "racs_agn_pchance.csv")
    out.to_csv(p1, index=False)

    wide = out.pivot_table(index="name", columns="variant",
                           values=["association_TS", "p_chance_racs"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    key = out[out.variant == "A_fit_plus_floor"].set_index("name")
    wide = key[["classification", "agn_group", "racs_ra", "racs_dec",
                "racs_sep_from_spt_arcsec", "racs_Ftot_mJy", "sigma_ra_arcsec",
                "sigma_dec_arcsec", "gaia_source_id", "gaia_sep_arcsec",
                "gaia_g_mag", "num_brighter"]].join(wide)

    # side-by-side with the SPT-position Galactic-null run
    spt = os.path.join(REPO, "outputs", "gaia_pchance_galactic.csv")
    if os.path.isfile(spt):
        g = pd.read_csv(spt)
        b = g.loc[g.groupby("spt_id")["association_TS"].idxmax()]
        wide = wide.join(b.set_index("spt_id")[
            ["association_TS", "p_galactic", "gaia_sep_arcsec", "source_id"]]
            .rename(columns={"association_TS": "association_TS_spt_pos",
                             "p_galactic": "p_galactic_spt_pos",
                             "gaia_sep_arcsec": "gaia_sep_arcsec_spt_pos",
                             "source_id": "gaia_source_id_spt_pos"}))
    p2 = os.path.join(REPO, "outputs", "racs_agn_pchance_wide.csv")
    wide.sort_values("p_chance_racs_A_fit_plus_floor").to_csv(p2)
    print("saved", p1, "and", p2)

    for v in VARIANTS:
        a = out[out.variant == v]
        for thr in (0.01, 0.02, 0.05):
            print(f"{v}: {int((a.p_chance_racs < thr).sum())}/{len(a)} "
                  f"with p < {thr}")


if __name__ == "__main__":
    main()
