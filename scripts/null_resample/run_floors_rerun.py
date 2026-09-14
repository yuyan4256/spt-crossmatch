#!/usr/bin/env python
"""Re-run the Galactic-field null and re-score the 73 sources with pointing
floors that include Wan+2025's relative-pointing term.

Wan+2025 (2509.08962) Sec 4.4 builds the positional uncertainty from three
terms in quadrature:

    sigma^2 = sigma_sys,abs^2 + sigma_sys,rel^2 + (theta_beam / SNR)^2

with sigma_sys,abs = 3.05" (xdecl.) / 4.42" (decl.), theta_beam = 50" / 49",
and sigma_sys,rel measured per flare event: Source 1 (1.67", 1.94"),
Source 2 (2.24", 2.05").  The production run so far used only the absolute
term as the floor (get_gaia_prob's pointing_floor_*), leaving the relative
term out because it is defined on a flare window.  This script adds it, in
two conventions, and keeps the current floors as the reference column:

    F0_abs        3.05 / 4.42   what outputs/gaia_pchance_galactic.csv used
    F1_abs_rel_mean             abs (+) mean of the two events' rel terms
    F2_abs_rel_max              abs (+) larger of the two events, per axis

Everything else (catalog, footprint sampling, 0.12 deg cone, TS grid, 250k
events per curve) is identical to make_galactic_mappings.py, so the F0
columns must reproduce the existing null -- that is the built-in check.

Outputs
    data/gaia_null/ts_floors_n250000.npz     best TS, (n_events, 39)
    outputs/floors_source_ts.csv             per-candidate TS for the 3 floors
Run:  python run_floors_rerun.py            (~15 min, 162M-star tree)
"""
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "scripts", "racs_pchance"))

from null_multi_sigma import MultiSigmaEngine, ARCSEC, DEG, vincenty_sep
from run_null_galactic import sample_footprint
from make_galactic_mappings import TS_GRID, N_PER_TS

DATADIR = os.path.join(REPO, "data", "gaia_null")
CAT = os.path.join(DATADIR, "gaia_galactic_g20.npz")
def _find_src():
    """The F0 scoring table; the file has been renamed by hand before
    ('gaia_pchance_galactic 3.05\'\' 4.42\'\'.csv'), so match on the stem."""
    import glob
    cands = sorted(glob.glob(os.path.join(REPO, "outputs",
                                          "gaia_pchance_galactic*.csv")))
    if not cands:
        raise FileNotFoundError("outputs/gaia_pchance_galactic*.csv")
    return cands[0]


SRC = _find_src()
NULL_OUT = os.path.join(DATADIR, "ts_floors_n%s.npz" %
                        os.environ.get("NFLOORS", "250000"))
SRC_OUT = os.path.join(REPO, "outputs", "floors_source_ts.csv")

REL = {"src1": (1.67, 1.94), "src2": (2.24, 2.05)}          # Wan+2025 Sec 4.4
ABS_RA, ABS_DEC = 3.05, 4.42
BEAM_RA, BEAM_DEC = 50.0, 49.0

FLOORS = {
    "F0_abs": (ABS_RA, ABS_DEC),
    "F1_abs_rel_mean": (
        float(np.hypot(ABS_RA, np.mean([REL["src1"][0], REL["src2"][0]]))),
        float(np.hypot(ABS_DEC, np.mean([REL["src1"][1], REL["src2"][1]]))),
    ),
    "F2_abs_rel_max": (
        float(np.hypot(ABS_RA, max(REL["src1"][0], REL["src2"][0]))),
        float(np.hypot(ABS_DEC, max(REL["src1"][1], REL["src2"][1]))),
    ),
    # control: the spt3g default / winter-field floor (Tandoi+2024 Eqs. 2-3),
    # i.e. what the shipped p-value table assumes. Not a floor for this field,
    # carried so the comparison "official floor + our null" exists.
    "F3_winter_default": (3.73, 4.35),
}

SEED = 3000                      # == the TS=50 member of the existing family
BG = 0.12 * DEG
SEARCH_ARCSEC = 10.0             # same candidate cut as gaia_official_10arcsec


def sigma_pairs():
    """(labels, pairs) for every (floor, event TS) combination, arcsec."""
    labels, pairs = [], []
    for name, (fra, fdec) in FLOORS.items():
        for ts in TS_GRID:
            labels.append((name, ts))
            pairs.append((np.hypot(fra, BEAM_RA / np.sqrt(ts)),
                          np.hypot(fdec, BEAM_DEC / np.sqrt(ts))))
    return labels, np.asarray(pairs, float)


def score_sources(eng, df):
    """Per-candidate association TS at the SPT centroid, exact full cone."""
    out = {name: np.full(len(df), np.nan) for name in FLOORS}
    out["n_cone"] = np.full(len(df), np.nan)
    out["n_brighter"] = np.full(len(df), np.nan)
    for sid, g in df.groupby("spt_id"):
        ra_e = np.deg2rad(g["spt_ra_deg"].iloc[0])
        dec_e = np.deg2rad(g["spt_dec_deg"].iloc[0])
        ts_ev = float(g["spt_snr_max"].iloc[0]) ** 2
        vec = np.array([[np.cos(dec_e) * np.cos(ra_e),
                         np.cos(dec_e) * np.sin(ra_e), np.sin(dec_e)]])
        cone = np.asarray(eng.tree.query_ball_point(
            vec, 2.0 * np.sin(BG / 2.0), workers=-1)[0], dtype=np.intp)
        if cone.size == 0:
            continue
        # candidates: the same Gaia stars the production run kept (within 10")
        sep = vincenty_sep(ra_e, dec_e, eng.ra[cone], eng.dec[cone])
        cand = cone[sep <= SEARCH_ARCSEC * ARCSEC]
        if cand.size == 0:
            continue
        rank = np.sort(eng.gmag[cone]).searchsorted(eng.gmag[cand], "left") + 1
        base = -np.log(rank) + np.log(4 * np.pi) + np.log(cone.size)
        ra_d = vincenty_sep(ra_e, dec_e, eng.ra[cand], dec_e)
        dec_d = np.abs(eng.dec[cand] - dec_e)
        # match our rows to the catalog candidates by sky position
        rows = g.index.values
        r_ra = np.deg2rad(g["ra"].values)
        r_dec = np.deg2rad(g["dec"].values)
        match = np.full(len(rows), -1, dtype=int)
        for j in range(len(rows)):
            d = vincenty_sep(r_ra[j], r_dec[j], eng.ra[cand], eng.dec[cand])
            k = int(d.argmin())
            if d[k] < 0.5 * ARCSEC:
                match[j] = k
        for j, ridx in enumerate(rows):
            if match[j] >= 0:
                out["n_cone"][ridx] = cone.size
                out["n_brighter"][ridx] = rank[match[j]]
        for name, (fra, fdec) in FLOORS.items():
            s_ra = np.hypot(fra, BEAM_RA / np.sqrt(ts_ev)) * ARCSEC
            s_dec = np.hypot(fdec, BEAM_DEC / np.sqrt(ts_ev)) * ARCSEC
            ts_cand = -0.5 * ((ra_d / s_ra) ** 2 + (dec_d / s_dec) ** 2) + base
            for j, ridx in enumerate(rows):
                if match[j] >= 0:
                    out[name][ridx] = ts_cand[match[j]]
    return out


def main():
    for name, (fra, fdec) in FLOORS.items():
        print(f"floor {name:16s} = {fra:.4f}\" / {fdec:.4f}\"", flush=True)

    d = np.load(CAT)
    print(f"catalog: {len(d['ra'])} stars", flush=True)
    t0 = time.time()
    eng = MultiSigmaEngine(d["ra"], d["dec"], d["g_mag"])
    print(f"tree built in {time.time()-t0:.0f}s", flush=True)

    # ---- 1. sources (cheap, do first so a crash in the null keeps it) ----
    df = pd.read_csv(SRC).reset_index(drop=True)
    t0 = time.time()
    cols = score_sources(eng, df)
    for name in FLOORS:
        df[f"association_TS_{name}"] = cols[name]
    df["n_cone_local"] = cols["n_cone"]
    df["n_brighter_local"] = cols["n_brighter"]
    dts = df["association_TS_F0_abs"] - df["association_TS"]
    print(f"sources scored in {time.time()-t0:.0f}s; F0 vs production TS: "
          f"median {np.nanmedian(dts):+.4f}, p95 |d| "
          f"{np.nanpercentile(np.abs(dts), 95):.4f}, "
          f"max |d| {np.nanmax(np.abs(dts)):.4f}, "
          f"{np.isnan(dts).sum()} unmatched", flush=True)
    df.to_csv(SRC_OUT, index=False)
    print("saved", SRC_OUT, flush=True)

    # ---- 2. null: one footprint sample, all 39 sigma pairs in one pass ----
    labels, pairs = sigma_pairs()
    rng = np.random.default_rng(SEED)
    n_ev = int(os.environ.get("NFLOORS", N_PER_TS))
    ev_ra, ev_dec = sample_footprint(n_ev, rng)
    t0 = time.time()

    def prog(i, n):
        el = time.time() - t0
        print(f"  {i}/{n} ({el:.0f}s, eta {el/max(i,1)*(n-i):.0f}s)", flush=True)

    best = eng.best_ts(ev_ra, ev_dec, pairs,
                       bg_radius=BG, cand_radius=60 * ARCSEC,
                       guarantee_radius=40 * ARCSEC, progress=prog)
    print(f"null done in {time.time()-t0:.0f}s", flush=True)

    # pruning check on a subset, exact vs pruned (verify_pruning.py logic)
    sub = slice(0, min(2000, n_ev))
    ex = eng.best_ts(ev_ra[sub], ev_dec[sub], pairs, bg_radius=BG, exact=True)
    dd = np.abs(ex - best[sub])
    print(f"pruning check (2000 events x {len(pairs)} sigmas): max |diff| "
          f"{np.nanmax(dd):.3e}, any pruned>exact: "
          f"{bool(np.nanmax(best[sub] - ex) > 1e-9)}", flush=True)

    np.savez(NULL_OUT, best_ts=best, ev_ra=ev_ra, ev_dec=ev_dec, seed=SEED,
             sigma_pairs=pairs,
             floor_names=np.array([l[0] for l in labels]),
             event_ts=np.array([l[1] for l in labels], float),
             bg_deg=0.12)
    print("saved", NULL_OUT, flush=True)

    # ---- 3. F0 must reproduce the existing TS=50 curve (same seed) ----
    print("\nF0 vs the existing 250k Galactic null (TS at fixed p):", flush=True)
    print(f"{'event TS':>9} {'p':>6} {'TS_new':>8} {'TS_old':>8} {'diff':>7}",
          flush=True)
    for tsv in TS_GRID:
        ref = os.path.join(DATADIR, f"ts_gal_ts{tsv:g}_n250000.npz")
        if not os.path.isfile(ref):
            continue
        old = np.sort(np.load(ref)["best_ts"])
        col = [i for i, l in enumerate(labels) if l == ("F0_abs", tsv)][0]
        new = np.sort(best[:, col][np.isfinite(best[:, col])])
        for p in (0.05, 0.02, 0.01):
            a = np.quantile(new, 1 - p)
            b = np.quantile(old, 1 - p)
            print(f"{tsv:9g} {p:6.2f} {a:8.3f} {b:8.3f} {a - b:+7.3f}",
                  flush=True)
    if n_ev == 250000 and SEED == 3000:
        old = np.load(os.path.join(DATADIR, "ts_gal_ts50_n250000.npz"))["best_ts"]
        col = [i for i, l in enumerate(labels) if l == ("F0_abs", 50)][0]
        print(f"bitwise (same positions): max |diff| "
              f"{np.nanmax(np.abs(best[:, col] - old)):.3e}", flush=True)


if __name__ == "__main__":
    main()
