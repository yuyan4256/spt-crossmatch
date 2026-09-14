#!/usr/bin/env python
"""Re-score the 73 sources' Gaia matches with the Galactic-field mapping
family, mirroring get_gaia_prob(ts_dependent=True) exactly:
nearest-in-log event-TS curve + cubic interpolation, fill (1, 0).

Compares against the shipped winter-field p_chance already in
outputs/gaia_official_10arcsec.csv. Writes outputs/gaia_pchance_galactic.csv.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.dirname(HERE)), "src"))
from pchance import survival_table
from make_galactic_mappings import TS_GRID, N_PER_TS

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATADIR = os.path.join(REPO, "data", "gaia_null")


def main():
    interps = {}
    for tsv in TS_GRID:
        d = np.load(os.path.join(DATADIR, f"ts_gal_ts{tsv:g}_n{N_PER_TS}.npz"))
        grid, p = survival_table(d["best_ts"])
        interps[tsv] = interp1d(grid, p, kind="cubic", fill_value=(1, 0),
                                bounds_error=False)
    tsvals = np.array(sorted(interps))

    df = pd.read_csv(os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv"))
    ev_ts = df["spt_snr_max"].values ** 2
    nearest = tsvals[np.argmin(
        np.abs(np.log(ev_ts)[:, None] - np.log(tsvals)[None, :]), axis=1)]
    df["event_ts"] = ev_ts
    df["mapping_ts"] = nearest
    df["p_galactic"] = [float(interps[t](a)) for t, a in
                        zip(nearest, df["association_TS"].values)]

    out = os.path.join(REPO, "outputs", "gaia_pchance_galactic.csv")
    df.to_csv(out, index=False)

    best = df.loc[df.groupby("spt_id")["association_TS"].idxmax()].copy()
    print(f"{len(best)} sources with a Gaia match within 10\"")
    for thr in (0.02, 0.05):
        w = best["p_chance"] < thr
        g = best["p_galactic"] < thr
        print(f"p<{thr}: winter-table {w.sum()} pass -> galactic {g.sum()} pass "
              f"({(w & ~g).sum()} lost, {(~w & g).sum()} gained)")
        for _, r in best[w != g].iterrows():
            print(f"  {r.spt_id}: TS={r.association_TS:.2f} evTS={r.event_ts:.0f} "
                  f"p_winter={r.p_chance:.4f} -> p_gal={r.p_galactic:.4f}")
    rat = best["p_galactic"] / best["p_chance"].clip(lower=1e-9)
    print(f"p_gal/p_winter (best matches): median {rat.median():.2f}, "
          f"range [{rat.min():.2f}, {rat.max():.2f}]")
    print("saved", out)


if __name__ == "__main__":
    main()
