#!/usr/bin/env python
"""WISE vs Gaia chance probability for the 73 sources, same field null.

Both columns come from the identical machinery (association_TS in a 0.12-deg
cone, sigma = repo pointing floor (+) 50"/49" / sqrt(TS), p from a
footprint-sampled null); the only difference is the counterpart catalog --
AllWISE ranked by W1 (this work, 400k null positions) vs Gaia DR3 G<=20
(outputs/gaia_pchance_galactic 3.05'' 4.42''.csv, 250k null positions per TS).

p_chance is uniform on (0,1) for random positions by construction, so a KS
test of the 73 p-values against U(0,1) measures whether the field's WISE
sources are counterparts at all, above the chance-coincidence background.

Writes outputs/wise_vs_gaia_pchance.csv
"""
import os

import pandas as pd
from scipy.stats import kstest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WISE = os.path.join(REPO, "outputs", "wise_pchance.csv")
GAIA = os.path.join(REPO, "outputs", "gaia_pchance_galactic 3.05'' 4.42''.csv")
OUT = os.path.join(REPO, "outputs", "wise_vs_gaia_pchance.csv")


def main():
    w = pd.read_csv(WISE)
    g = pd.read_csv(GAIA)
    best = g.loc[g.groupby("spt_id")["association_TS"].idxmax()][
        ["spt_id", "association_TS", "gaia_sep_arcsec", "p_galactic"]]
    best = best.rename(columns={"association_TS": "association_TS_gaia",
                                "p_galactic": "p_gaia_gal"})
    m = w.merge(best, on="spt_id", how="left")
    m.to_csv(OUT, index=False)
    print("saved", OUT)

    p = m["p_wise_gal"].dropna().to_numpy()
    ks = kstest(p, "uniform")
    print(f"\n{len(p)} sources with a WISE cone")
    for thr in (0.02, 0.05, 0.1):
        print(f"  p_wise < {thr}: {(p < thr).sum()} observed, "
              f"{thr*len(p):.1f} expected by chance")
    print(f"  KS vs U(0,1): D={ks.statistic:.3f}, p={ks.pvalue:.3f}")

    both = m.dropna(subset=["p_gaia_gal"])
    print(f"\n{len(both)} sources scored in both catalogs")
    for thr in (0.02, 0.05):
        a = both["p_wise_gal"] < thr
        b = both["p_gaia_gal"] < thr
        print(f"  p<{thr}: WISE {a.sum()}, Gaia {b.sum()}, both {(a & b).sum()}, "
              f"either {(a | b).sum()}")
    r = both[["p_wise_gal", "p_gaia_gal"]].corr(method="spearman").iloc[0, 1]
    print(f"  Spearman rho(p_wise, p_gaia) = {r:.3f}")


if __name__ == "__main__":
    main()
