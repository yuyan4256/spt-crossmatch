"""Identify which of the 73 sources correspond to points read off the
Distance vs nu*L_nu (150 GHz) figure.

Method: nu*L_nu = 4 pi d^2 * nu * S_nu.  Two independent tests per source:
  A. distance test  -- |log10(d_src / d_point)|, uses Gaia parallax only
  B. flux test      -- S required to sit at the point's (d, nuL) vs max_value
A source is a match only if BOTH agree.  A is the robust one (x-axis ticks
read cleanly); B carries the y-axis reading error, ~0.2 dex.
"""
import numpy as np
import pandas as pd

PC_CM = 3.0857e18
NU = 150e9          # Hz
MJY = 1e-26         # erg s^-1 cm^-2 Hz^-1

# (label, distance_pc, nuLnu_erg_s) read off the figure axes
PLOT_POINTS = [
    ("red#1 (short)",   971, 1.05e31),
    ("red#2 (short)",  1496, 3.50e31),
    ("green (long)",   3048, 2.40e31),
]


def nu_L_nu(d_pc, s_mjy):
    return 4 * np.pi * (d_pc * PC_CM) ** 2 * NU * s_mjy * MJY


def s_required(d_pc, nuL):
    return nuL / (4 * np.pi * (d_pc * PC_CM) ** 2 * NU * MJY)


def load():
    cat = pd.read_csv("data/centroids 4_yr.cat", sep="|", skipinitialspace=True)
    cat.columns = [c.strip() for c in cat.columns]
    cat = cat.loc[:, ~cat.columns.str.contains("^Unnamed")]
    cat["id"] = cat["id"].str.strip()

    gaia = pd.read_csv("outputs/gaia_official_10arcsec.csv")
    gaia = gaia.sort_values("p_chance").groupby("spt_id", as_index=False).first()

    cls = pd.read_excel("data/gp_long_transient_match.xlsx")
    cls = cls.rename(columns={"Unnamed: 6": "subtype"})[
        ["name", "classification", "subtype"]
    ].drop_duplicates("name")

    t = cat.merge(gaia, left_on="id", right_on="spt_id", how="left")
    t = t.merge(cls, left_on="id", right_on="name", how="left")
    t["dist_pc"] = np.where(t.parallax > 0, 1000.0 / t.parallax, np.nan)
    t["nuLnu"] = nu_L_nu(t.dist_pc, t.max_value.abs())
    return t


def main():
    t = load()
    # Galactic-field points must be confirmed Galactic, not AGN candidates
    galactic = ~t.classification.astype(str).str.contains(
        "AGN|artifact|Blazar|nan", case=False, na=True
    )
    pd.set_option("display.width", 250)

    for label, d0, l0 in PLOT_POINTS:
        s_need = s_required(d0, l0)
        print(f"\n=== {label}:  d={d0} pc, nuL={l0:.2e}  ->  needs S150 = {s_need:.1f} mJy"
              f"  (plx = {1000/d0:.3f} mas) ===")
        c = t[t.dist_pc.notna()].copy()
        c["d_test"] = np.abs(np.log10(c.dist_pc / d0))
        c["flux_ratio"] = c.max_value.abs() / s_need
        c["galactic"] = galactic[c.index]
        c = c.nsmallest(4, "d_test")
        print(c[["id", "parallax", "parallax_over_error", "dist_pc", "d_test",
                 "max_value", "flux_ratio", "nuLnu", "classification", "subtype",
                 "galactic"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\n=== all confirmed-Galactic sources, with and without a Gaia distance ===")
    g = t[galactic][["id", "classification", "subtype", "max_value",
                     "parallax", "parallax_over_error", "dist_pc", "nuLnu"]]
    print(g.sort_values("dist_pc").to_string(index=False, float_format=lambda x: f"{x:.4g}"))


if __name__ == "__main__":
    main()
