#!/usr/bin/env python
"""Summary figure for the official-Gaia crossmatch report.

Left panel: selection funnel 73 -> 63 (>=1 Gaia in 10") -> 42 (best p_chance
< 0.02) -> 29 (astrometric motion >=5 sigma, i.e. Galactic stars).
Right panel: best-match p_chance vs Gaia G, marker filled if the match has
significant astrometry (parallax_over_error>5 or per-component pm/err>5 —
raw Gaia DR3 columns only). Key sources annotated.

Input : outputs/gaia_official_10arcsec.csv, outputs/v4_crossmatch_table.csv
Output: outputs/images/gaia_verdict_summary.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G_CSV = os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv")
V4_CSV = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
OUT = os.path.join(REPO, "outputs", "images", "gaia_verdict_summary.png")

BLUE, ORANGE, GRAY = "#3b6fb6", "#d97b29", "#8a8f98"

g = pd.read_csv(G_CSV)
v4 = pd.read_csv(V4_CSV)
g["pmra_snr"] = (g.pmra / g.pmra_error).abs()
g["pmdec_snr"] = (g.pmdec / g.pmdec_error).abs()
g["moving"] = (g.parallax_over_error.abs() > 5) | (g.pmra_snr > 5) | (g.pmdec_snr > 5)

best = g.sort_values("p_chance").groupby("spt_id").first().reset_index()
rel = best[best.p_chance < 0.02]
n_all, n_gaia, n_rel = len(v4), len(best), len(rel)
n_star = int(rel.moving.sum())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2), width_ratios=[1, 1.35])

# ---- funnel ----
steps = [
    (f"SPT 4-yr variable sources", n_all),
    (f'$\\geq$1 Gaia DR3 within 10"', n_gaia),
    ("best $p_{\\rm chance}<0.02$\n(reliable counterpart)", n_rel),
    ("parallax or PM $\\geq 5\\sigma$\n($\\Rightarrow$ Galactic star)", n_star),
]
y = np.arange(len(steps))[::-1]
ax1.barh(y, [s[1] for s in steps], height=0.62, color=[GRAY, GRAY, BLUE, BLUE])
for yi, (lab, n) in zip(y, steps):
    ax1.text(n + 1.2, yi, str(n), va="center", fontsize=12, fontweight="bold")
    ax1.text(1.5, yi, lab, va="center", fontsize=9.5, color="white")
ax1.set_xlim(0, 82)
ax1.set_yticks([])
ax1.set_xlabel("number of SPT sources")
ax1.set_title("Gaia DR3 official crossmatch: selection funnel", fontsize=11)
for s in ("top", "right", "left"):
    ax1.spines[s].set_visible(False)

# ---- p_chance vs G ----
mv, st = best[best.moving], best[~best.moving]
ax2.scatter(mv.p_chance, mv.phot_g_mean_mag, s=42, c=BLUE, edgecolors="none",
            alpha=0.85, label="parallax or PM $\\geq5\\sigma$ (star)")
ax2.scatter(st.p_chance, st.phot_g_mean_mag, s=48, facecolors="none",
            edgecolors=ORANGE, linewidths=1.6,
            label="no significant astrometry")
ax2.axvline(0.02, color="0.35", ls=":", lw=1)
ax2.text(0.017, 10.9, "$p_{\\rm chance}=0.02$ ", fontsize=9, color="0.35",
         ha="right", va="top")

marks = {
    "SPT3G_J173508.4-293000.8": ("J173508", (-10, -14)),
    "SPT3G_J171736.0-334209.1": ("TXS 1714-336", (-10, 14)),
    "SPT3G_J171310.1-341825.7": ("J171310", (10, 10)),
    "SPT3G_J170919.0-352512.8": ("J170919", (8, -22)),
}
for sid, (lab, off) in marks.items():
    row = best[best.spt_id == sid]
    if len(row):
        r = row.iloc[0]
        ax2.annotate(lab, (r.p_chance, r.phot_g_mean_mag), textcoords="offset points",
                     xytext=off, fontsize=8, color="0.15",
                     arrowprops=dict(arrowstyle="-", lw=0.7, color="0.5"))

ax2.set_xscale("log")
ax2.invert_yaxis()
ax2.set_xlabel("best-match $p_{\\rm chance}$ (SPT-3G official, TS + MC table)")
ax2.set_ylabel("Gaia $G$ (mag)")
ax2.set_title("Best Gaia match per source: chance probability vs brightness",
              fontsize=11)
ax2.legend(loc="lower left", fontsize=8.5, frameon=False)
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
print(f"funnel: {n_all} -> {n_gaia} -> {n_rel} -> {n_star}")
wd = best[best.spt_id == "SPT3G_J173508.4-293000.8"].iloc[0]
print(f"J173508 p_chance = {wd.p_chance:.2e}")
