#!/usr/bin/env python
"""p_chance from AllWISE vs from Gaia, one point per source scored in both.

Same null machinery, same field, same floors (3.05"/4.42") -- the two axes are
independent counterpart catalogs, so agreement in the lower-left quadrant is a
genuine cross-validation and off-diagonal points flag sources detected in one
band only (e.g. optically extincted but infrared-bright).

Writes outputs/images/wise_vs_gaia_pchance.png
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "outputs", "wise_vs_gaia_pchance.csv")
OUT = os.path.join(REPO, "outputs", "images", "wise_vs_gaia_pchance.png")

FLOOR = 2.5e-6
LO = 1e-5
THR = 0.02
LABEL = {
    "SPT3G_J173508.4-293000.8": "J173508",
    "SPT3G_J174135.6-300634.5": "J174135",
    "SPT3G_J180600.4-281706.4": "J180600",
    "SPT3G_J180350.1-281905.2": "J180350",
}

d = pd.read_csv(SRC).dropna(subset=["p_wise_gal", "p_gaia_gal"])
x = np.clip(d["p_wise_gal"].to_numpy(), max(FLOOR, LO * 1.05), 1)
y = np.clip(d["p_gaia_gal"].to_numpy(), max(FLOOR, LO * 1.05), 1)

plt.rcParams.update({"font.family": "serif", "font.size": 12,
                     "mathtext.fontset": "cm"})
fig, ax = plt.subplots(figsize=(5.8, 5.4))
# shade only the quadrant where BOTH catalogues call the match reliable
frac = np.log10(THR / LO) / np.log10(1.0 / LO)
ax.axvspan(LO, THR, ymin=0, ymax=frac, color="#f2e2e2", zorder=0)
ax.plot([LO, 1], [LO, 1], color="0.7", lw=1.0, ls="-", zorder=1)
both = (x < THR) & (y < THR)
ax.scatter(x[~both], y[~both], s=26, facecolor="none", edgecolor="#3b6fb6",
           lw=1.1, zorder=3)
ax.scatter(x[both], y[both], s=48, color="#b03030", zorder=4)
ax.axhline(THR, color="k", ls=":", lw=1.1)
ax.axvline(THR, color="k", ls=":", lw=1.1)
for _, r in d.iterrows():
    lab = LABEL.get(r["spt_id"])
    if lab:
        ax.annotate(lab, (max(r["p_wise_gal"], LO * 1.05),
                          max(r["p_gaia_gal"], LO * 1.05)),
                    textcoords="offset points", xytext=(7, 4), fontsize=9)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(LO, 1)
ax.set_ylim(LO, 1)
ax.set_xlabel(r"$p_{\rm chance}$  (AllWISE, W1)")
ax.set_ylabel(r"$p_{\rm chance}$  (Gaia DR3, $G$)")
ax.text(THR * 0.75, LO * 1.6, "reliable in both", fontsize=9, color="#8a3b3b",
        ha="right")
ax.set_title(f"{len(d)} sources scored in both catalogues", fontsize=11)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=180)
print("saved", OUT, f"| {int(both.sum())} sources below {THR} in both")
