#!/usr/bin/env python
"""Histogram of best-match Gaia p_chance for the 4-yr SPT variable sources.

One entry per SPT source with >=1 Gaia DR3 match within 10" (63 of 73):
the lowest p_chance among its candidates, i.e. the same "best match per
source" convention as plot_gaia_verdict.py. Log-spaced bins; the reliable-
counterpart threshold p_chance < 0.02 splits the coloring.

Input : outputs/gaia_official_10arcsec.csv, outputs/v4_crossmatch_table.csv
Output: outputs/images/gaia_pchance_hist.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
G_CSV = os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv")
V4_CSV = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
OUT = os.path.join(REPO, "outputs", "images", "gaia_pchance_hist.png")

BLUE, GRAY = "#3b6fb6", "#8a8f98"
P_REL = 0.02

g = pd.read_csv(G_CSV)
v4 = pd.read_csv(V4_CSV)
best = g.sort_values("p_chance").groupby("spt_id").first().reset_index()

n_rel = int((best.p_chance < P_REL).sum())
n_unrel = len(best) - n_rel
n_nomatch = len(v4) - len(best)

fig, ax = plt.subplots(figsize=(7.2, 4.6))

bins = np.logspace(np.log10(best.p_chance.min()) - 0.1, np.log10(0.5), 26)
# split the bin edges exactly at the threshold so the two colors never share a bin
bins = np.sort(np.append(bins[~np.isclose(bins, P_REL)], P_REL))

ax.hist(best.p_chance[best.p_chance < P_REL], bins=bins, color=BLUE,
        edgecolor="white", linewidth=0.6)
ax.hist(best.p_chance[best.p_chance >= P_REL], bins=bins, color=GRAY,
        edgecolor="white", linewidth=0.6)

ax.set_ylim(0, ax.get_ylim()[1] * 1.18)

ax.axvline(P_REL, color="0.25", ls=":", lw=1.2)
ax.text(P_REL, ax.get_ylim()[1] * 0.99, " $p_{\\rm chance}=0.02$",
        fontsize=9, color="0.25", ha="left", va="top")

ax.text(0.03, 0.86,
        f"reliable counterpart  $p<0.02$:  {n_rel} sources",
        transform=ax.transAxes, fontsize=10, color=BLUE, fontweight="bold")
ax.text(0.03, 0.78,
        f"ambiguous  $p\\geq 0.02$:  {n_unrel} sources",
        transform=ax.transAxes, fontsize=10, color=GRAY, fontweight="bold")

ax.set_xscale("log")
ax.set_xlabel("best-match $p_{\\rm chance}$ (SPT-3G official, TS + MC table)")
ax.set_ylabel("number of SPT sources")
ax.set_title(f"Gaia DR3 chance-association probability "
             f"({len(best)}/{len(v4)} sources with a Gaia match within 10\")",
             fontsize=11)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
print(f"best-match sources: {len(best)}  (no Gaia in 10\": {n_nomatch})")
print(f"p<{P_REL}: {n_rel}   p>={P_REL}: {n_unrel}")
