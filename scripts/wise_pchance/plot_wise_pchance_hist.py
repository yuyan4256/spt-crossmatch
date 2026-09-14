#!/usr/bin/env python
"""p_chance histogram for the AllWISE run, drawn on the same log-log axes as
outputs/images/our_pchance_hist.png (Tandoi+2024 Fig. 6 layout) with the
Galactic-field Gaia result overlaid for comparison.

Both histograms come from the identical machinery and the same field null --
only the counterpart catalog differs (AllWISE ranked by W1, 400k null
positions; Gaia DR3 G<=20, 250k per event-TS curve), floors 3.05"/4.42".
The dashed line is what a sample with no real counterpart must give
(p uniform on [0,1]), scaled to the WISE sample size.

Writes outputs/images/wise_pchance_hist.png
"""
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WISE = os.path.join(REPO, "outputs", "wise_pchance.csv")
GAIA = os.path.join(REPO, "outputs", "gaia_pchance_galactic 3.05'' 4.42''.csv")
OUT = os.path.join(REPO, "outputs", "images", "wise_pchance_hist.png")

FLOOR = 2.5e-6   # 1 / 400000 null draws: nothing is resolvable below this
LO = 1e-5

w = pd.read_csv(WISE)["p_wise_gal"].dropna().to_numpy()
g = pd.read_csv(GAIA)
g = g.loc[g.groupby("spt_id")["association_TS"].idxmax(), "p_galactic"].to_numpy()

plt.rcParams.update({"font.family": "serif", "font.size": 12,
                     "mathtext.fontset": "cm"})
bins = np.logspace(np.log10(LO), 0, 25)
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ax.hist(np.clip(w, max(FLOOR, LO * 1.05), 1), bins=bins, histtype="step",
        color="#b03030", lw=1.8, label=f"AllWISE, W1  ({len(w)} sources)")
ax.hist(np.clip(g, max(FLOOR, LO * 1.05), 1), bins=bins, histtype="step",
        color="#3b6fb6", lw=1.4, ls="-", label=f"Gaia DR3, $G$  ({len(g)} sources)")
ax.step(bins[1:], len(w) * np.diff(bins), where="pre", color="0.45", ls="--",
        lw=1.2, label="uniform (no real counterpart)")
ax.axvline(0.02, color="k", ls=":", lw=1.2)
ax.text(0.023, 40, "0.02", fontsize=9, rotation=90, va="top")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(LO, 1)
ax.set_ylim(0.7, 60)
ax.set_xlabel("chance probability $p$")
ax.set_ylabel("count")
ax.legend(fontsize=9, loc="upper left", frameon=False)
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=180)
print("saved", OUT)
for c in (0.01, 0.02, 0.05):
    print(f"  WISE p<{c}: {int((w < c).sum())}/{len(w)} (expected {c*len(w):.1f})"
          f" | Gaia {int((g < c).sum())}/{len(g)} (expected {c*len(g):.1f})")
print(f"  min p: WISE {w.min():.2e}, Gaia {g.min():.2e}")
