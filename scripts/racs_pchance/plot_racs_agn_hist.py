#!/usr/bin/env python
"""p_chance histogram for the radio-anchored search on the AGN candidates,
drawn on exactly the same axes as outputs/images/our_pchance_hist.png so the
two can be put side by side.

One entry per AGN-marked source (45), scored at its RACS-mid position against
this field's Gaia null (sigma = fit (+) floor 0.65"/0.97", Duchesne+2024).
Run with --uniform to overlay what a sample with no stellar counterpart at
all must give (p uniform on [0, 1]); the default is the plain histogram, in
the same format as the SPT-centroid one.

Writes outputs/images/racs_agn_pchance_hist.png
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "outputs", "racs_agn_pchance_wide.csv")
OUT = os.path.join(REPO, "outputs", "images", "racs_agn_pchance_hist.png")

p = pd.read_csv(SRC)["p_chance_racs_A_fit_plus_floor"].values
n = len(p)

plt.rcParams.update({"font.family": "serif", "font.size": 12,
                     "mathtext.fontset": "cm"})
bins = np.logspace(-3, 0, 19)
fig, ax = plt.subplots(figsize=(6.4, 4.8))
ax.hist(np.clip(p, 1.1e-3, 1), bins=bins, histtype="step", color="k", lw=1.6,
        label=f"{n} AGN candidates, RACS-mid position")
if "--uniform" in sys.argv:
    ax.step(bins[1:], n * np.diff(bins), where="pre", color="0.45", ls="--",
            lw=1.2, label="uniform (no stellar counterpart)")
    ax.legend(fontsize=9, loc="upper left", frameon=False)
ax.axvline(0.02, color="k", ls=":", lw=1.2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-3, 1)
ax.set_ylim(0.7, 40)
ax.set_xlabel("p")
ax.set_ylabel("count")
fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=180)
print("saved", OUT)
for c in (0.01, 0.02, 0.05):
    print(f"  p<{c}: {int((p < c).sum())}/{n}  (expected {c * n:.1f})")
print(f"  min p = {p.min():.4f}")
