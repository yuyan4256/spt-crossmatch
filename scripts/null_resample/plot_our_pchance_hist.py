#!/usr/bin/env python
"""p_chance histogram of our own matched sources, in the layout used by
Tandoi et al. 2024 Fig. 6 (log p, log count, dotted line at the 0.02 cut).

One entry per SPT source that has a Gaia DR3 star within 10" (63 of 73),
taking that source's best match. p comes from this field's own null
(400k random positions per event-TS curve, floors 3.05"/4.42").

Writes outputs/images/our_pchance_hist.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SRC = os.path.join(REPO, "outputs", "gaia_pchance_floors.csv")
OUT = os.path.join(REPO, "outputs", "images", "our_pchance_hist.png")

d = pd.read_csv(SRC)
p = d.loc[d.groupby("spt_id")["association_TS_F0_abs"].idxmax(),
          "p_F0_abs"].values
n = len(p)
FLOOR = 2.5e-6
p = np.clip(p, FLOOR, 1.0)

plt.rcParams.update({"font.family": "serif", "font.size": 12,
                     "mathtext.fontset": "cm"})
bins = np.logspace(np.log10(1e-3), 0, 19)
fig, ax = plt.subplots(figsize=(6.4, 4.8))
ax.hist(p, bins=bins, histtype="step", color="k", lw=1.6)
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
print("saved", OUT, f"| {n} sources, {(p < 0.02).sum()} below 0.02, "
      f"min p = {p.min():.4f}")
