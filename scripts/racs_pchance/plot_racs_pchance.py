#!/usr/bin/env python
"""Summary figure for the RACS-mid-position Gaia p_chance re-run of the
AGN-marked sources.

Left : empirical CDF of p_chance against the uniform expectation -- if none of
       the AGN candidates has a real Gaia counterpart the p-values must be
       uniform on [0, 1].
Right: p_chance from the SPT centroid (Galactic null) vs from the RACS-mid
       position, same null family, for the sources present in both runs.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
W = os.path.join(REPO, "outputs", "racs_agn_pchance_wide.csv")
OUT = os.path.join(REPO, "outputs", "images", "racs_agn_pchance.png")

w = pd.read_csv(W)
p = np.sort(w["p_chance_racs_A_fit_plus_floor"].values)
n = len(p)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))

ax[0].step(np.concatenate(([0], p)), np.arange(n + 1) / n, where="post",
           color="C0", lw=2, label=f"RACS-mid position ({n} AGN-marked)")
ax[0].plot([0, 1], [0, 1], "k--", lw=1, label="uniform (no counterpart)")
ax[0].set_xlabel(r"$p_{\rm chance}$ (Gaia DR3, Galactic-plane null)")
ax[0].set_ylabel("cumulative fraction of sources")
ax[0].set_xlim(0, 1)
ax[0].set_ylim(0, 1)
ax[0].legend(loc="lower right", fontsize=9)
ax[0].set_title(r"$\sigma = \sqrt{\sigma_{\rm fit}^2+\sigma_{\rm floor}^2}$,"
                " floor 0.65''/0.97'' (Duchesne+2024)", fontsize=9)

m = w["p_galactic_spt_pos"].notna()
sc = ax[1].scatter(w.loc[m, "p_galactic_spt_pos"],
                   w.loc[m, "p_chance_racs_A_fit_plus_floor"],
                   c=w.loc[m, "gaia_sep_arcsec"], cmap="viridis", s=45,
                   edgecolor="k", linewidth=0.4)
ax[1].plot([1e-3, 1], [1e-3, 1], "k--", lw=1)
for thr in (0.02,):
    ax[1].axhline(thr, color="crimson", lw=0.8, ls=":")
    ax[1].axvline(thr, color="crimson", lw=0.8, ls=":")
ax[1].set_xscale("log")
ax[1].set_yscale("log")
ax[1].set_xlabel(r"$p_{\rm chance}$ from SPT centroid")
ax[1].set_ylabel(r"$p_{\rm chance}$ from RACS-mid position")
plt.colorbar(sc, ax=ax[1], label="Gaia separation from RACS position (arcsec)")
for _, r in w.loc[m].nsmallest(2, "p_chance_racs_A_fit_plus_floor").iterrows():
    ax[1].annotate(r["name"].replace("SPT3G_", ""),
                   (r["p_galactic_spt_pos"],
                    r["p_chance_racs_A_fit_plus_floor"]),
                   textcoords="offset points", xytext=(-8, -10), fontsize=7,
                   ha="right")

fig.tight_layout()
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT, dpi=160)
print("saved", OUT)
