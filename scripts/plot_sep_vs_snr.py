#!/usr/bin/env python
"""Association diagnostic in the style of the SPT-3G x AT20G thesis figure:
match separation vs |SNR|, with the positional-uncertainty model overlaid.

Curves: radial sigma_pos(SNR) from the Wan+2025 (arXiv:2509.08962) model,
calibrated on this Galactic-plane field (Sec 4.4, Eqs 4-6; sigma_sys,rel
per-event term pending, 0 for now):
    sigma_RA  = sqrt(3.05^2 + (50/|SNR|)^2)   [arcsec]
    sigma_Dec = sqrt(4.42^2 + (49/|SNR|)^2)
    sigma_r   = sqrt(sigma_RA^2 + sigma_Dec^2)
drawn at 1x (dashed), 2x (solid) and 3x (dotted). Matches above the 3x
curve are likely chance alignments; the curves flatten at high SNR where
the pointing floor dominates — same behavior as the reference figure.

Input : outputs/v4_crossmatch_table.csv, outputs/gaia_official_10arcsec.csv
Output: outputs/images/sep_vs_snr/sep_vs_snr_<catalog>.png
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V4     = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
GA     = os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv")
OUTDIR = os.path.join(REPO, "outputs", "images", "sep_vs_snr")
os.makedirs(OUTDIR, exist_ok=True)

BLUE, ORANGE, GRAY = "#3b6fb6", "#d97b29", "#9aa0a8"

KEY = {
    "SPT3G_J171736.0-334209.1": "TXS 1714-336",
    "SPT3G_J171310.1-341825.7": "J171310",
    "SPT3G_J170919.0-352512.8": "J170919",
    "SPT3G_J173508.4-293000.8": "J173508",
}


def sigma_r(snr):
    """Radial 1-sigma positional uncertainty (arcsec), Wan+2025 Eqs 4-6."""
    s = np.abs(snr)
    return np.sqrt(3.05**2 + (50.0/s)**2 + 4.42**2 + (49.0/s)**2)


def draw(ax, x, y, n, title, xlab):
    grid = np.logspace(np.log10(5), np.log10(300), 200)
    for k, ls in ((1, "--"), (2, "-"), (3, ":")):
        ax.plot(grid, k*sigma_r(grid), color=ORANGE, ls=ls, lw=1.5,
                label=f"{k}σ_pos(SNR)")
    ax.scatter(x, y, s=26, c=BLUE, alpha=0.85, edgecolors="none", zorder=3)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(5, 300)
    ax.set_xlabel(xlab)
    ax.set_ylabel("angular separation (arcsec)")
    ax.set_title(title, fontsize=12)
    ax.text(0.03, 0.05, f"n = {n}", transform=ax.transAxes, fontsize=9,
            color="0.3")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def annotate(ax, sub, xcol, ycol):
    return  # per-source annotations disabled (2026-08-03 request)
    for sid, lab in KEY.items():
        row = sub[sub["id"] == sid]
        if len(row):
            r = row.iloc[0]
            ax.scatter(r[xcol], r[ycol], s=60, facecolors="none",
                       edgecolors=ORANGE, linewidths=1.8, zorder=4)
            ax.annotate(lab, (r[xcol], r[ycol]), textcoords="offset points",
                        xytext=(6, 5), fontsize=8, color="0.15")


v4 = pd.read_csv(V4)
v4["abs_snr"] = v4["snr_max"].abs()

CATS = ["RACS-mid", "VLASS", "NVSS", "AT20G", "AllWISE"]
for cat in CATS:
    sub = v4[(v4[f"{cat}_match"] == True) & v4[f"{cat}_sep_arcsec"].notna()].copy()
    fig, ax = plt.subplots(figsize=(7, 5.6))
    draw(ax, sub["abs_snr"], sub[f"{cat}_sep_arcsec"], len(sub),
         f"SPT-3G × {cat} associations",
         "SPT-3G |snr_max| (4-yr variability)")
    annotate(ax, sub, "abs_snr", f"{cat}_sep_arcsec")
    fig.tight_layout()
    out = os.path.join(OUTDIR, f"sep_vs_snr_{cat}.png")
    fig.savefig(out, dpi=200); plt.close(fig)
    print("wrote", out)

# Gaia: all pairs, colored by official reliability
g = pd.read_csv(GA)
g["abs_snr"] = g["spt_snr_max"].abs()
g["id"] = g["spt_id"]
fig, ax = plt.subplots(figsize=(7.5, 5.8))
rel, oth = g[g.p_chance < 0.02], g[g.p_chance >= 0.02]
ax.scatter(oth.abs_snr, oth.gaia_sep_arcsec, s=12, c=GRAY, alpha=0.45,
           edgecolors="none", label="p_chance ≥ 0.02", zorder=2)
draw(ax, rel.abs_snr, rel.gaia_sep_arcsec, len(g),
     "SPT-3G × Gaia DR3 associations (all pairs within 10″)",
     "SPT-3G |snr_max| (4-yr variability)")
best = g.sort_values("p_chance").groupby("spt_id").first().reset_index()
best["id"] = best["spt_id"]
annotate(ax, best, "abs_snr", "gaia_sep_arcsec")
h, l = ax.get_legend_handles_labels()
ax.legend(h, ["p_chance ≥ 0.02" if x == "p_chance ≥ 0.02" else x for x in l],
          fontsize=8, loc="lower left", frameon=False)
fig.tight_layout()
out = os.path.join(OUTDIR, "sep_vs_snr_Gaia.png")
fig.savefig(out, dpi=200); plt.close(fig)
print("wrote", out)
