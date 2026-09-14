#!/usr/bin/env python
"""SPT-3G 150 GHz peak flux vs AllWISE W1 flux — one point per source (73).

AGN-vs-Galactic framing: sources with a reliable Gaia counterpart
(official p_chance < 0.02, adopted 2026-08-03) are Galactic stars; the
rest have no reliable optical counterpart. Key sources annotated.

Units (no raw map units on axes):
- SPT: `max_value` from `centroids 4_yr.cat` = photutils peak in the
  discovery-year yearly map, mJy per the Wan+2025 detection pipeline
  (matched-filtered maps are mJy-calibrated; cross-checked against
  Wan+2025 Table: J173508 daily peak 86.6 mJy in 2024 vs yearly-diluted
  8.52 here). Flux scale pending confirmation with Yujie.
- WISE: W1 Vega mag -> mJy with the official zero point F0(W1)=309.540 Jy
  (WISE Explanatory Supplement, same constant as src/units.py); right
  axis shows the equivalent Vega mag.

Input : outputs/v4_crossmatch_table.csv, outputs/gaia_official_10arcsec.csv
Output: outputs/images/spt_vs_wise_flux.png
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))
from units import UNWISE_F0_JY

V4  = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
GA  = os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv")
OUT = os.path.join(REPO, "outputs", "images", "spt_vs_wise_flux.png")

BLUE, ORANGE, GRAY = "#3b6fb6", "#d97b29", "#9aa0a8"
F0_W1_MJY = UNWISE_F0_JY["W1"] * 1e3          # 309540 mJy, Vega zero point

KEY = {
    "SPT3G_J171736.0-334209.1": "TXS 1714-336 (blazar)",
    "SPT3G_J171310.1-341825.7": "J171310 (QSO?)",
    "SPT3G_J173508.4-293000.8": "J173508 (Wan+2025)",
}

v4 = pd.read_csv(V4)
gaia = pd.read_csv(GA)

# official per-source reliability: min p_chance over all pairs within 10"
pmin = gaia.groupby("spt_id")["p_chance"].min()
v4["gaia_reliable"] = v4["id"].map(pmin).lt(0.02).fillna(False)

v4["w1_flux_mJy"] = F0_W1_MJY * 10 ** (-v4["AllWISE_flux_mag"] / 2.5)
v4["close"] = v4["AllWISE_sep_arcsec"] <= 10.0   # beyond 10": likely chance

fig, ax = plt.subplots(figsize=(7.5, 5.8))
for rel, color in ((True, BLUE), (False, GRAY)):
    for close, fc in ((True, color), (False, "none")):
        sub = v4[(v4.gaia_reliable == rel) & (v4.close == close)]
        ax.scatter(sub.max_value, sub.w1_flux_mJy, s=34, facecolors=fc,
                   edgecolors=color, linewidths=1.2,
                   alpha=0.9 if close else 0.7)

from matplotlib.lines import Line2D
handles = []
for rel, color, gtxt in ((True, BLUE, "Gaia star (p < 0.02)"),
                         (False, GRAY, "no reliable Gaia")):
    for close, fc, wtxt in ((True, None, "W1 ≤ 10″"),
                            (False, "none", "W1 10–30″")):
        n = int(((v4.gaia_reliable == rel) & (v4.close == close)).sum())
        handles.append(Line2D([], [], ls="", marker="o",
                              mfc=color if fc is None else "none", mec=color,
                              label=f"{gtxt}, {wtxt}, n={n}"))
ax.legend(handles=handles, fontsize=8, loc="upper right", frameon=False)

KEY_ANNOTATIONS = False  # per-source annotations disabled (2026-08-10 request)
if KEY_ANNOTATIONS:
    for sid, lab in KEY.items():
        row = v4[v4.id == sid]
        if len(row):
            r = row.iloc[0]
            ax.scatter(r.max_value, r.w1_flux_mJy, s=90, facecolors="none",
                       edgecolors=ORANGE, linewidths=1.8, zorder=3)
            right_edge = r.max_value > 100
            ax.annotate(lab, (r.max_value, r.w1_flux_mJy),
                        textcoords="offset points",
                        xytext=(-7, 5) if right_edge else (7, 5),
                        ha="right" if right_edge else "left",
                        fontsize=8, color="0.15")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("SPT-3G 150 GHz peak flux, discovery-year map (mJy)")
ax.set_ylabel("AllWISE W1 flux (mJy, Vega-based)")

mag_axis = ax.secondary_yaxis(
    "right",
    functions=(lambda f: -2.5 * np.log10(np.maximum(f, 1e-12) / F0_W1_MJY),
               lambda m: F0_W1_MJY * 10 ** (-m / 2.5)))
from matplotlib.ticker import FixedLocator, NullLocator, FormatStrFormatter
mag_axis.yaxis.set_major_locator(FixedLocator(np.arange(0, 16, 2)))
mag_axis.yaxis.set_minor_locator(NullLocator())
mag_axis.yaxis.set_major_formatter(FormatStrFormatter("%g"))
mag_axis.set_ylabel("W1 (Vega mag)")

ax.set_title("SPT 150 GHz vs WISE W1 flux — 73 SPT-3G variables", fontsize=12)
ax.text(0.99, 0.01,
        "open markers: W1 match at 10–30″ (chance alignment likely)\n"
        "SPT peak per Wan+2025 pipeline (yearly matched-filtered map)",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=7, color="0.45")
for s in ("top",):
    ax.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(OUT, dpi=200)
print("wrote", OUT)
