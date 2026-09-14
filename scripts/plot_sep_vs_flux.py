#!/usr/bin/env python
"""Match separation vs catalog flux — one standalone figure per catalog,
plus a 2x2 WISE W1-W4 figure and a Gaia all-pairs figure.

Diagnostic for counterpart reliability: genuine counterparts sit in the
bright/close corner; chance alignments are faint and spread in separation.
Brightness increases to the RIGHT in every panel: radio fluxes (mJy, span
~4 dex) get a log axis; magnitude axes are linear-in-mag but inverted —
mag is already logarithmic in flux, so a linear mag axis IS a log-flux axis.

Input : outputs/v4_crossmatch_table.csv, outputs/gaia_official_10arcsec.csv
Output: outputs/images/sep_vs_flux/sep_vs_flux_<catalog>.png
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V4     = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
GA     = os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv")
OUTDIR = os.path.join(REPO, "outputs", "images", "sep_vs_flux")
os.makedirs(OUTDIR, exist_ok=True)

BLUE, ORANGE, GRAY = "#3b6fb6", "#d97b29", "#9aa0a8"

KEY = {  # sources to annotate wherever they appear
    "SPT3G_J171736.0-334209.1": "TXS 1714-336",
    "SPT3G_J171310.1-341825.7": "J171310",
    "SPT3G_J170919.0-352512.8": "J170919",
    "SPT3G_J173508.4-293000.8": "J173508",
}

v4 = pd.read_csv(V4)
med_sigma = v4["sigma_pos_arcsec_heuristic"].median()


def style(ax, n):
    # thresholds don't depend on catalog flux -> horizontal lines at the
    # sample-median sigma_pos (per-source values span 6.1-11.6").
    for k, ls in ((1, "--"), (2, "-"), (3, ":")):
        ax.axhline(k * med_sigma, color=ORANGE, ls=ls, lw=1.5,
                   label=f"{k}×median σ_pos")
    ax.set_yscale("log")
    ax.set_ylabel("angular separation (arcsec)")
    ax.text(0.03, 0.05, f"n = {n}", transform=ax.transAxes,
            fontsize=9, va="bottom", color="0.3")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def annotate(ax, sub, xcol, scale=1.0):
    return  # per-source annotations disabled (2026-08-03 request)
    for sid, lab in KEY.items():
        row = sub[sub.id == sid]
        if len(row):
            r = row.iloc[0]
            ax.scatter(r[xcol] * scale, r["sep"], s=60, facecolors="none",
                       edgecolors=ORANGE, linewidths=1.8, zorder=3)
            ax.annotate(lab, (r[xcol] * scale, r["sep"]),
                        textcoords="offset points", xytext=(6, 6),
                        fontsize=8, color="0.15")


# ── one standalone figure per radio catalog (log flux axis, mJy) ──
RADIO = [  # (title, cat prefix, flux column, to-mJy factor)
    ("RACS-mid 1.4 GHz", "RACS-mid", "RACS-mid_flux_mJy", 1.0),
    ("VLASS 3 GHz",      "VLASS",    "VLASS_flux_Jy",     1000.0),
    ("NVSS 1.4 GHz",     "NVSS",     "NVSS_flux_mJy",     1.0),
    ("AT20G 20 GHz",     "AT20G",    "AT20G_flux_Jy",     1000.0),
]
for title, cat, fcol, scale in RADIO:
    sub = v4[(v4[f"{cat}_match"] == True) & v4[fcol].notna()].copy()
    sub["sep"] = sub[f"{cat}_sep_arcsec"]
    fig, ax = plt.subplots(figsize=(7, 5.4))
    ax.scatter(sub[fcol] * scale, sub["sep"], s=30, c=BLUE, alpha=0.85,
               edgecolors="none")
    annotate(ax, sub, fcol, scale)
    ax.set_xscale("log")
    ax.set_xlabel("flux density (mJy)")
    ax.set_title(f"{title}: separation vs flux", fontsize=12)
    style(ax, len(sub))
    fig.tight_layout()
    out = os.path.join(OUTDIR, f"sep_vs_flux_{cat}.png")
    fig.savefig(out, dpi=200); plt.close(fig)
    print("wrote", out)

# ── WISE W1-W4, 2x2 in one figure (mag axes: linear-in-mag = log-in-flux,
#     inverted so brightness increases rightward) ──
WISE = [("W1 (3.4 μm)", "AllWISE_flux_mag"),
        ("W2 (4.6 μm)", "AllWISE_W2mag"),
        ("W3 (12 μm)",  "AllWISE_W3mag"),
        ("W4 (22 μm)",  "AllWISE_W4mag")]
fig, axes = plt.subplots(2, 2, figsize=(12, 9.5))
for ax, (label, col) in zip(axes.flat, WISE):
    sub = v4[(v4["AllWISE_match"] == True) & v4[col].notna()].copy()
    sub["sep"] = sub["AllWISE_sep_arcsec"]
    ax.scatter(sub[col], sub["sep"], s=28, c=BLUE, alpha=0.85, edgecolors="none")
    annotate(ax, sub, col)
    ax.invert_xaxis()
    ax.set_xlabel(f"{label.split()[0]} (mag, Vega) — brighter →")
    ax.set_title(f"AllWISE {label}", fontsize=11)
    style(ax, len(sub))
fig.suptitle("AllWISE W1–W4: separation vs brightness "
             "(same match, four bands; W3/W4 mags may include upper limits)",
             fontsize=12, y=0.995)
fig.tight_layout()
out = os.path.join(OUTDIR, "sep_vs_flux_WISE_W1-W4.png")
fig.savefig(out, dpi=200); plt.close(fig)
print("wrote", out)

# ── Gaia: all 792 pairs, reliability from official p_chance ──
g = pd.read_csv(GA)
fig, ax = plt.subplots(figsize=(7.5, 5.6))
rel, oth = g[g.p_chance < 0.02], g[g.p_chance >= 0.02]
ax.scatter(oth.phot_g_mean_mag, oth.gaia_sep_arcsec, s=14, c=GRAY, alpha=0.5,
           edgecolors="none", label="p_chance ≥ 0.02")
ax.scatter(rel.phot_g_mean_mag, rel.gaia_sep_arcsec, s=26, c=BLUE, alpha=0.9,
           edgecolors="none", label="p_chance < 0.02")
ax.invert_xaxis()
ax.set_xlabel("Gaia G (mag) — brighter →")
ax.set_title("Gaia DR3 (all pairs within 10″): separation vs brightness",
             fontsize=12)
ax.legend(fontsize=8, loc="lower left", frameon=False)
style(ax, len(g))
fig.tight_layout()
out = os.path.join(OUTDIR, "sep_vs_flux_Gaia.png")
fig.savefig(out, dpi=200); plt.close(fig)
print("wrote", out)
