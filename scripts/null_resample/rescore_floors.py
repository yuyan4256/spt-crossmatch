#!/usr/bin/env python
"""Turn the multi-floor null (run_floors_rerun.py) into p_chance for the 73
sources, one column per pointing-floor variant.

Reads   data/gaia_null/ts_floors_n{N}.npz       null best-TS, one column
                                               per (floor, event TS)
        outputs/floors_source_ts.csv           per-candidate TS, all floors
Set NFLOORS to pick the null size (default 250000; the production run is
400000, matching the official null).
Writes  outputs/gaia_pchance_floors.csv        all floors, per candidate
        outputs/images/floors_pchance.png      what the extra sigma costs

p is read exactly as rescore_73_galactic.py does it: empirical survival
function on a 4000-point grid, cubic interpolation, fill (1, 0), and the
event-TS curve chosen nearest in log -- i.e. get_gaia_prob(ts_dependent=True).
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "src"))
from pchance import survival_interp
from run_floors_rerun import FLOORS, SRC_OUT, NULL_OUT

OUT_CSV = os.path.join(REPO, "outputs", "gaia_pchance_floors.csv")
OUT_PNG = os.path.join(REPO, "outputs", "images", "floors_pchance.png")
THRESH = (0.01, 0.02, 0.05)



def main():
    d = np.load(NULL_OUT)
    best, names, ev_ts = d["best_ts"], d["floor_names"], d["event_ts"]
    interps = {}
    for j in range(best.shape[1]):
        interps[(str(names[j]), float(ev_ts[j]))] = survival_interp(best[:, j])
    grid_ts = np.array(sorted({float(t) for t in ev_ts}))
    print(f"null: {best.shape[0]} events, {best.shape[1]} curves, "
          f"TS grid {grid_ts.min():g}..{grid_ts.max():g}")

    df = pd.read_csv(SRC_OUT)
    ev = df["spt_snr_max"].values ** 2
    near = grid_ts[np.argmin(np.abs(np.log(ev)[:, None] -
                                    np.log(grid_ts)[None, :]), axis=1)]
    df["event_ts"] = ev
    df["mapping_ts"] = near
    for name in FLOORS:
        ts = df[f"association_TS_{name}"].values
        df[f"p_{name}"] = [float(interps[(name, t)](a)) for t, a in zip(near, ts)]

    keys = list(FLOORS)
    best_rows = {}
    for name in keys:
        idx = df.groupby("spt_id")[f"association_TS_{name}"].idxmax()
        best_rows[name] = df.loc[idx].set_index("spt_id")

    print(f"\n{len(best_rows[keys[0]])} sources with a Gaia star within 10\"")
    print(f"{'floor':18s} {'sigma_RA/Dec':>14s} " +
          "".join(f"  p<{t:<6g}" for t in THRESH))
    for name in keys:
        b = best_rows[name]
        fr, fd = FLOORS[name]
        print(f"{name:18s} {fr:6.3f}/{fd:5.3f}\" " +
              "".join(f"{int((b['p_' + name] < t).sum()):9d}" for t in THRESH))

    prod = "p_galactic"
    if prod in df.columns:
        pb = df.loc[df.groupby("spt_id")["association_TS"].idxmax()] \
               .set_index("spt_id")
        floor_label = '3.050/4.420"'   # no backslash inside the f-string (py<3.12)
        print(f"{'production (F0)':18s} {floor_label:>14s} " +
              "".join(f"{int((pb[prod] < t).sum()):9d}" for t in THRESH))

    print("\nsources at p<0.05 in any variant:")
    ids = sorted(set().union(*[set(best_rows[n].index[best_rows[n][f"p_{n}"]
                                                      < 0.05]) for n in keys]))
    for sid in ids:
        row = " ".join(f"{n.split('_')[0]} TS={best_rows[n].loc[sid, f'association_TS_{n}']:5.2f} "
                       f"p={best_rows[n].loc[sid, f'p_{n}']:.4f}" for n in keys)
        print(f"  {sid}  {row}")

    df.to_csv(OUT_CSV, index=False)
    print("\nsaved", OUT_CSV)

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.5))
    f0 = best_rows[keys[0]]
    FLOOR_P = 1e-4
    palette = ("#C44E52", "#4C72B0", "#1F6F72", "#9B59B6")
    marks = ("o", "s", "^", "D")
    for name, c, m in zip(keys[1:], palette, marks):
        b = best_rows[name]
        ax[0].scatter(np.clip(f0[f"p_{keys[0]}"], FLOOR_P, 1),
                      np.clip(b.loc[f0.index, f"p_{name}"], FLOOR_P, 1),
                      s=42, marker=m, color=c, edgecolor="k", linewidth=0.4,
                      alpha=0.85,
                      label=f"{name} ({FLOORS[name][0]:.2f}\"/{FLOORS[name][1]:.2f}\")")
    ax[0].plot([FLOOR_P, 1], [FLOOR_P, 1], "k--", lw=1)
    for t in (0.02,):
        ax[0].axhline(t, color="crimson", lw=0.8, ls=":")
        ax[0].axvline(t, color="crimson", lw=0.8, ls=":")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$p_{\rm chance}$, floor = $\sigma_{\rm sys,abs}$ only"
                     " (3.05\"/4.42\")")
    ax[0].set_ylabel(r"$p_{\rm chance}$, alternative floor")
    ax[0].legend(fontsize=8.5, loc="upper left")
    ax[0].set_title("Changing the pointing floor moves nothing", fontsize=10)

    for name, c in zip(keys, ("#333333",) + palette):
        b = np.sort(best_rows[name][f"p_{name}"].values)
        ax[1].step(np.concatenate(([1e-4], b)), np.arange(len(b) + 1) / len(b),
                   where="post", lw=2, color=c,
                   label=f"{name}: {int((b < 0.02).sum())} at p<0.02")
    ax[1].axvline(0.02, color="crimson", lw=0.8, ls=":")
    ax[1].set_xscale("log")
    ax[1].set_xlim(1e-4, 1)
    ax[1].set_xlabel(r"$p_{\rm chance}$ of the best match")
    ax[1].set_ylabel("cumulative fraction of the 63 sources")
    ax[1].legend(fontsize=8.5, loc="upper left")
    ax[1].set_title("Same 63 sources, four pointing floors", fontsize=10)

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    fig.savefig(OUT_PNG, dpi=170)
    print("saved", OUT_PNG)


if __name__ == "__main__":
    main()
