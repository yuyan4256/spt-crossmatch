#!/usr/bin/env python
"""CSC2-candidate cutout figures: the standard per-source panels plus the
Chandra position and its 95% error circle.

This is `src.source_figure.plot_csc2_source` (= `plot_source` with two extra markers) — the figure
layout, the panel fetching, the SPT contour convention and the error handling
all live in that one module. Until 2026-09 this script carried its own copies
of `_load_spt`, `_reproject_spt`, `_off_source_rms` and `_draw_spt_contours`,
which is why a contour tweak used to have to be made in four places.
`csc2_marker` now lives in src/source_figure.py too (the notebook uses it).

CHANGED 2026-09: W3/W4, DECaPS r and RACS-mid now come from the calibrated
FITS cache instead of HiPS preview tiles, so every colorbar is quantitative
(2026-07 spec). Regenerating a figure therefore does not reproduce the old
one pixel-for-pixel — it replaces it with the version whose colorbars mean
something.

Usage:
  python scripts/plot_csc2_cutout_unwise.py                     # all CSC2 matches
  python scripts/plot_csc2_cutout_unwise.py SPT3G_J172051.9-354642.8
  python scripts/plot_csc2_cutout_unwise.py --out-dir /tmp/check
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from source_figure import csc2_marker, plot_csc2_source

TABLE = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')
OUT_DIR = os.path.join(REPO, 'outputs', 'images', 'v4_csc2_candidates')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ids', nargs='*', help='SPT ids; default = every CSC2 match')
    ap.add_argument('--out-dir', default=OUT_DIR)
    ap.add_argument('--fov', type=float, default=45.0)
    args = ap.parse_args()

    df = pd.read_csv(TABLE)
    if args.ids:
        rows = df[df['id'].isin(args.ids)]
        missing = set(args.ids) - set(rows['id'])
        if missing:
            print(f'!! not in the table: {sorted(missing)}')
    else:
        rows = df[df['CSC2_cat_ra'].notna()]
    if rows.empty:
        raise SystemExit('nothing to plot')

    os.makedirs(args.out_dir, exist_ok=True)
    print(f'{len(rows)} source(s) -> {args.out_dir}')
    for n, (_, row) in enumerate(rows.iterrows(), 1):
        sid = row['id']
        if csc2_marker(row) is None:
            print(f'[{n}/{len(rows)}] {sid}  (no CSC2 position — SPT marker only)')
        out = os.path.join(args.out_dir, f'{sid}.png')
        print(f'[{n}/{len(rows)}] {sid}', flush=True)
        _, failures = plot_csc2_source(row, out_path=out, fov_arcsec=args.fov)
        for key, err in failures.items():
            print(f'    [{key}] {err[:90]}')


if __name__ == '__main__':
    main()
