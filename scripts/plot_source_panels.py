#!/usr/bin/env python
"""Run src.source_figure.plot_source over one source, a few, or the whole list.

This script is only the loop: read the table, call plot_source once per row,
report what failed. All figure logic lives in src/source_figure.py, so a
change to how a panel looks is a change in one place and every caller —
this script, a notebook, an ad-hoc one-source check — gets it.

Usage
-----
  python scripts/plot_source_panels.py --all
  python scripts/plot_source_panels.py SPT3G_J173508.4-293000.8
  python scripts/plot_source_panels.py --all --skip-existing
  python scripts/plot_source_panels.py --all --panels W1,W2,W3,W4 --big DECaPS-r,NVSS

Panels are read from data/fits_cache/, so run scripts/precache_fits.py first;
after that a full run needs no network. Failures are listed on stdout and
written to outputs/logs/plot_source_panels_errors.csv.
"""
import argparse
import os
import sys
import time
import warnings

import matplotlib
matplotlib.use('Agg')
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from source_figure import DEFAULT_BIG, DEFAULT_SMALL, PANELS, plot_source

TABLE = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')
OUT_DIR = os.path.join(REPO, 'outputs', 'images', 'source_panels')
ERR_LOG = os.path.join(REPO, 'outputs', 'logs', 'plot_source_panels_errors.csv')


def _panel_list(arg, default):
    if not arg:
        return tuple(default)
    keys = tuple(k.strip() for k in arg.split(',') if k.strip())
    unknown = [k for k in keys if k not in PANELS]
    if unknown:
        raise SystemExit(f'unknown panel(s) {unknown}; '
                         f'available: {sorted(PANELS)}')
    return keys


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('ids', nargs='*', help='SPT ids; omit with --all')
    ap.add_argument('--all', action='store_true', help='every source in the table')
    ap.add_argument('--table', default=TABLE)
    ap.add_argument('--out-dir', default=OUT_DIR)
    ap.add_argument('--panels', help=f'small 2x2 block, comma-separated '
                                     f'(default {",".join(DEFAULT_SMALL)})')
    ap.add_argument('--big', help=f'full-height right column '
                                  f'(default {",".join(DEFAULT_BIG)})')
    ap.add_argument('--fov', type=float, default=45.0, help='arcsec (default 45)')
    ap.add_argument('--dpi', type=int, default=150)
    ap.add_argument('--skip-existing', action='store_true',
                    help='leave figures that are already on disk alone')
    ap.add_argument('--limit', type=int, help='stop after N sources (a dry run)')
    ap.add_argument('--quiet-warnings', action='store_true',
                    help='mute astropy WCS/FITS warnings (they are diagnostic '
                         '— only use this once a run is known clean)')
    args = ap.parse_args()

    if args.quiet_warnings:
        warnings.filterwarnings('ignore')
    small = _panel_list(args.panels, DEFAULT_SMALL)
    big = _panel_list(args.big, DEFAULT_BIG)

    df = pd.read_csv(args.table)
    if args.all:
        rows = df
    elif args.ids:
        rows = df[df['id'].isin(args.ids)]
        missing = set(args.ids) - set(rows['id'])
        if missing:
            print(f'!! not in {os.path.basename(args.table)}: {sorted(missing)}')
        if rows.empty:
            raise SystemExit('nothing to plot')
    else:
        raise SystemExit('give source ids or --all (see --help)')
    if args.limit:
        rows = rows.head(args.limit)

    os.makedirs(args.out_dir, exist_ok=True)
    print(f'{len(rows)} source(s) -> {args.out_dir}')
    print(f'panels: {list(small)} + {list(big)}, fov {args.fov:.0f}"\n')

    problems, t0 = [], time.time()
    for n, (_, row) in enumerate(rows.iterrows(), 1):
        sid = row['id']
        out_path = os.path.join(args.out_dir, f'{sid}.png')
        if args.skip_existing and os.path.isfile(out_path):
            print(f'[{n}/{len(rows)}] {sid}  (exists, skipped)')
            continue
        print(f'[{n}/{len(rows)}] {sid}', flush=True)
        try:
            _, failures = plot_source(row, out_path=out_path, small=small,
                                      big=big, fov_arcsec=args.fov,
                                      dpi=args.dpi)
        except Exception as e:
            print(f'    FAILED: {type(e).__name__}: {e}')
            problems.append(dict(spt_id=sid, panel='(whole figure)',
                                 error=f'{type(e).__name__}: {e}'))
            continue
        for key, err in failures.items():
            print(f'    [{key}] {err[:90]}')
            problems.append(dict(spt_id=sid, panel=key, error=err))

    dt = time.time() - t0
    print(f'\ndone in {dt:.0f}s ({dt / max(len(rows), 1):.1f}s per source)')
    if problems:
        os.makedirs(os.path.dirname(ERR_LOG), exist_ok=True)
        pd.DataFrame(problems).to_csv(ERR_LOG, index=False)
        bad = {p['spt_id'] for p in problems}
        print(f'{len(problems)} panel failure(s) across {len(bad)} source(s) '
              f'-> {ERR_LOG}')
    else:
        print('every panel drawn, no failures')


if __name__ == '__main__':
    main()
