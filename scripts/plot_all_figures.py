#!/usr/bin/env python
"""Every figure in one command: `spt figures 60` draws the whole set at fov 60".

Two kinds of figure live in this repo and this script runs both:

  per-source panels   one figure per source, cut out of the FITS cache at a
                      field of view you choose. Each batch lands in its own
                      folder, named for that fov (source_panels_fov60/, ...),
                      so a rerun overwrites only figures drawn at the same fov.

  field figures       one figure for the whole sample — p_chance histograms,
                      separation diagnostics, the sky overlays. They have no
                      fov, so they are simply redrawn in place; skip them with
                      --panels-only when you are only re-cutting panels.

A step that fails does not stop the run: every failure is listed at the end and
the exit status is non-zero, so a batch of 73 is never lost to one bad panel.

Usage
-----
  spt figures 60                 every figure, panels at fov 60"
  spt figures 60 --panels-only   just the per-source panels
  spt figures                    same, at the default fov 45"
"""
import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

# (label, script, args-that-do-not-depend-on-fov); '--fov N' is appended to the
# per-source steps only.
PER_SOURCE = [
    ('source panels (all 73)', 'plot_source_panels.py', ['--all', '--quiet-warnings']),
    ('CSC2 candidates',        'plot_csc2_cutout_unwise.py', []),
]
FIELD = [
    ('Gaia p_chance histogram', 'plot_gaia_pchance_hist.py', []),
    ('Gaia verdict summary',    'plot_gaia_verdict.py', []),
    ('SPT vs WISE flux',        'plot_spt_vs_wise_flux.py', []),
    ('separation vs flux',      'plot_sep_vs_flux.py', []),
    ('separation vs SNR',       'plot_sep_vs_snr.py', []),
    ('SPT source overlay',      'plot_spt_sources_overlay.py', []),
    ('SPT blazar overlay',      'plot_spt_blazar_overlay.py', []),
]


def run(label, script, args):
    cmd = [sys.executable, os.path.join(HERE, script)] + args
    print(f'\n=== {label} — {script} {" ".join(args)}', flush=True)
    t0 = time.time()
    rc = subprocess.call(cmd)
    print(f'--- {label}: {"ok" if rc == 0 else f"FAILED (exit {rc})"} '
          f'in {time.time() - t0:.0f}s', flush=True)
    return rc


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('fov', nargs='?', type=float, default=45.0,
                    help='arcsec for the per-source panels (default 45)')
    ap.add_argument('--panels-only', action='store_true',
                    help='skip the field figures, which have no fov')
    ap.add_argument('--dpi', type=int, help='passed to the per-source panels')
    args = ap.parse_args()

    steps = [(lab, s, a + ['--fov', f'{args.fov:g}']
              + (['--dpi', str(args.dpi)] if args.dpi and s == 'plot_source_panels.py' else []))
             for lab, s, a in PER_SOURCE]
    if not args.panels_only:
        steps += FIELD

    print(f'{len(steps)} step(s); panels at fov {args.fov:g}"')
    t0, failed = time.time(), []
    for lab, script, a in steps:
        if run(lab, script, a) != 0:
            failed.append(lab)

    print(f'\n{"=" * 60}\nall figures: {len(steps) - len(failed)}/{len(steps)} steps ok '
          f'in {time.time() - t0:.0f}s')
    print(f'per-source panels -> outputs/images/*_fov{args.fov:g}/')
    if failed:
        print('FAILED: ' + ', '.join(failed))
        raise SystemExit(1)


if __name__ == '__main__':
    main()
