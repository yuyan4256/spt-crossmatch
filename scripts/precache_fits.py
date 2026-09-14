#!/usr/bin/env python
"""Pre-download every figure panel for all 73 sources, once, at full size.

Each survey is fetched at fits_cutouts.DOWNLOAD_FOV_ARCSEC (786″) at its
native pixel scale into data/fits_cache/ (unWISE into data/unwise/); figures
crop to whatever fov they ask for. Resumable: fetchers skip files already in
the cache. RACS-mid needs CASDA_USERNAME set (see fits_cutouts.fetch_racs_fits).

  python scripts/precache_fits.py                 # everything
  python scripts/precache_fits.py --only W1,W2    # a subset of panels
"""
import argparse, os, sys, time, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))

import pandas as pd
from fits_cutouts import (DOWNLOAD_FOV_ARCSEC, fetch_allwise_atlas,
                          fetch_decaps_fits, fetch_nvss_fits, fetch_racs_fits)
from unwise import UNWISE_PIXSCALE, fetch_unwise
import numpy as np

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
v4 = pd.read_csv(os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv'))
print(f'{len(v4)} sources')

UNWISE_PX = int(np.ceil(DOWNLOAD_FOV_ARCSEC / UNWISE_PIXSCALE))

JOBS = [('W1',    lambda ra, dec, sid: fetch_unwise(ra, dec, source_id=sid, size_px=UNWISE_PX, band=1)),
        ('W2',    lambda ra, dec, sid: fetch_unwise(ra, dec, source_id=sid, size_px=UNWISE_PX, band=2)),
        ('W3',    lambda ra, dec, sid: fetch_allwise_atlas(ra, dec, 'W3')),
        ('W4',    lambda ra, dec, sid: fetch_allwise_atlas(ra, dec, 'W4')),
        ('NVSS',  lambda ra, dec, sid: fetch_nvss_fits(ra, dec)),
        ('DECaPS', lambda ra, dec, sid: fetch_decaps_fits(ra, dec)),
        ('RACS',   lambda ra, dec, sid: fetch_racs_fits(ra, dec))]

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--only', help='comma-separated subset of ' + ','.join(n for n, _ in JOBS))
args = ap.parse_args()
if args.only:
    keep = set(args.only.split(','))
    JOBS = [j for j in JOBS if j[0] in keep]
print(f'panels {[n for n, _ in JOBS]} at {DOWNLOAD_FOV_ARCSEC:.0f}″')

fails = []
for i, r in v4.iterrows():
    ra, dec = float(r['ra_deg']), float(r['dec_deg'])
    for name, fn in JOBS:
        ok = False
        for attempt in (1, 2, 3):
            try:
                fn(ra, dec, r['id'])
                ok = True
                break
            except Exception as e:
                if attempt < 3:
                    time.sleep(10 * attempt)
                else:
                    print(f'[{i+1}/{len(v4)}] {r["id"]} {name}: '
                          f'{type(e).__name__}: {str(e)[:60]}')
        if not ok:
            fails.append((r['id'], name))
    print(f'[{i+1}/{len(v4)}] {r["id"]} done')

print(f'\n{len(fails)} failures')
for sid, name in fails[:20]:
    print(f'  {sid} {name}')
