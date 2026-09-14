#!/usr/bin/env python
"""Pre-download all figure FITS into data/fits_cache/ (resumable by design:
fetchers skip existing cache files). Run before figure notebooks so plot
runs are fast and offline-capable."""
import os, sys, time, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))

import pandas as pd
from fits_cutouts import (fetch_allwise_atlas, fetch_decaps_fits,
                          fetch_nvss_fits, fetch_racs_fits)

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
v4 = pd.read_csv(os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv'))
print(f'{len(v4)} sources')

JOBS = [('W3',    lambda ra, dec: fetch_allwise_atlas(ra, dec, 'W3')),
        ('W4',    lambda ra, dec: fetch_allwise_atlas(ra, dec, 'W4')),
        ('NVSS',  lambda ra, dec: fetch_nvss_fits(ra, dec)),
        ('DECaPS', lambda ra, dec: fetch_decaps_fits(ra, dec)),
        ('RACS',   lambda ra, dec: fetch_racs_fits(ra, dec))]

fails = []
for i, r in v4.iterrows():
    ra, dec = float(r['ra_deg']), float(r['dec_deg'])
    for name, fn in JOBS:
        ok = False
        for attempt in (1, 2, 3):
            try:
                fn(ra, dec)
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
