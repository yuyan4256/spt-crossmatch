"""Project-relative paths so scripts/notebooks work from any CWD."""
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(ROOT, 'outputs')

SPT_CUTOUT_DIR = os.path.join(DATA, 'v3')        # v3 (32 dual-band sources) cutouts
V4_CUTOUT_DIR = os.path.join(DATA, 'v4')         # 90/150 GHz cutouts of the 3-yr/61 run (ids != the 73)
COADD_DIR = os.path.join(DATA, 'coadd')          # big galaxy_3yr_pc_gc_v2_{90,150}.fits
UNWISE_DIR = os.path.join(DATA, 'unwise')        # unWISE cache
IMAGES_DIR = os.path.join(OUT, 'images')
LEGACY_DIR = os.path.join(OUT, 'legacy')
