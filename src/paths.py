"""Project-relative paths so scripts/notebooks work from any CWD."""
import os
import re

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

DATA = os.path.join(ROOT, 'data')
OUT = os.path.join(ROOT, 'outputs')

SPT_CUTOUT_DIR = os.path.join(DATA, 'v3')        # v3 (32 dual-band sources) cutouts
V4_CUTOUT_DIR = os.path.join(DATA, 'v4')         # 90/150 GHz cutouts of the 3-yr/61 run (ids != the 73)
COADD_DIR = os.path.join(DATA, 'coadd')          # big galaxy_3yr_pc_gc_v2_{90,150}.fits
UNWISE_DIR = os.path.join(DATA, 'unwise')        # unWISE cache
IMAGES_DIR = os.path.join(OUT, 'images')
LEGACY_DIR = os.path.join(OUT, 'legacy')


def fov_dir(base, fov_arcsec):
    """outputs/images/<base>_fov<NN> — a batch of figures names its fov.

    Every batch goes in its own folder, so a rerun can only overwrite figures
    drawn at the same fov; a different fov lands next to them, not on top.
    """
    return os.path.join(IMAGES_DIR, f'{base}_fov{fov_arcsec:g}')


def check_fov_dir(out_dir, fov_arcsec):
    """Refuse an --out-dir whose name says a different fov than the one asked."""
    m = re.search(r'_fov(\d+(?:\.\d+)?)$', os.path.basename(os.path.normpath(out_dir)))
    if m and float(m.group(1)) != float(fov_arcsec):
        raise SystemExit(f'--out-dir {out_dir!r} is a fov {m.group(1)}" folder but '
                         f'--fov is {fov_arcsec:g}"; pick one (figures at different '
                         f'fov must not share a folder)')
    return out_dir
