"""
Independent blazar cross-check with Roma-BZCAT 5th ed (VizieR VII/274) — a
dedicated all-sky blazar catalog — against Milliquas (which flagged exactly 1
blazar in the SPT-3G field).

Finds every Roma-BZCAT blazar inside the 90 GHz weight>0 footprint and its
nearest of the 61 SPT variable sources. Saves outputs/roma_bzcat_footprint.csv.

Result: Roma-BZCAT has 3 blazars in the field (vs Milliquas's 1), but only ONE
coincides with a source — 5BZB J1717-3342 = TXS 1714-336 (0.9"), the same object
Milliquas found. So the single blazar counterpart is confirmed independently.
The other 2 Roma blazars sit ~1 deg from any source. Neither Roma nor Milliquas
contains the SIMBAD QSO SPT3G_J171310.1-341827.9 (Galactic-plane radio AGN).
"""
import os
import numpy as np, pandas as pd
from astroquery.vizier import Vizier
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings; warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FITS = os.path.join(REPO, 'data', 'coadd', 'galaxy_3yr_pc_gc_v2_90.fits')
SRC  = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')
OUT  = os.path.join(REPO, 'outputs', 'roma_bzcat_footprint.csv')
MATCH_R = 60.0  # arcsec

# ── footprint ──
with fits.open(FITS) as hd:
    weight = np.asarray(hd[1].data, float); wcs = WCS(hd[0].header).celestial
ny, nx = weight.shape
ys, xs = np.where(weight > 0)
samp = slice(None, None, max(1, len(xs)//20000))
sky = wcs.pixel_to_world(xs[samp], ys[samp])
ramin, ramax = sky.ra.deg.min(), sky.ra.deg.max()
demin, demax = sky.dec.deg.min(), sky.dec.deg.max()

# ── Roma-BZCAT VII/274 (all-sky) ──
v = Vizier(columns=['**'], row_limit=-1); v.TIMEOUT = 120
tl = v.get_catalogs('VII/274')
df = tl[next(n for n in tl.keys() if 'RAJ2000' in tl[n].colnames)].to_pandas()
print(f'Roma-BZCAT: {len(df)} blazars all-sky | {dict(df.Class.value_counts())}')
coo = SkyCoord(df.RAJ2000.astype(str).values, df.DEJ2000.astype(str).values,
               unit=(u.hourangle, u.deg))
df['ra_deg'], df['dec_deg'] = coo.ra.deg, coo.dec.deg

# ── inside footprint ──
m = df.ra_deg.between(ramin-0.3, ramax+0.3) & df.dec_deg.between(demin-0.3, demax+0.3)
box = df[m].copy()
px, py = wcs.world_to_pixel(SkyCoord(box.ra_deg.values*u.deg, box.dec_deg.values*u.deg))
pxi, pyi = np.round(px).astype(int), np.round(py).astype(int)
ins = (pxi>=0)&(pxi<nx)&(pyi>=0)&(pyi<ny)
ins[ins] &= weight[pyi[ins], pxi[ins]] > 0
ftp = box[ins].copy()

# ── nearest SPT source ──
src = pd.read_csv(SRC)
csrc = SkyCoord(src.ra_deg.values*u.deg, src.dec_deg.values*u.deg)
c = SkyCoord(ftp.ra_deg.values*u.deg, ftp.dec_deg.values*u.deg)
idx, sep, _ = c.match_to_catalog_sky(csrc)
ftp['nearest_source_id'] = src.id.values[idx]
ftp['sep_arcsec'] = np.round(sep.arcsec, 2)
ftp['coincident'] = ftp.sep_arcsec <= MATCH_R

n_over = int(ftp.coincident.sum())
print(f'in footprint: {len(ftp)} Roma blazars; {n_over} coincide with a source (<= {MATCH_R:.0f}")')
keep = ['Name', 'Class', 'ra_deg', 'dec_deg', 'z', 'Rmag',
        'nearest_source_id', 'sep_arcsec', 'coincident']
out = ftp[keep].sort_values('sep_arcsec')
print(out.to_string(index=False))
out.to_csv(OUT, index=False)
print(f'\nsaved -> {OUT}')
