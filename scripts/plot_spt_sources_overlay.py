"""
Reproduce Yujie's seed figure locally (no sptgrid / spt3g_software needed).

Yujie's snippet (on sptgrid):
    map90 = frame['T'] / core.G3Units.uK          # after RemoveWeights
    ax.imshow(map90, origin='lower', vmin=-10000, vmax=10000)
    Circle((x[i], y[i]), radius=20, facecolor='none')

Local equivalent: our FITS HDU0 is the already-weight-removed T in mK_CMB,
so spt_to_uk() (x1000) == frame['T']/core.G3Units.uK. Verified: peak of
SPT3G_J170918.8-352514.0 = 3740.1 uK in both.

This builds ONLY the two fully-specified layers (90 GHz background + the 61
source circles), plus two layers I own and don't need to guess at:
  - RA/DEC (WCS) axes  - white survey-footprint exterior (HDU1 weight > 0)
NOT included (pending Yujie): black contours, classification text labels.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from astropy.io import fits
from astropy.wcs import WCS
import warnings
warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))
from units import spt_to_uk  # noqa: E402

BAND   = '90'                 # Yujie used 90GHz
VMIN, VMAX = -10000, 10000    # uK_CMB, Yujie's stretch
RADIUS = 20                   # pixels, Yujie's value
FITS_PATH = os.path.join(REPO, 'data', 'coadd', f'galaxy_3yr_pc_gc_v2_{BAND}.fits')
CSV_PATH  = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')
OUT_PATH  = os.path.join(REPO, 'outputs', 'images', f'spt{BAND}_sources_overlay.png')

# ── load map (HDU0 = weight-removed T in mK) + weight (HDU1) ──
with fits.open(FITS_PATH) as hd:
    T_mK   = np.asarray(hd[0].data, dtype=float)
    weight = np.asarray(hd[1].data, dtype=float)
    wcs    = WCS(hd[0].header).celestial

map_uk = spt_to_uk(T_mK)                    # mK -> uK_CMB  (== frame['T']/G3Units.uK)
map_uk[weight <= 0] = np.nan                # survey-footprint exterior -> white

# ── sources ──
df = pd.read_csv(CSV_PATH)
print(f'{len(df)} sources from {os.path.basename(CSV_PATH)}; map {map_uk.shape} ({BAND} GHz)')

# ── figure ──
cmap = plt.get_cmap('viridis').copy()
cmap.set_bad('white')

fig = plt.figure(figsize=(13, 8))
ax  = fig.add_subplot(111, projection=wcs)
im  = ax.imshow(map_uk, origin='lower', vmin=VMIN, vmax=VMAX, cmap=cmap)

for r in df.itertuples():
    ax.add_patch(Circle((r.xcentroid, r.ycentroid), radius=RADIUS,
                        edgecolor='red', facecolor='none', lw=1.2, alpha=0.85))

lon, lat = ax.coords[0], ax.coords[1]
lon.set_axislabel('RA');  lat.set_axislabel('DEC')
lon.set_major_formatter('hh:mm'); lat.set_major_formatter('dd')
ax.set_title(f'SPT-3G galaxy_3yr_pc_gc_v2  {BAND} GHz  +  {len(df)} variable sources',
             fontsize=11)

cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, extend='both')
cbar.set_label(r'T ($\mu$K$_{CMB}$)')

fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
fig.savefig(OUT_PATH, dpi=140, facecolor='white', bbox_inches='tight')
print(f'saved -> {OUT_PATH}')
