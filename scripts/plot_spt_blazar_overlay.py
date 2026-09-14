"""
Overplot the Milliquas blazar catalog on the SPT-3G 90 GHz source map, and
mark how many blazars coincide with the 61 variable sources.

Result (verified on the weight>0 footprint): the field holds only 6 Milliquas
AGN of ANY type and exactly ONE blazar (TXS 1714-336), which sits 2.1" from
SPT3G_J171736.0-334209.6. Roma-BZCAT (3 field blazars) and SIMBAD (2 AGN among
the 61) independently agree: only TXS 1714-336 coincides with a source. The
field is on the Galactic plane (l~0), a blazar zone of avoidance -- so the
single coincidence is not a chance alignment.

Layers:
  - 90 GHz background (Yujie's stretch, +-10000 uK), footprint masked white
  - 61 SPT variable sources  (red circles, as in plot_spt_sources_overlay.py)
  - Milliquas BLAZARS in footprint (Type contains 'B')  -> cyan star, labelled
  - other Milliquas AGN in footprint (context)          -> orange diamond
  - SIMBAD QSO / BL Lac among the 61 sources            -> magenta square
    (surfaces the ICRF QSO SPT3G_J171310.1-341827.9 that Milliquas omits)
  - Roma-BZCAT 5th-ed blazars in footprint              -> black ring
    (reads outputs/roma_bzcat_footprint.csv, written by check_roma_bzcat.py)
Builds on plot_spt_sources_overlay.py.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.lines import Line2D
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings; warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))
from units import spt_to_uk  # noqa: E402

BAND, VMIN, VMAX, RADIUS = '90', -10000, 10000, 20
FITS_PATH = os.path.join(REPO, 'data', 'coadd', f'galaxy_3yr_pc_gc_v2_{BAND}.fits')
CSV_PATH  = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')
MQ_PATH   = os.path.join(REPO, 'data', 'milliquas_v8.csv')
ROMA_PATH = os.path.join(REPO, 'outputs', 'roma_bzcat_footprint.csv')
OUT_PATH  = os.path.join(REPO, 'outputs', 'images', f'spt{BAND}_blazar_overlay.png')
MATCH_R   = 60.0   # arcsec: "coincident with a source" threshold

def is_blazar(t):
    return isinstance(t, str) and 'B' in t.strip().upper()

# ── map + weight ──
with fits.open(FITS_PATH) as hd:
    T_mK   = np.asarray(hd[0].data, float)
    weight = np.asarray(hd[1].data, float)
    wcs    = WCS(hd[0].header).celestial
ny, nx = weight.shape
map_uk = spt_to_uk(T_mK)
map_uk[weight <= 0] = np.nan

# ── 61 sources ──
src = pd.read_csv(CSV_PATH)
csrc = SkyCoord(src.ra_deg.values*u.deg, src.dec_deg.values*u.deg)

# ── Milliquas in footprint ──
ys, xs = np.where(weight > 0)
samp = slice(None, None, max(1, len(xs)//20000))
sky  = wcs.pixel_to_world(xs[samp], ys[samp])
ramin, ramax = sky.ra.deg.min(), sky.ra.deg.max()
demin, demax = sky.dec.deg.min(), sky.dec.deg.max()

mq = pd.read_csv(MQ_PATH, usecols=['RAJ2000','DEJ2000','Name','Type','Rmag','z'])
m = ((mq.RAJ2000.between(ramin-0.3, ramax+0.3)) & (mq.DEJ2000.between(demin-0.3, demax+0.3)))
box = mq[m].copy()
px, py = wcs.world_to_pixel(SkyCoord(box.RAJ2000.values*u.deg, box.DEJ2000.values*u.deg))
pxi, pyi = np.round(px).astype(int), np.round(py).astype(int)
ins = (pxi>=0)&(pxi<nx)&(pyi>=0)&(pyi<ny)
ins[ins] &= weight[pyi[ins], pxi[ins]] > 0
box['px'], box['py'], box['ins'] = px, py, ins
box['blz'] = box.Type.apply(is_blazar)
ftp   = box[box.ins].copy()
blz   = ftp[ftp.blz]
other = ftp[~ftp.blz]

# overlap count
cblz = SkyCoord(blz.RAJ2000.values*u.deg, blz.DEJ2000.values*u.deg)
_, sep, _ = cblz.match_to_catalog_sky(csrc)
n_over = int((sep.arcsec <= MATCH_R).sum())
print(f'footprint: {len(ftp)} Milliquas AGN, {len(blz)} blazar(s); '
      f'{n_over}/{len(blz)} blazar(s) within {MATCH_R:.0f}" of a source')

# ── SIMBAD QSO / BL Lac among the 61 sources (independent classification) ──
# SIMBAD otype: long label from the old script API ('BLLac'), short code from
# the TAP client used by astroquery>=0.4.8 ('BLL'). Accept both.
SIM_AGN = {'QSO', 'BLLac', 'BLL'}
sim = src[src.simbad_otype.isin(SIM_AGN)].copy()
sim_in_mq = sim.Milliquas_match.fillna(False).astype(bool)
print(f'SIMBAD QSO/BL Lac among the 61 sources: {len(sim)} '
      f'({int(sim_in_mq.sum())} also in Milliquas, {int((~sim_in_mq).sum())} Milliquas-missed)')

# ── Roma-BZCAT blazars in footprint (from check_roma_bzcat.py; no network here) ──
ROMA_ABBR = {'QSO RLoud flat radio sp.': 'FSRQ', 'Blazar Uncertain type': 'uncertain',
             'BL Lac Candidate': 'BL Lac', 'BL Lac': 'BL Lac', 'BL Lac-galaxy dominated': 'BL Lac (gal)'}
if os.path.exists(ROMA_PATH):
    roma = pd.read_csv(ROMA_PATH)
    rpx, rpy = wcs.world_to_pixel(SkyCoord(roma.ra_deg.values*u.deg, roma.dec_deg.values*u.deg))
    roma['px'], roma['py'] = rpx, rpy
    print(f'Roma-BZCAT in footprint: {len(roma)} blazar(s); '
          f'{int(roma.coincident.sum())} coincide with a source')
else:
    roma = None
    print(f'(!) {ROMA_PATH} missing — run scripts/check_roma_bzcat.py to add the Roma layer')

# ── figure ──
cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('white')
fig = plt.figure(figsize=(13, 8))
ax  = fig.add_subplot(111, projection=wcs)
im  = ax.imshow(map_uk, origin='lower', vmin=VMIN, vmax=VMAX, cmap=cmap)

for r in src.itertuples():
    ax.add_patch(Circle((r.xcentroid, r.ycentroid), radius=RADIUS,
                        edgecolor='red', facecolor='none', lw=1.0, alpha=0.7))
ax.scatter(other.px, other.py, marker='D', s=60, facecolor='orange',
           edgecolor='black', lw=0.6, zorder=5)
ax.scatter(blz.px, blz.py, marker='*', s=420, facecolor='cyan',
           edgecolor='black', lw=1.0, zorder=6)
ax.scatter(sim.xcentroid, sim.ycentroid, marker='s', s=620, facecolor='none',
           edgecolor='magenta', lw=2.0, zorder=7)
if roma is not None:
    ax.scatter(roma.px, roma.py, marker='o', s=900, facecolor='none',
               edgecolor='black', lw=2.2, zorder=4)
    for r in roma.itertuples():        # name the Roma field blazars that miss every source
        if not r.coincident:
            right = r.px > 0.62 * nx
            ax.annotate(f'{r.Name} ({ROMA_ABBR.get(r.Class, r.Class)})',
                        xy=(r.px, r.py), xytext=(r.px + (-110 if right else 110), r.py + 150),
                        fontsize=7.5, color='black', ha='right' if right else 'left', va='center')

# label the two SIMBAD AGN sources (one is the Milliquas-missed ICRF QSO)
offs = {'BLLac': (-260, 360), 'BLL': (-260, 360), 'QSO': (170, -430)}
for r in sim.itertuples():
    in_mq = bool(r.Milliquas_match) if not pd.isna(r.Milliquas_match) else False
    tag = f'in Milliquas ({r.Milliquas_Name})' if in_mq else 'NOT in Milliquas'
    ec  = '0.35' if in_mq else 'magenta'
    otype = 'BL Lac' if r.simbad_otype in ('BLLac', 'BLL') else r.simbad_otype
    dx, dy = offs.get(r.simbad_otype, (200, 200))
    ax.annotate(f'{r.simbad_name}\n{otype} · {tag}',
                xy=(r.xcentroid, r.ycentroid), xytext=(r.xcentroid+dx, r.ycentroid+dy),
                fontsize=8.5, ha='left', va='center',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=ec,
                          lw=1.6 if not in_mq else 0.8, alpha=0.93),
                arrowprops=dict(arrowstyle='->', color=ec, lw=1.1))

handles = [
    Line2D([],[], marker='o', ls='none', mfc='none', mec='red', ms=10, label=f'SPT variable sources ({len(src)})'),
    Line2D([],[], marker='*', ls='none', mfc='cyan', mec='black', ms=16, label=f'Milliquas blazar in field ({len(blz)})'),
    Line2D([],[], marker='s', ls='none', mfc='none', mec='magenta', mew=2.0, ms=13, label=f'SIMBAD QSO / BL Lac source ({len(sim)})'),
    Line2D([],[], marker='o', ls='none', mfc='none', mec='black', mew=2.0, ms=14, label=f'Roma-BZCAT blazar in field ({0 if roma is None else len(roma)})'),
    Line2D([],[], marker='D', ls='none', mfc='orange', mec='black', ms=8, label=f'Other Milliquas AGN in field ({len(other)})'),
]
ax.legend(handles=handles, loc='lower left', fontsize=9, framealpha=0.92).set_zorder(10)

lon, lat = ax.coords[0], ax.coords[1]
lon.set_axislabel('RA'); lat.set_axislabel('DEC')
lon.set_major_formatter('hh:mm'); lat.set_major_formatter('dd')
ax.set_title(f'Blazar/AGN cross-check over the SPT-3G {BAND} GHz field: only TXS 1714-336 '
             f'coincides with a source (Milliquas, Roma-BZCAT & SIMBAD agree)', fontsize=10)
cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, extend='both')
cbar.set_label(r'T ($\mu$K$_{CMB}$)')

fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
fig.savefig(OUT_PATH, dpi=140, facecolor='white', bbox_inches='tight')
print(f'saved -> {OUT_PATH}')
