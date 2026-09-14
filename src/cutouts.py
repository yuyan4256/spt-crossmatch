"""SPT cutout I/O + reprojection + aperture photometry."""
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import map_coordinates

from paths import SPT_CUTOUT_DIR, COADD_DIR


# ─────────── per-source cutouts (data/v3/) ───────────
def load_spt_cutout(source_id, band, cutout_dir=None):
    """Read individual SPT-3G cutout. Returns (data, wcs) or (None, None)."""
    cutout_dir = cutout_dir or SPT_CUTOUT_DIR
    path = os.path.join(cutout_dir, f'{source_id}_{band}.fits')
    if not os.path.isfile(path):
        return None, None
    with fits.open(path) as hdul:
        data = hdul['SCI_1'].data.copy()
        wcs = WCS(hdul['SCI_1'].header)
    return data, wcs


def reproject_to(src_data, src_wcs, tgt_wcs, tgt_shape, order=1):
    """Resample src_data onto target WCS grid (for contour overlays)."""
    yy, xx = np.mgrid[:tgt_shape[0], :tgt_shape[1]]
    world = tgt_wcs.pixel_to_world(xx, yy)
    src_xx, src_yy = src_wcs.world_to_pixel(world)
    return map_coordinates(src_data, [src_yy, src_xx],
                           order=order, mode='constant', cval=np.nan)


# ─────────── big coadd cutout (data/coadd/) ───────────
def open_coadd(band, coadd_dir=None):
    """Open galaxy_3yr_pc_gc_v2_{band}.fits and return its primary HDU data.

    Returns a numpy copy so the FITS file handle is closed before return.
    """
    coadd_dir = coadd_dir or COADD_DIR
    path = os.path.join(coadd_dir, f'galaxy_3yr_pc_gc_v2_{band}.fits')
    with fits.open(path) as hdul:
        return hdul[0].data.copy()


def coadd_cutout(arr, x, y, half=25):
    """Square cutout of 2*half px centered on (x, y); pads with NaN at edge."""
    ny, nx = arr.shape
    x, y = int(round(x)), int(round(y))
    out = np.full((2 * half, 2 * half), np.nan)
    x0, x1 = max(0, x - half), min(nx, x + half)
    y0, y1 = max(0, y - half), min(ny, y + half)
    out[y0 - (y - half):y1 - (y - half),
        x0 - (x - half):x1 - (x - half)] = arr[y0:y1, x0:x1]
    return out


# ─────────── aperture photometry ───────────
def aperture_phot(data, wcs, ra, dec, aper_arcsec=5.0, pixscale=2.75):
    """Aperture photometry at (ra, dec). Returns dict with bg-subtracted flux."""
    coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg)
    px, py = wcs.world_to_pixel(coord)
    px, py = float(px), float(py)
    ny, nx = data.shape
    yy, xx = np.mgrid[:ny, :nx]
    r_px = aper_arcsec / pixscale

    aper = ((xx - px) ** 2 + (yy - py) ** 2) <= r_px ** 2
    annulus = (((xx - px) ** 2 + (yy - py) ** 2) > (3 * r_px) ** 2) & \
              (((xx - px) ** 2 + (yy - py) ** 2) <= (5 * r_px) ** 2)

    if aper.sum() == 0:
        return None

    aper_flux = np.nansum(data[aper])
    n_pix = int(aper.sum())
    if annulus.sum() > 0:
        bg_per_pix = np.nanmedian(data[annulus])
        bg_std = np.nanstd(data[annulus])
        flux_sub = aper_flux - bg_per_pix * n_pix
        noise = bg_std * np.sqrt(n_pix)
        snr = flux_sub / noise if noise > 0 else 0.0
    else:
        flux_sub = aper_flux
        snr = 0.0

    return {'flux': flux_sub, 'snr': snr, 'n_pix': n_pix}
