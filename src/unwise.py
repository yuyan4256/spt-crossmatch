"""unWISE W1/W2 cutouts (neo6, Vega nanomaggies, 2.75″/px) with a disk cache.

Each (band, size) request is extracted into its own folder and only that
band's `-w{band}-img-m` tiles are mosaicked. Until 2026-09-14 W1 and W2 tiles
were extracted into one folder and globbed together, so every cached "W2"
mosaic was the W1+W2 average; caches from that scheme (unwise_w{1,2}.fits
without a size tag) must not be used.
"""
import os, tarfile, glob
import urllib.request
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd

UNWISE_CACHE = os.path.join(os.path.dirname(__file__), '..', 'data', 'unwise')
UNWISE_PIXSCALE = 2.75    # arcsec/px
UNWISE_MAX_PX = 1024      # cutout_fits caps size here (checked 2026-09-14)


def fetch_unwise(ra, dec, source_id=None, size_px=100, band=1, version='neo6',
                 cache_dir=None):
    """
    Download unWISE W1 cutout and return (data_2d, wcs).

    Parameters
    ----------
    ra, dec : float
        J2000 coordinates in degrees.
    source_id : str, optional
        Used for caching. If None, uses "ra_dec".
    size_px : int
        Cutout size in unWISE pixels (2.75"/px). Default 100 ~ 4.6'.
    band : int
        1 = W1 (3.4 um), 2 = W2 (4.6 um).
    version : str
        unWISE version (default 'neo6').
    cache_dir : str, optional
        Override default cache directory.

    Returns
    -------
    data : 2D ndarray
    wcs : astropy.wcs.WCS
    """
    if cache_dir is None:
        cache_dir = os.path.abspath(UNWISE_CACHE)

    if source_id is None:
        source_id = f'{ra:.4f}_{dec:.4f}'

    if size_px > UNWISE_MAX_PX:
        raise ValueError(f'size_px={size_px}: unwise.me caps cutouts at {UNWISE_MAX_PX} px')
    src_dir = os.path.join(cache_dir, source_id)
    mosaic_path = os.path.join(src_dir, f'unwise_w{band}_{size_px}px.fits')
    tile_dir = os.path.join(src_dir, f'tiles_w{band}_{size_px}px')

    # Return cached mosaic if exists
    if os.path.isfile(mosaic_path):
        with fits.open(mosaic_path) as hdul:
            return hdul[0].data.copy(), WCS(hdul[0].header)

    # Download
    os.makedirs(tile_dir, exist_ok=True)
    url = (f'https://unwise.me/cutout_fits?version={version}'
           f'&ra={ra}&dec={dec}&size={size_px}&bands={band}')
    tar_path = os.path.join(tile_dir, 'cutout.tar.gz')

    print(f'  Downloading unWISE: {url[:80]}...')
    urllib.request.urlretrieve(url, tar_path)

    # Extract FITS files
    with tarfile.open(tar_path, 'r:gz') as tf:
        tf.extractall(tile_dir)
    os.remove(tar_path)

    # This band's image tiles only.
    fits_files = glob.glob(os.path.join(tile_dir, '**', f'*-w{band}-img-m.fits'),
                           recursive=True)
    if not fits_files:
        raise FileNotFoundError(f'No w{band} image tiles found in {tile_dir}')

    # Load HDUs
    hdus = []
    for fp in fits_files:
        hdu = fits.open(fp)[0]
        if hdu.data is not None and hdu.data.ndim == 2:
            hdus.append(hdu)

    if len(hdus) == 0:
        raise ValueError(f'No valid 2D FITS data in {src_dir}')

    if len(hdus) == 1:
        data = hdus[0].data.astype(float)
        wcs = WCS(hdus[0].header)
    else:
        # Mosaic multiple tiles
        wcs_out, shape_out = find_optimal_celestial_wcs(hdus)
        data, footprint = reproject_and_coadd(
            hdus, wcs_out, shape_out=shape_out,
            reproject_function=reproject_interp,
            match_background=True)
        wcs = wcs_out

    # Build a clean TAN WCS centered on the request coordinates.
    # The original tile WCS has CRPIX far outside the cutout, which
    # confuses Cutout2D.  We know the unWISE pixel scale (2.75"/px)
    # and that the cutout is centered on (ra, dec).
    from astropy.wcs import WCS as _WCS
    ny, nx = data.shape
    clean_wcs = _WCS(naxis=2)
    clean_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    clean_wcs.wcs.crval = [ra, dec]
    clean_wcs.wcs.crpix = [nx / 2.0 + 0.5, ny / 2.0 + 0.5]
    clean_wcs.wcs.cdelt = [-2.75 / 3600.0, 2.75 / 3600.0]   # 2.75"/px
    wcs = clean_wcs

    # Save mosaic for caching
    header = wcs.to_header()
    fits.writeto(mosaic_path, data.astype(np.float32), header, overwrite=True)
    print(f'  Cached mosaic: {mosaic_path}')

    return data, wcs
