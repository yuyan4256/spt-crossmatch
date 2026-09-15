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
from reproject.mosaicking import reproject_and_coadd

UNWISE_CACHE = os.path.join(os.path.dirname(__file__), '..', 'data', 'unwise')
UNWISE_PIXSCALE = 2.75    # arcsec/px
UNWISE_MAX_PX = 1024      # cutout_fits caps size here (checked 2026-09-14)


def _icrs_wcs(header):
    """WCS of an unWISE tile, in ICRS.

    unWISE headers carry `EPOCH = -1` ("epoch number" of the coadd), which
    astropy reads as an FK4 equinox of -1: any frame-aware step (SkyCoord
    markers, reprojection) then lands in the wrong place or returns NaN. The
    astrometry is Gaia-tied, i.e. ICRS.
    """
    hdr = header.copy()
    for key in ('EPOCH', 'EQUINOX'):
        hdr.remove(key, ignore_missing=True)
    hdr['RADESYS'] = 'ICRS'
    return WCS(hdr)


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

    fits_files = glob.glob(os.path.join(tile_dir, '**', f'*-w{band}-img-m.fits'),
                           recursive=True)
    if not fits_files:                       # not downloaded yet
        os.makedirs(tile_dir, exist_ok=True)
        url = (f'https://unwise.me/cutout_fits?version={version}'
               f'&ra={ra}&dec={dec}&size={size_px}&bands={band}')
        tar_path = os.path.join(tile_dir, 'cutout.tar.gz')
        print(f'  Downloading unWISE: {url[:80]}...')
        urllib.request.urlretrieve(url, tar_path)
        with tarfile.open(tar_path, 'r:gz') as tf:
            tf.extractall(tile_dir)
        os.remove(tar_path)
        # This band's image tiles only.
        fits_files = glob.glob(os.path.join(tile_dir, '**', f'*-w{band}-img-m.fits'),
                               recursive=True)
    if not fits_files:
        raise FileNotFoundError(f'No w{band} image tiles found in {tile_dir}')

    tiles = []
    for fp in sorted(fits_files):
        with fits.open(fp) as hdul:
            if hdul[0].data is not None and hdul[0].data.ndim == 2:
                tiles.append((hdul[0].data.astype(float), _icrs_wcs(hdul[0].header)))
    if not tiles:
        raise ValueError(f'No valid 2D FITS data in {tile_dir}')

    # The cutout headers carry a correct TAN WCS (tile tangent point, CRPIX in
    # tile coordinates). Use it. Until 2026-09-14 this function replaced it by
    # a WCS that put (ra, dec) at the array centre; unwise.me actually places
    # the requested position ~0.85 px off-centre in x and y, so every unWISE
    # panel was shifted by ~3″ (measured by cross-correlation with DECaPS).
    if len(tiles) == 1:
        data, wcs = tiles[0]                 # native pixels, no resampling
    else:
        # Position straddles unWISE tiles: coadd onto a north-up TAN grid at
        # the native 2.75″ scale whose centre pixel IS (ra, dec). No
        # background matching — the tiles are already sky-subtracted and
        # calibrated, and matching would shift their flux zero points.
        n = int(size_px)
        wcs = WCS(naxis=2)
        wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
        wcs.wcs.crval = [ra, dec]
        wcs.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
        wcs.wcs.cd = [[-2.75 / 3600.0, 0.0], [0.0, 2.75 / 3600.0]]
        data, _ = reproject_and_coadd(tiles, wcs, shape_out=(n, n),
                                      reproject_function=reproject_interp,
                                      match_background=False)
        data[_ == 0] = np.nan

    fits.writeto(mosaic_path, np.asarray(data, dtype=np.float32), wcs.to_header(),
                 overwrite=True)
    print(f'  Cached mosaic: {mosaic_path}')
    return data, wcs
