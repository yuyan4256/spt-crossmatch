"""Calibrated-FITS cutout fetchers with a local disk cache.

Replaces HiPS preview tiles for figure panels (2026-07 spec: plots must be
drawn from downloaded FITS). Every fetcher returns (data, wcs) with data in
the stated physical unit, and caches the raw download under data/fits_cache/
so reruns are offline-capable.

Surveys
-------
AllWISE W3/W4 : SkyView 'WISE 12'/'WISE 22' — AllWISE Atlas mosaics, native
                DN; converted to mJy (Vega zeropoints, same constants as the
                AllWISE catalog convention).
DECaPS r      : legacysurvey.org fits-cutout, layer=decaps2 — calibrated
                nanomaggies (AB); converted to μJy. Quantitative, unlike the
                DECaPS HiPS preview.
NVSS          : SkyView 'NVSS' — Jy/beam; converted to mJy/beam.
RACS-mid      : CASDA SODA cutout of the restored Stokes-I tile, 2″/px —
                Jy/beam; converted to mJy/beam. Needs an OPAL login
                (see fetch_racs_fits).

unWISE W1/W2 are handled by src/unwise.py (already FITS + cache).

Download once, crop at plot time
--------------------------------
Every fetcher downloads DOWNLOAD_FOV_ARCSEC (786″) at the survey's NATIVE
pixel scale, whatever the figure will show; `source_figure._crop` then trims
to the requested fov. 786″ is the largest DECaPS field the legacysurvey
cutout server returns in one request at 0.262″/px (3000 px cap, checked
2026-09-14). The cache filename carries the download size, so files from the
old fixed-45″ scheme are never mistaken for full-size ones.
"""
import os

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

_HERE      = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR  = os.path.join(_HERE, '..', 'data', 'fits_cache')

DOWNLOAD_FOV_ARCSEC = 786.0

# native pixel scales, arcsec/px
DECAPS_PIXSCALE = 0.262      # legacysurvey cutout server cap: 3000 px
ALLWISE_PIXSCALE = 1.375     # AllWISE Atlas via SkyView
NVSS_PIXSCALE = 15.0         # NVSS via SkyView
DECAPS_MAX_PX = 3000

# AllWISE Atlas DN -> mJy (Vega), same constants used for the catalog mags.
ALLWISE_DN_TO_MJY = {'W1': 1.935e-3, 'W2': 2.7048e-3,
                     'W3': 1.8326e-3, 'W4': 0.052269}
SKYVIEW_NAME = {'W1': 'WISE 3.4', 'W2': 'WISE 4.6',
                'W3': 'WISE 12',  'W4': 'WISE 22'}

# Legacy Survey / DECaPS nanomaggies are the AB kind: 1 nmgy = 3631 Jy ×
# 10^-9 = 3.631 μJy. Do NOT reuse this for unWISE — unWISE nanomaggies are
# VEGA (F0 = 309.540 / 171.787 Jy; see src/units.py).
NMGY_TO_UJY = 3.631


def _cache_path(tag, ra_deg, dec_deg, fov_arcsec):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR,
                        f'{tag}_{fov_arcsec:.0f}as_{ra_deg:.5f}_{dec_deg:.5f}.fits')


def _skyview_pixels(fov_arcsec, native_arcsec):
    """SkyView resamples to `pixels` (default 300) whatever the width, so ask
    for the native pixel count explicitly."""
    return int(np.ceil(fov_arcsec / native_arcsec))


def _north_up(data, wcs):
    """Resample a 2-D cutout onto a north-up / east-left TAN grid.

    Cutout servers deliver arbitrary pixel-grid orientations (the
    legacysurvey DECaPS cutouts e.g. come axis-swapped), so multi-survey
    panel figures end up with mixed orientations. Regrid everything onto a
    canonical WCS: same sky center, same mean pixel scale,
    CD = [[-s, 0], [0, s]] (RA increases leftward, Dec upward).
    Images already in canonical orientation are returned untouched.
    """
    from astropy.wcs.utils import proj_plane_pixel_scales
    from scipy.ndimage import map_coordinates

    ny, nx = data.shape
    s = float(np.mean(proj_plane_pixel_scales(wcs)))
    target_cd = np.array([[-s, 0.0], [0.0, s]])
    if np.allclose(wcs.pixel_scale_matrix, target_cd, atol=s * 1e-3):
        return data, wcs

    center = wcs.pixel_to_world((nx - 1) / 2, (ny - 1) / 2)
    n = max(nx, ny)
    tgt = WCS(naxis=2)
    tgt.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    tgt.wcs.crval = [center.ra.deg, center.dec.deg]
    tgt.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    tgt.wcs.cd = target_cd
    yy, xx = np.mgrid[:n, :n]
    sx, sy = wcs.world_to_pixel(tgt.pixel_to_world(xx, yy))
    out = map_coordinates(np.asarray(data, float), [sy, sx], order=1,
                          mode='constant', cval=np.nan)
    return out, tgt


def _read(path):
    with fits.open(path) as hdul:
        hdu = next(h for h in hdul if h.data is not None)
        data, wcs = np.squeeze(hdu.data.astype(float)), WCS(hdu.header).celestial
        return _north_up(data, wcs)


def _validate_or_remove(path):
    """Open-check a downloaded FITS; delete and raise if corrupt so a bad
    partial download can never poison the cache."""
    try:
        with fits.open(path) as hdul:
            next(h for h in hdul if h.data is not None)
    except Exception as e:
        os.remove(path)
        raise RuntimeError(f'corrupt download removed: {e}')


def fetch_allwise_atlas(ra_deg, dec_deg, band, fov_arcsec=DOWNLOAD_FOV_ARCSEC):
    """AllWISE Atlas cutout via SkyView at 1.375″/px. Returns (data_mJy_Vega, wcs)."""
    path = _cache_path(f'allwise_{band}', ra_deg, dec_deg, fov_arcsec)
    if not os.path.isfile(path):
        from astroquery.skyview import SkyView
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        imgs = SkyView.get_images(
            position=SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg),
            survey=[SKYVIEW_NAME[band]],
            width=fov_arcsec*u.arcsec, height=fov_arcsec*u.arcsec,
            pixels=_skyview_pixels(fov_arcsec, ALLWISE_PIXSCALE))
        imgs[0].writeto(path, overwrite=True)
        print(f'  Downloaded AllWISE {band}: {os.path.basename(path)}')
    data_dn, wcs = _read(path)
    return data_dn * ALLWISE_DN_TO_MJY[band], wcs


def fetch_decaps_fits(ra_deg, dec_deg, band='r', fov_arcsec=DOWNLOAD_FOV_ARCSEC,
                      pixscale=DECAPS_PIXSCALE):
    """DECaPS DR2 calibrated cutout via the Legacy Survey cutout server.

    Returns (data_uJy_AB, wcs). Pixel values are real fluxes (nanomaggies
    converted to μJy) — colorbar-safe, unlike the DECaPS HiPS preview.
    """
    path = _cache_path(f'decaps_{band}', ra_deg, dec_deg, fov_arcsec)
    if not os.path.isfile(path):
        size = int(round(fov_arcsec / pixscale))
        if size > DECAPS_MAX_PX:
            raise ValueError(f'{fov_arcsec}″ at {pixscale}″/px is {size} px; the '
                             f'cutout server silently caps at {DECAPS_MAX_PX} px')
        url = ('https://www.legacysurvey.org/viewer/fits-cutout'
               f'?ra={ra_deg}&dec={dec_deg}&layer=decaps2&pixscale={pixscale}'
               f'&bands={band}&size={size}')
        # requests with an IPv4-only adapter: this network has no IPv6 route
        # and urllib's happy-eyeballs order makes it die with errno 51.
        import requests, socket
        from urllib3.util import connection as _uc
        _orig = _uc.allowed_gai_family
        _uc.allowed_gai_family = lambda: socket.AF_INET
        try:
            resp = requests.get(url, timeout=60)
            resp.raise_for_status()
            with open(path, 'wb') as fh:
                fh.write(resp.content)
            _validate_or_remove(path)
            print(f'  Downloaded DECaPS {band} (legacysurvey): '
                  f'{os.path.basename(path)}')
        except Exception:
            # legacysurvey (NERSC-hosted) is often in maintenance; fall back
            # to NOIRLab Astro Data Lab SIA — DECaPS's home institution.
            _fetch_decaps_datalab(ra_deg, dec_deg, band, fov_arcsec, path)
        finally:
            _uc.allowed_gai_family = _orig
    data_nmgy, wcs = _read(path)
    return data_nmgy * NMGY_TO_UJY, wcs


def _fetch_decaps_datalab(ra_deg, dec_deg, band, fov_arcsec, path):
    """DECaPS DR2 cutout via NOIRLab Astro Data Lab SIA (fallback host)."""
    import pyvo, requests
    svc = pyvo.dal.SIAService('https://datalab.noirlab.edu/sia/decaps_dr2')
    fov_deg = fov_arcsec / 3600.0
    res = svc.search(pos=(ra_deg, dec_deg), size=fov_deg).to_table()
    # image rows for the requested band, prefer stacked/coadd products
    mask = [(str(r.get('obs_bandpass', '')).strip().startswith(band) and
             str(r.get('prodtype', '')) == 'image') for r in res]
    cand = res[np.asarray(mask, dtype=bool)]
    if len(cand) == 0:
        raise RuntimeError(f'Data Lab SIA: no DECaPS {band} image here')
    url = str(cand[0]['access_url'])
    resp = requests.get(url, timeout=120)
    resp.raise_for_status()
    with open(path, 'wb') as fh:
        fh.write(resp.content)
    _validate_or_remove(path)
    print(f'  Downloaded DECaPS {band} (Data Lab): {os.path.basename(path)}')


def fetch_racs_fits(ra_deg, dec_deg, fov_arcsec=DOWNLOAD_FOV_ARCSEC, username=None):
    """RACS-mid Stokes-I cutout via CASDA SODA. Returns (data_mJy_per_beam, wcs).

    Requires a (free) CSIRO OPAL account. First-time setup — run once in a
    terminal so the password goes straight into the macOS keyring:

        python -c "from astroquery.casda import Casda; \
                   Casda().login(username='YOUR_OPAL_EMAIL')"

    After that this function reads the credential from the keyring.
    Set username here or export CASDA_USERNAME.
    """
    path = _cache_path('racsmid', ra_deg, dec_deg, fov_arcsec)
    if not os.path.isfile(path):
        from astroquery.casda import Casda
        from astroquery.utils.tap.core import TapPlus
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        user = username or os.environ.get('CASDA_USERNAME')
        if not user:
            raise RuntimeError(
                'RACS cutouts need CASDA (OPAL) credentials: register free at '
                'https://opal.atnf.csiro.au/ then run the one-time login in a '
                'terminal (see fetch_racs_fits docstring), and set '
                'CASDA_USERNAME=<your OPAL email>.')
        casda = Casda()
        casda.login(username=user)          # password from keyring

        tap = TapPlus(url='https://casda.csiro.au/casda_vo_tools/tap')
        # CASDA filenames carry no low/mid/high tag — select the band by
        # wavelength: RACS-mid spans 1295.5-1439.5 MHz → em_min ≈ 0.208 m.
        # Stokes I restored continuum images only. Order tiles by distance
        # from the source to the tile CENTRE: tile-edge cutouts are blanked
        # (all-NaN), so the most-interior tile is the one that has data.
        q = ("SELECT TOP 5 * FROM ivoa.obscore "
             "WHERE dataproduct_subtype = 'cont.restored.t0' "
             "AND filename LIKE 'image.i.%RACS%' "
             "AND em_min BETWEEN 0.20 AND 0.22 "
             f"AND 1 = CONTAINS(POINT('ICRS', {ra_deg}, {dec_deg}), s_region)")
        tbl = tap.launch_job(q).get_results()
        if len(tbl) == 0:
            raise RuntimeError('no RACS-mid image covers this position')
        # Sort client-side by source-to-tile-centre distance (CASDA ADQL
        # rejects geometry functions in ORDER BY): the most-interior tile
        # is the one whose cutout is least likely to be a blanked edge.
        _c0 = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
        _cc = SkyCoord(ra=np.asarray(tbl['s_ra'], dtype=float)*u.deg,
                       dec=np.asarray(tbl['s_dec'], dtype=float)*u.deg)
        tbl = tbl[np.argsort(_c0.separation(_cc).deg)]

        coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
        got = False
        for k in range(len(tbl)):
            urls  = casda.cutout(tbl[k:k+1], coordinates=coord,
                                 radius=(fov_arcsec/2*1.05)*u.arcsec)
            files = casda.download_files(urls, savedir=CACHE_DIR)
            fits_files = [f for f in files if f.endswith('.fits')]
            if not fits_files:
                continue
            try:
                _validate_or_remove(fits_files[0])
            except RuntimeError:
                continue
            d, _ = _read(fits_files[0])
            if np.isfinite(d).mean() < 0.2:      # blanked tile edge
                for f in files:
                    if os.path.isfile(f): os.remove(f)
                print(f'  RACS tile {k} blank at this position, trying next')
                continue
            os.replace(fits_files[0], path)
            for f in files:                      # other cutouts + .checksum
                if f != fits_files[0] and os.path.isfile(f): os.remove(f)
            got = True
            break
        if not got:
            raise RuntimeError('all covering RACS-mid tiles blank here')
        print(f'  Downloaded RACS-mid: {os.path.basename(path)}')
    data_jyb, wcs = _read(path)
    return data_jyb * 1000.0, wcs


def fetch_nvss_fits(ra_deg, dec_deg, fov_arcsec=DOWNLOAD_FOV_ARCSEC):
    """NVSS cutout via SkyView at 15″/px. Returns (data_mJy_per_beam, wcs)."""
    path = _cache_path('nvss', ra_deg, dec_deg, fov_arcsec)
    if not os.path.isfile(path):
        from astroquery.skyview import SkyView
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        imgs = SkyView.get_images(
            position=SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg),
            survey=['NVSS'],
            width=fov_arcsec*u.arcsec, height=fov_arcsec*u.arcsec,
            pixels=_skyview_pixels(fov_arcsec, NVSS_PIXSCALE))
        imgs[0].writeto(path, overwrite=True)
        print(f'  Downloaded NVSS: {os.path.basename(path)}')
    data_jyb, wcs = _read(path)
    return data_jyb * 1000.0, wcs
