"""Wide + zoom optical figure for ONE source (Wan+2025-style layout).

LEFT  — DECaPS r grayscale, 3′ field, SPT 1/2/3σ_pos rings, zoom footprint box
RIGHT — DECaPS g/r/i Lupton RGB, 30″ field, the same rings
Both panels carry a 10″ scale bar and, optionally, a Chandra CSC2 circle.

⚠ Both panels are drawn from HiPS tiles (CDS DECaPS DR2), not calibrated FITS,
so there is deliberately no colorbar: the images show morphology only. This is
an exception to the "panels from calibrated FITS" rule that source_figure.py
follows.
"""
import os

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from astropy.wcs import WCS

from units import HIPS

WIDE_FOV_ARCSEC = 180.0
ZOOM_FOV_ARCSEC = 30.0


def _fetch_hips(hips_id, ra, dec, fov_arcsec, width=300, height=300):
    from astroquery.hips2fits import hips2fits
    hdu = hips2fits.query(hips=hips_id, ra=ra * u.deg, dec=dec * u.deg,
                          fov=fov_arcsec * u.arcsec, width=width, height=height,
                          projection='TAN', format='fits')[0]
    return np.squeeze(hdu.data), WCS(hdu.header).celestial


def decaps_rgb(coord, fov_arcsec, w=300, h=300):
    """DECaPS g/r/i -> (Lupton RGB image, wcs), stretch set from the image.

    stretch = 1.5 × the 95th percentile of the positive mean intensity, which
    keeps faint structure visible without saturating typical stars.
    """
    from astropy.visualization import make_lupton_rgb
    bands, wcs = {}, None
    for b in ('g', 'r', 'i'):
        d, wcs = _fetch_hips(HIPS[f'DECaPS2-{b}'], coord.ra.deg, coord.dec.deg,
                             fov_arcsec, w, h)
        bands[b] = np.nan_to_num(d, nan=0.0)
    intensity = (bands['g'] + bands['r'] + bands['i']) / 3.0
    pos = intensity[intensity > 0]
    p95 = float(np.nanpercentile(pos, 95)) if pos.size else 1.0
    rgb = make_lupton_rgb(bands['i'], bands['r'], bands['g'],
                          Q=10, stretch=max(p95 * 1.5, 1e-3), minimum=0)
    return rgb, wcs


def add_sigma_rings(ax, ra, dec, sigma_arcsec, color, n=3, max_radius_arcsec=None):
    """Circles at k·σ, k = 1..n, skipping any larger than `max_radius_arcsec`."""
    for k in range(1, n + 1):
        r = k * sigma_arcsec
        if max_radius_arcsec and r > max_radius_arcsec:
            continue
        ax.add_patch(SphericalCircle(
            (ra * u.deg, dec * u.deg), r * u.arcsec,
            transform=ax.get_transform('icrs'),
            edgecolor=color, facecolor='none', linewidth=1.3))


def add_circle(ax, ra, dec, r_arcsec, color='royalblue'):
    """One solid circle, e.g. a CSC2 95% error radius."""
    ax.add_patch(SphericalCircle(
        (ra * u.deg, dec * u.deg), r_arcsec * u.arcsec,
        transform=ax.get_transform('icrs'),
        edgecolor=color, facecolor='none', linewidth=1.4))


def add_scale_bar(ax, length_arcsec=10, color='black', label='10″'):
    """Horizontal scale bar in the lower right; pixel scale read from ax.wcs."""
    imgs = ax.get_images()
    if not imgs:
        return
    h, w = imgs[0].get_array().shape[:2]   # (H,W) or (H,W,3)
    c0 = ax.wcs.pixel_to_world(0, h // 2)
    c1 = ax.wcs.pixel_to_world(1, h // 2)
    bar_px = length_arcsec / c0.separation(c1).arcsec
    x_right, y_bot = w * 0.92, h * 0.08
    ax.plot([x_right - bar_px, x_right], [y_bot, y_bot], color=color, linewidth=2.5)
    ax.text(x_right - bar_px / 2, y_bot + h * 0.02, label, ha='center',
            va='bottom', color=color, fontsize=9, fontweight='bold')


def add_zoom_box(ax, coord, zoom_fov_arcsec, color='gray'):
    """Rectangle on the wide panel marking the zoom panel's footprint."""
    half = (zoom_fov_arcsec / 2) / 3600.0
    ra, dec = coord.ra.deg, coord.dec.deg
    dra = half / np.cos(np.radians(dec))
    ax.plot([ra - dra, ra + dra, ra + dra, ra - dra, ra - dra],
            [dec - half, dec - half, dec + half, dec + half, dec - half],
            transform=ax.get_transform('icrs'), color=color, linewidth=1.0)


def _style_axes(ax, ra_fmt, dec_fmt):
    ax.coords['ra'].set_axislabel('R.A.', fontsize=9)
    ax.coords['dec'].set_axislabel('Decl.', fontsize=9)
    ax.coords['ra'].set_major_formatter(ra_fmt)
    ax.coords['dec'].set_major_formatter(dec_fmt)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=10, loc='upper right', borderpad=0.6,
              facecolor='white', framealpha=0.92)


def plot_wide_zoom(source, csc2=None, out_path=None, title_prefix='', dpi=150):
    """Build the wide + zoom figure for one source.

    source : mapping with 'id', 'ra_deg', 'dec_deg', 'snr_max', 'abs_b_deg',
             'sigma_pos_arcsec_heuristic' (a row of v4_crossmatch_table.csv).
    csc2   : optional (ra_deg, dec_deg, r95_arcsec) to overplot.
    Returns (figure or out_path, {panel: error}).
    """
    import matplotlib.pyplot as plt

    sid = source['id']
    coord = SkyCoord(ra=float(source['ra_deg']) * u.deg,
                     dec=float(source['dec_deg']) * u.deg)
    sigma = float(source['sigma_pos_arcsec_heuristic'])

    fig = plt.figure(figsize=(12.5, 6))
    fig.suptitle(f"{title_prefix}{sid}  (SNR={float(source['snr_max']):.1f}, "
                 f"|b|={float(source['abs_b_deg']):.1f}°)", fontsize=11, y=0.99)
    failures = {}

    def overlays(ax, fov_arcsec, ring_color, bar_color):
        add_sigma_rings(ax, coord.ra.deg, coord.dec.deg, sigma, ring_color,
                        max_radius_arcsec=fov_arcsec / 2)
        ax.plot([], [], color=ring_color, linewidth=1.3, label='SPT (1, 2, 3 σ)')
        if csc2 is not None:
            add_circle(ax, csc2[0], csc2[1], csc2[2])
            ax.plot([], [], color='royalblue', linewidth=1.4,
                    label=f'CSC2 95% (r={csc2[2]:.2f}″)')
        add_scale_bar(ax, color=bar_color)

    try:
        data, wcs = _fetch_hips(HIPS['DECaPS2-r'], coord.ra.deg, coord.dec.deg,
                                WIDE_FOV_ARCSEC, 400, 400)
        ax = fig.add_subplot(1, 2, 1, projection=wcs)
        v1, v99 = np.nanpercentile(data, [1, 99])
        ax.imshow(data, origin='lower', cmap='gray_r', vmin=v1, vmax=v99)
        overlays(ax, WIDE_FOV_ARCSEC, 'red', 'black')
        add_zoom_box(ax, coord, ZOOM_FOV_ARCSEC, color='dimgray')
        _style_axes(ax, 'hh:mm:ss', 'dd:mm')
    except Exception as e:
        failures['wide'] = f'{type(e).__name__}: {e}'
        ax = fig.add_subplot(1, 2, 1)
        ax.text(0.5, 0.5, f"wide panel failed\n{failures['wide'][:80]}",
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')

    try:
        rgb, wcs = decaps_rgb(coord, ZOOM_FOV_ARCSEC)
        ax = fig.add_subplot(1, 2, 2, projection=wcs)
        ax.imshow(rgb, origin='lower')
        overlays(ax, ZOOM_FOV_ARCSEC, 'cyan', 'white')
        _style_axes(ax, 'hh:mm:ss.s', 'dd:mm:ss')
    except Exception as e:
        failures['zoom'] = f'{type(e).__name__}: {e}'
        ax = fig.add_subplot(1, 2, 2)
        ax.text(0.5, 0.5, f"zoom panel failed\n{failures['zoom'][:80]}",
                ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if out_path is None:
        return fig, failures
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return out_path, failures
