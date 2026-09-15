"""One multi-panel cutout figure for ONE source.

`plot_source()` is the whole per-source figure. It takes an explicit source
record — not a row index into some global DataFrame — so it can be called for
a single source, from a notebook, or from a loop over the whole list without
any hidden state. Running the full list is `scripts/plot_source_panels.py`,
which is nothing but that loop plus argument parsing and error logging.

Panels come from the calibrated FITS in `data/fits_cache/` (2026-07 spec:
plots are drawn from downloaded FITS, never from HiPS preview tiles), so a
figure run is offline once `scripts/precache_fits.py` has been run. Every
colorbar is in a physical unit — mJy, μJy, mJy/beam — never raw DN, and every
panel uses a LINEAR percentile stretch: a colorbar under an asinh stretch
misstates the pixel values it claims to label.

Contour convention (the mentor convention, as in the figures already in
outputs/images/): levels step down from the in-FOV peak by 2σ, floored at 3σ,
at most 6 levels — so a faint source still gets contours and a bright one does
not get 40 of them. σ is a MAD estimate on the outer ring of the *unsmoothed*
SPT cutout; full-cutout MAD inflates it ~3× on the brightest sources because
the matched-filter PSF then dominates the cutout. Dashed contours are the
opposite sign at 3σ/5σ, to make sidelobes and noise blips visible.
"""
import os
import re

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.visualization.wcsaxes import SphericalCircle
from scipy.ndimage import gaussian_filter

from cutouts import load_spt_cutout, reproject_to
from fits_cutouts import DOWNLOAD_FOV_ARCSEC
from paths import OUT, V4_CUTOUT_DIR
from units import unwise_dn_to_mjy

# ─────────────────────────── panel registry ───────────────────────────
# Each entry: fetch(ra_deg, dec_deg) -> (data in the stated unit, wcs).
# Adding a survey to a figure = adding one entry here and naming it in
# `panels=`; nothing else in this module is survey-aware.
MJY_VEGA = r'$F_\nu$ (mJy, Vega)'


def _unwise(band):
    def fetch(ra, dec, source_id=None):
        from unwise import UNWISE_PIXSCALE, fetch_unwise
        size_px = int(np.ceil(DOWNLOAD_FOV_ARCSEC / UNWISE_PIXSCALE))
        data_dn, wcs = fetch_unwise(ra, dec, source_id=source_id,
                                    size_px=size_px, band=int(band[1]))
        return unwise_dn_to_mjy(data_dn, band), wcs
    return fetch


def _allwise(band):
    def fetch(ra, dec, source_id=None):
        from fits_cutouts import fetch_allwise_atlas
        return fetch_allwise_atlas(ra, dec, band)
    return fetch


def _decaps(ra, dec, source_id=None):
    from fits_cutouts import fetch_decaps_fits
    return fetch_decaps_fits(ra, dec)


def _racs(ra, dec, source_id=None):
    from fits_cutouts import fetch_racs_fits
    return fetch_racs_fits(ra, dec)


def _nvss(ra, dec, source_id=None):
    from fits_cutouts import fetch_nvss_fits
    return fetch_nvss_fits(ra, dec)


PANELS = {
    'W1':       dict(title='unWISE W1 (3.4 μm)',   cbar=MJY_VEGA, fetch=_unwise('W1')),
    'W2':       dict(title='unWISE W2 (4.6 μm)',   cbar=MJY_VEGA, fetch=_unwise('W2')),
    'W3':       dict(title='AllWISE W3 (12 μm)',   cbar=MJY_VEGA, fetch=_allwise('W3')),
    'W4':       dict(title='AllWISE W4 (22 μm)',   cbar=MJY_VEGA, fetch=_allwise('W4')),
    'DECaPS-r': dict(title='DECaPS2 r',            cbar=r'$F_\nu$ (μJy, AB)', fetch=_decaps),
    'RACS-mid': dict(title='RACS-mid 1.37 GHz',    cbar='I (mJy/beam)', fetch=_racs),
    'NVSS':     dict(title='NVSS 1.4 GHz',         cbar='I (mJy/beam)', fetch=_nvss),
}

# 2 rows × 3 cols: four small panels on the left, two tall ones on the right.
DEFAULT_SMALL = ('W1', 'W2', 'W3', 'W4')
DEFAULT_BIG = ('DECaPS-r', 'RACS-mid')

SPT_BANDS = (('90GHz', 'blue'), ('150GHz', 'gold'))
STEP_SIGMA = 2.0
MAX_LEVELS = 6
MIN_SIGMA = 3.0
DEFAULT_CUTOUT_TOL = 5.0   # arcsec, id-to-position cutout lookup


# ─────────────────────────── SPT contours ───────────────────────────
def off_source_rms(data, border=6):
    """MAD noise from the outer `border`-wide ring of an SPT cutout.

    The ring excludes the source; full-cutout MAD is biased high (~3×) when a
    strong matched-filter PSF covers a large fraction of the cutout. Falls
    back to full-cutout MAD if the ring has too few finite pixels.
    """
    mask = np.zeros(data.shape, dtype=bool)
    mask[:border, :] = mask[-border:, :] = True
    mask[:, :border] = mask[:, -border:] = True
    ring = data[mask & np.isfinite(data)]
    if ring.size < 50:
        flat = data[np.isfinite(data)]
        return float(1.4826 * np.median(np.abs(flat - np.median(flat))))
    return float(1.4826 * np.median(np.abs(ring - np.median(ring))))


def _parse_spt_id(source_id):
    """(ra_deg, dec_deg) encoded in an SPT3G_JHHMMSS.s+DDMMSS.s name."""
    m = re.match(r'SPT3G_J(\d{2})(\d{2})(\d{2}(?:\.\d+)?)'
                 r'([+-]\d{2})(\d{2})(\d{2}(?:\.\d+)?)', str(source_id))
    if not m:
        return None
    h, mi, s, d, dm, ds = m.groups()
    ra = (int(h) + int(mi) / 60 + float(s) / 3600) * 15.0
    dec = float(d) + np.sign(float(d)) * (int(dm) / 60 + float(ds) / 3600)
    return ra, dec


def find_spt_cutout_id(source_id, ra_deg, dec_deg, cutout_dir=None,
                       tol_arcsec=DEFAULT_CUTOUT_TOL):
    """Which file in `cutout_dir` holds this source's SPT cutout.

    The cutout filenames are named after the centroid run that produced them,
    and data/v4/ currently holds the 3-yr / 61-source run while the source
    list is the 4-yr / 73-source one — the same sky position gets a name that
    differs by a few arcsec, so an exact id lookup finds nothing for any of
    the 73. So: exact id first, then the nearest cutout within `tol_arcsec`.

    Returns (cutout_id, separation_arcsec); (None, None) if nothing is close.
    The separation is returned so the caller can SAY which file it used —
    quietly contouring a source with a neighbour's cutout would be worse than
    drawing no contours at all.
    """
    cutout_dir = cutout_dir or V4_CUTOUT_DIR
    if not os.path.isdir(cutout_dir):
        return None, None
    ids = sorted({f.rsplit('_', 1)[0] for f in os.listdir(cutout_dir)
                  if f.endswith('.fits')})
    if source_id in ids:
        return source_id, 0.0
    best, best_sep = None, np.inf
    for cid in ids:
        pos = _parse_spt_id(cid)
        if pos is None:
            continue
        sep = np.hypot((pos[0] - ra_deg) * np.cos(np.deg2rad(dec_deg)),
                       pos[1] - dec_deg) * 3600.0
        if sep < best_sep:
            best, best_sep = cid, sep
    if best is None or best_sep > tol_arcsec:
        return None, None
    return best, float(best_sep)


def load_spt_contour_set(source_id, ra_deg=None, dec_deg=None, cutout_dir=None,
                         bands=SPT_BANDS, tol_arcsec=DEFAULT_CUTOUT_TOL):
    """({band: {data, wcs, rms, color}}, note) for one source.

    `note` is a short human-readable string for the figure title, saying which
    cutout file the contours came from — empty when the id matched exactly.
    Bands whose FITS file is missing are skipped (no coverage in that band).
    """
    cutout_dir = cutout_dir or V4_CUTOUT_DIR
    cid, sep = source_id, 0.0
    if ra_deg is not None and dec_deg is not None:
        cid, sep = find_spt_cutout_id(source_id, ra_deg, dec_deg,
                                      cutout_dir=cutout_dir,
                                      tol_arcsec=tol_arcsec)
    if cid is None:
        return {}, f'no SPT cutout within {tol_arcsec:.0f}″ — contours omitted'
    out = {}
    for band, color in bands:
        data, wcs = load_spt_cutout(cid, band, cutout_dir=cutout_dir)
        if data is not None:
            out[band] = dict(data=data, wcs=wcs, rms=off_source_rms(data),
                             color=color)
    if not out:
        return {}, 'no SPT cutout — contours omitted'
    note = '' if cid == source_id else f'contours from {cid} ({sep:.1f}″ away)'
    return out, note


def draw_spt_contours(ax, contour_sets, tgt_wcs, tgt_shape, is_neg,
                      source_id='', smooth_spt_px=0.6, label=True):
    """Overlay every SPT band's contours on one panel axis.

    Per-band failures are reported and skipped, never silenced: one band
    without coverage must not cost the whole figure.
    """
    for band, s in contour_sets.items():
        try:
            src = gaussian_filter(s['data'], sigma=smooth_spt_px)
            reproj = reproject_to(src, s['wcs'], tgt_wcs, tgt_shape, order=3)
            finite = reproj[np.isfinite(reproj)]
            if finite.size == 0:
                continue
            rms = s['rms']
            peak_sigma = float(np.nanmax(np.abs(finite))) / rms
            if peak_sigma < MIN_SIGMA:
                continue
            sig = np.array([peak_sigma - i * STEP_SIGMA
                            for i in range(MAX_LEVELS)])
            sig = sig[sig >= MIN_SIGMA]
            if sig.size == 0:
                continue
            sign = -1.0 if is_neg else 1.0
            primary = np.sort(sign * np.sort(sig) * rms)
            cs = ax.contour(reproj, levels=primary, colors=s['color'],
                            linewidths=1.2, linestyles='solid')
            if label:
                ax.clabel(cs, fmt={lv: f'{abs(lv) / rms:.0f}$\\sigma$'
                                   for lv in primary},
                          fontsize=7, inline=True)
            sec = np.sort(-sign * np.array([3.0, 5.0]) * rms)
            ax.contour(reproj, levels=sec, colors=s['color'],
                       linewidths=0.8, linestyles='dashed')
        except Exception as e:
            print(f'  [contour] {source_id} {band}: {type(e).__name__}: {e}')


# ─────────────────────────── one panel ───────────────────────────
def _cut_to_fov(data, wcs, coord, fov_arcsec, margin_px=2):
    """(data, wcs) trimmed to a square of `fov_arcsec` (+ a few pixels) around
    `coord`; pixels outside the download are NaN."""
    from astropy.nddata import Cutout2D
    from astropy.wcs.utils import proj_plane_pixel_scales
    scale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
    n = int(np.ceil(fov_arcsec / scale)) + 2 * margin_px
    if n >= max(data.shape):
        return data, wcs
    cut = Cutout2D(data, coord, (n, n), wcs=wcs, mode='partial', fill_value=np.nan)
    return cut.data, cut.wcs


def _crop(ax, wcs, coord, fov_arcsec):
    """Limit the displayed field to `fov_arcsec`, centred on the source.

    Panels are downloaded once at fits_cutouts.DOWNLOAD_FOV_ARCSEC (786″) at
    native pixel scale; this is where a figure picks how much of that to show.
    A fov beyond the download size just shows everything there is.
    """
    from astropy.wcs.utils import proj_plane_pixel_scales
    scale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
    cx, cy = wcs.world_to_pixel(coord)
    half = (fov_arcsec / 2.0) / scale
    ny, nx = ax.get_images()[0].get_array().shape
    ax.set_xlim(max(cx - half, -0.5), min(cx + half, nx - 0.5))
    ax.set_ylim(max(cy - half, -0.5), min(cy + half, ny - 0.5))


def _draw_panel(fig, gridspec_cell, spec, source, coord, contour_sets,
                markers, fov_arcsec, is_neg, with_legend):
    """Draw one survey panel. Returns None on success, else the error text."""

    sid = source['id']
    try:
        data, wcs = spec['fetch'](float(source['ra_deg']),
                                  float(source['dec_deg']), source_id=sid)
        # Cut the full-size download down to the fov first, so the colour
        # stretch describes the pixels on display and contours are reprojected
        # onto a small grid.
        data, wcs = _cut_to_fov(data, wcs, coord, fov_arcsec)
        if not np.any(np.isfinite(data)):
            raise ValueError('all NaN')
        # Linear percentile stretch: the colorbar has to mean what it says.
        vmin, vmax = np.nanpercentile(data, [1.0, 99.5])
        ax = fig.add_subplot(gridspec_cell, projection=wcs)
        im = ax.imshow(data, origin='lower', cmap='gray_r',
                       vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.08)
        cbar.set_label(spec['cbar'], fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        draw_spt_contours(ax, contour_sets, wcs, data.shape, is_neg, sid)
        _draw_markers(ax, coord, markers)
        _crop(ax, wcs, coord, fov_arcsec)

        ax.coords[0].set_axislabel('R.A.', fontsize=8)
        ax.coords[1].set_axislabel('Decl.', fontsize=8)
        ax.coords[0].set_major_formatter('d.ddd')
        ax.coords[1].set_major_formatter('d.ddd')
        ax.tick_params(labelsize=7)
        ax.set_title(spec['title'], fontsize=10)
        if with_legend:
            for band, color in SPT_BANDS:
                if band in contour_sets:
                    ax.plot([], [], color=color, lw=1.2,
                            label=f'SPT {band.replace("GHz", " GHz")}')
            ax.legend(fontsize=7, loc='upper right', facecolor='white',
                      framealpha=0.85)
        return None
    except Exception as e:
        msg = f'{type(e).__name__}: {e}'
        ax = fig.add_subplot(gridspec_cell)
        ax.set_facecolor('#dde')
        ax.text(0.5, 0.5, f"{spec['title']}\n⚠ {msg[:70]}", ha='center',
                va='center', transform=ax.transAxes, fontsize=8,
                fontweight='bold')
        ax.axis('off')
        return msg


def _draw_markers(ax, coord, markers):
    """SPT centroid + any extra positions (counterparts, error circles)."""
    tr = ax.get_transform('icrs')
    for m in markers:
        ax.plot(m['ra'], m['dec'], transform=tr, marker=m.get('marker', '+'),
                color=m.get('color', 'red'), markersize=m.get('size', 10),
                markeredgewidth=1.5, linestyle='none',
                label=m.get('label'))
        r = m.get('radius_arcsec')
        if r and np.isfinite(r):
            ax.add_patch(SphericalCircle(
                (m['ra'] * u.deg, m['dec'] * u.deg), r * u.arcsec,
                transform=tr, edgecolor=m.get('color', 'red'),
                facecolor='none', linewidth=1.0, linestyle='--'))


def spt_markers(source):
    """The default marker set: SPT centroid + its positional uncertainty."""
    sigma = source.get('sigma_pos_arcsec_heuristic', np.nan)
    m = dict(ra=float(source['ra_deg']), dec=float(source['dec_deg']),
             marker='+', color='red', size=12, label='SPT centroid')
    if sigma is not None and np.isfinite(float(sigma)):
        m['radius_arcsec'] = float(sigma)
        m['label'] = f'SPT centroid (σ={float(sigma):.1f}″)'
    return [m]


def csc2_marker(row):
    """The Chandra position + its 95% positional error circle, or None."""
    ra = float(row.get('CSC2_cat_ra', np.nan))
    dec = float(row.get('CSC2_cat_dec', np.nan))
    if not (np.isfinite(ra) and np.isfinite(dec)):
        return None
    r = np.nanmax([float(row.get('CSC2_r0', np.nan)),
                   float(row.get('CSC2_r1', np.nan))])
    sep = float(row.get('CSC2_sep_arcsec', np.nan))
    label = 'CSC2'
    if np.isfinite(sep):
        label += f' (sep={sep:.1f}″)'
    m = dict(ra=ra, dec=dec, marker='x', color='lime', size=10, label=label)
    if np.isfinite(r):
        m['radius_arcsec'] = r
        m['label'] = label + f', 95% r={r:.2f}″'
    return m


def crossmatch_row(name):
    """One row of outputs/v4_crossmatch_table.csv, looked up by source id."""
    import pandas as pd
    table = pd.read_csv(os.path.join(OUT, 'v4_crossmatch_table.csv'))
    rows = table[table['id'] == name]
    if rows.empty:
        raise KeyError(f'{name!r} is not in v4_crossmatch_table.csv')
    return rows.iloc[0]


def plot_source_by_name(name, fov=45.0, out_path=None):
    """The standard panels for ONE source of the 73, by source id.

    name : SPT source id, e.g. 'SPT3G_J174423.2-311650.6'
    fov  : displayed field of view in arcsec, any value up to 786″, cropped
           offline from FITS cached at 786″ (native pixels; fill the cache with
           scripts/precache_fits.py). Larger values warn and show 786″. SPT
           contours only cover ~570″ (38 × 15″ cutouts).
    out_path : save the PNG here; None returns the figure (notebook use).

    Works for every source, CSC2-matched or not. Markers = SPT centroid +
    σ_pos circle; sources with `has_csc2_match_within_3sigma` also get the
    Chandra position and its 95% circle. A CSC2 source that was recorded but
    fell outside 3σ_pos is NOT marked. Returns (figure or out_path,
    {panel: error}), as `plot_source` does.
    """
    if fov > DOWNLOAD_FOV_ARCSEC:
        import warnings
        warnings.warn(f'fov={fov}″ is larger than the cached {DOWNLOAD_FOV_ARCSEC:.0f}″ '
                      'download; the panels show the full cached field only.')
    source = crossmatch_row(name)
    markers = spt_markers(source)
    if bool(source.get('has_csc2_match_within_3sigma', False)):
        csc2 = csc2_marker(source)
        if csc2:
            markers.append(csc2)
    return plot_source(source, out_path=out_path, markers=markers, fov_arcsec=fov)


# ─────────────────────────── the figure ───────────────────────────
def plot_source(source, out_path=None, small=DEFAULT_SMALL, big=DEFAULT_BIG,
                markers=None, fov_arcsec=45.0, cutout_dir=None, dpi=150,
                figsize=(18, 12), close=True):
    """Build the multi-panel cutout figure for ONE source.

    Parameters
    ----------
    source : mapping
        Needs 'id', 'ra_deg', 'dec_deg', 'snr_max'; uses
        'sigma_pos_arcsec_heuristic' for the uncertainty circle if present.
        A row of outputs/v4_crossmatch_table.csv satisfies this, but so does
        a plain dict — the function never reads a table itself.
    out_path : str or None
        Where to save. None returns the figure without saving (notebook use).
    small, big : sequences of PANELS keys
        `small` fills the left 2×2 block, `big` the two full-height cells on
        the right. Pass fewer keys for a smaller figure.
    markers : list of dicts or None
        Extra positions to mark: {'ra', 'dec', 'label', 'color', 'marker',
        'size', 'radius_arcsec'}. None uses `spt_markers(source)`; pass
        `spt_markers(source) + [...]` to add counterparts to the default.

    Returns
    -------
    (out_path_or_figure, failures) where `failures` is {panel_key: error}
    for the panels that could not be drawn — empty dict means a clean figure.
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    sid = source['id']
    coord = SkyCoord(ra=float(source['ra_deg']) * u.deg,
                     dec=float(source['dec_deg']) * u.deg)
    gal = coord.galactic
    snr = float(source['snr_max'])
    is_neg = snr < 0
    if markers is None:
        markers = spt_markers(source)
    contour_sets, cutout_note = load_spt_contour_set(
        sid, ra_deg=coord.ra.deg, dec_deg=coord.dec.deg, cutout_dir=cutout_dir)

    n_small_rows = 2 if len(small) > 2 else 1
    fig = plt.figure(figsize=figsize)
    fig.suptitle(
        f'{sid}   SNR={snr:.1f},  l={gal.l.deg:.2f}°, b={gal.b.deg:+.2f}°'
        + (f'   [{cutout_note}]' if cutout_note else ''),
        fontsize=12, y=0.99)
    gs = gridspec.GridSpec(max(n_small_rows, len(big)), 3, figure=fig,
                           left=0.05, right=0.97, top=0.93, bottom=0.07,
                           wspace=0.30, hspace=0.30)

    failures = {}
    for i, key in enumerate(small):
        cell = gs[i // 2, i % 2]
        err = _draw_panel(fig, cell, PANELS[key], source, coord, contour_sets,
                          markers, fov_arcsec, is_neg, with_legend=False)
        if err:
            failures[key] = err
    for i, key in enumerate(big):
        err = _draw_panel(fig, gs[i, 2], PANELS[key], source, coord,
                          contour_sets, markers, fov_arcsec, is_neg,
                          with_legend=True)
        if err:
            failures[key] = err

    if out_path is None:
        return fig, failures
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    if close:
        plt.close(fig)
    return out_path, failures
