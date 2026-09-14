"""Figure helpers that are not tied to one source: the coadd cutout grid.

Per-source SPT contours and markers live in `source_figure.py`.
"""
import numpy as np
import matplotlib.pyplot as plt


# ─────────── coadd cutout grid ───────────
def coadd_grid(maps_dict, df, x_col='xcentroid', y_col='ycentroid',
               id_col='id', extra_label_col='snr_max',
               half=25, vlim=1000, n_cols=8, cmap='RdBu_r',
               unit_divisor=1e-3, unit_label='μK'):
    """Plot a grid of cutouts from one coadd map.

    Parameters
    ----------
    maps_dict : single 2-D array or dict of {band: array}
        If dict, returns one figure per band (returns dict of figs).
    df : DataFrame with x_col, y_col, id_col, optionally extra_label_col
    unit_divisor : data is divided by this (default 1e-3 = G3Units.uK)
    """

    if isinstance(maps_dict, dict):
        return {b: _single_grid(arr, df, x_col, y_col, id_col,
                                extra_label_col, half, vlim, n_cols, cmap,
                                unit_divisor, unit_label, band_label=b)
                for b, arr in maps_dict.items()}
    return _single_grid(maps_dict, df, x_col, y_col, id_col,
                        extra_label_col, half, vlim, n_cols, cmap,
                        unit_divisor, unit_label)


def _single_grid(arr, df, x_col, y_col, id_col, extra_col, half, vlim,
                 n_cols, cmap, unit_divisor, unit_label, band_label=None):
    from cutouts import coadd_cutout

    n = len(df)
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.0 * n_cols, 2.2 * n_rows))
    axes = np.atleast_2d(axes).ravel()
    for i, (_, row) in enumerate(df.iterrows()):
        cut = coadd_cutout(arr, row[x_col], row[y_col], half=half) / unit_divisor
        ax = axes[i]
        ax.imshow(cut, vmin=-vlim, vmax=vlim, origin='lower', cmap=cmap)
        ax.axhline(half, color='k', lw=0.3)
        ax.axvline(half, color='k', lw=0.3)
        sid = str(row[id_col]).split('_')[-1] if '_' in str(row[id_col]) else str(row[id_col])
        if extra_col and extra_col in row:
            title = f"{sid}\n{extra_col}={row[extra_col]:.0f}"
        else:
            title = sid
        ax.set_title(title, fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
    for j in range(n, len(axes)):
        axes[j].axis('off')
    suptitle = f'cutouts (±{vlim} {unit_label})'
    if band_label:
        suptitle = f'{band_label} ' + suptitle
    fig.suptitle(suptitle, y=1.0)
    fig.tight_layout()
    return fig
