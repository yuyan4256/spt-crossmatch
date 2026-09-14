"""Unit conversions + colorbar labels.

Hard rule: never display raw detector/image units on a colorbar. Convert to a
physical unit (μJy, mJy, mag, μK) before plotting. See `unwise_dn_to_*` and
`cbar_label`.
"""
import numpy as np

# ─────────── unWISE (W1 / W2) ───────────
# unWISE coadds are in **Vega** nanomaggies. The zeropoint relation
#   m = 22.5 - 2.5 log10(DN)
# does NOT by itself imply AB — the magnitude system comes from the data
# product docs, and for unWISE it is VEGA (unwise.me).
#
# Vega zero-magnitude flux densities (WISE Explanatory Supplement):
#   W1: 309.540 Jy    W2: 171.787 Jy
# → 1 DN = F0 × 1e-9 Jy   (W1: 0.30954 μJy, W2: 0.171787 μJy)
#
# HISTORY (2026-07-29): this block previously applied the AB-nanomaggy
# factor 3.631 μJy/DN to unWISE, overestimating W1 flux ~12×. The AB
# factor is only correct for AB-nanomaggy products (e.g. Legacy Survey /
# DECaPS — see src/fits_cutouts.NMGY_TO_UJY).
UNWISE_MAGZP = 22.5                                  # VEGA system
UNWISE_F0_JY     = {'W1': 309.540,    'W2': 171.787}
UNWISE_DN_TO_UJY = {'W1': 309.540e-3, 'W2': 171.787e-3}   # μJy per DN
AB_TO_VEGA = {'W1': 2.699, 'W2': 3.339}              # m_AB = m_Vega + offset
# Saturation: ~50,000 DN → W1 ≈ 15.5 mJy ≈ Vega 10.7 mag
UNWISE_SATURATION_DN = 50_000

def unwise_dn_to_ujy(dn, band):
    """unWISE DN → μJy (Vega-based flux). `band` is 'W1' or 'W2'."""
    return np.asarray(dn, dtype=float) * UNWISE_DN_TO_UJY[band]

def unwise_dn_to_mjy(dn, band):
    """unWISE DN → mJy. `band` is 'W1' or 'W2'."""
    return unwise_dn_to_ujy(dn, band) * 1e-3

def unwise_dn_to_vega_mag(dn):
    """unWISE DN → Vega mag (native system). Returns NaN where dn ≤ 0."""
    dn = np.asarray(dn, dtype=float)
    out = np.full_like(dn, np.nan)
    pos = dn > 0
    out[pos] = UNWISE_MAGZP - 2.5 * np.log10(dn[pos])
    return out

def unwise_dn_to_ab_mag(dn, band):
    """unWISE DN → AB mag (Vega mag + band offset). `band` is 'W1' or 'W2'."""
    return unwise_dn_to_vega_mag(dn) + AB_TO_VEGA[band]


# ─────────── SPT 3G coadd ───────────
# FITS header UNITS=Tcmb. Raw values are mK_CMB (base unit = mK).
# core.G3Units.uK = 1e-3 in this dataset → divide by 1e-3 to get μK.
G3UNITS_UK = 1e-3

def spt_to_uk(data):
    """SPT coadd raw values (mK_CMB) → μK_CMB."""
    return np.asarray(data, dtype=float) / G3UNITS_UK


# ─────────── Colorbar labels ───────────
# Use this dict so every plot ends up with a physical unit.
# Key = (survey, unit). Pick the physical unit you want; never request raw
# image units such as DN/counts/ADU.
CBAR_LABEL = {
    ('unWISE-W1', 'mJy'):  'F$_\\nu$ (mJy, Vega)',
    ('unWISE-W1', 'μJy'):  'F$_\\nu$ (μJy, Vega)',
    ('unWISE-W1', 'mag'):  'W1 (Vega mag)',
    ('unWISE-W2', 'mJy'):  'F$_\\nu$ (mJy, Vega)',
    ('unWISE-W2', 'mag'):  'W2 (Vega mag)',
    ('DECaPS2',   'mag'):  'r (mag)',
    ('RACS-mid',  'mJy/beam'): 'I (mJy/beam)',
    ('NVSS',      'mJy/beam'): 'I (mJy/beam)',
    ('VLASS',     'mJy/beam'): 'I (mJy/beam)',
    ('SPT-90',    'μK'):   'T (μK$_{CMB}$)',
    ('SPT-150',   'μK'):   'T (μK$_{CMB}$)',
}

RAW_IMAGE_UNITS = {'dn', 'count', 'counts', 'adu', 'data number', 'data numbers'}

def cbar_label(survey, unit):
    """Return matplotlib-ready colorbar label.

    Raises if a raw detector/image unit is requested.
    """
    if unit is None:
        return None
    if str(unit).strip().lower() in RAW_IMAGE_UNITS:
        raise ValueError(
            f'{unit!r} is not a physical unit. Convert to μJy/mJy/mag/μK '
            'before adding a colorbar.')
    key = (survey, unit)
    if key not in CBAR_LABEL:
        raise KeyError(f'No colorbar label registered for {key!r}. '
                       f'Add it to CBAR_LABEL in src/units.py.')
    return CBAR_LABEL[key]


# ─────────── HiPS URLs (kept here so all plotting code uses one source) ───────────
HIPS = {
    'DECaPS2-g':  'CDS/P/DECaPS/DR2/g',
    'DECaPS2-r':  'CDS/P/DECaPS/DR2/r',
    'DECaPS2-i':  'CDS/P/DECaPS/DR2/i',
    'DECaPS2-z':  'CDS/P/DECaPS/DR2/z',
    # Registered HiPS ID (not the alasky CASDA mirror URL): the mirror has
    # stale/empty tiles at some Galactic-plane positions; the registered ID
    # resolves to a working backend.
    'RACS-mid':   'CSIRO/P/RACS/mid/I',
    'NVSS':       'CDS/P/NVSS',
    'VLASS-epoch1': 'CDS/P/VLASS/epoch1',
}
