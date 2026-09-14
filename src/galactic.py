"""Galactic coords + |b| binning. b<10 deg = SPT Galactic-plane field."""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

ABS_B_BINS = [0, 1, 3, 5, 6, 10, 90]
ABS_B_LABELS = ['<1 deg', '1-3 deg', '3-5 deg', '5-6 deg', '6-10 deg', '>10 deg']


def add_galactic_coords(df, ra_col='ra_deg', dec_col='dec_deg'):
    """Add l_deg, b_deg, abs_b_deg, abs_b_bin to df in place."""
    gc = SkyCoord(ra=df[ra_col].values * u.deg,
                  dec=df[dec_col].values * u.deg,
                  frame='icrs').galactic
    df['l_deg'] = gc.l.deg
    df['b_deg'] = gc.b.deg
    df['abs_b_deg'] = np.abs(gc.b.deg)
    df['abs_b_bin'] = pd.cut(df['abs_b_deg'], bins=ABS_B_BINS,
                             labels=ABS_B_LABELS, right=False)
    return df


def bin_counts(df):
    """Return source counts per |b| bin.

    Falls back to computing bins from abs_b_deg when the abs_b_bin column
    is absent (it is dropped from the saved crossmatch CSV).
    """
    if 'abs_b_bin' in df.columns:
        s = df['abs_b_bin']
    else:
        s = pd.cut(df['abs_b_deg'], bins=ABS_B_BINS,
                   labels=ABS_B_LABELS, right=False)
    return s.value_counts().reindex(ABS_B_LABELS, fill_value=0)
