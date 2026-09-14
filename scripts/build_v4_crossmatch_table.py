#!/usr/bin/env python
"""Build the multi-catalog cross-match table from the 4-yr source list.

Input : data/centroids 4_yr.cat  (73 sources)
Steps : 1. per-source sigma_pos, Wan+2025 Eqs 4-6 (src/positions.py)
        2. SIMBAD, nearest object within 30″
        3. VizieR catalogs at the SPT centroid (src/crossmatch.py):
           CSC2 search radius = 3·sigma_pos per source, every other catalog 30″
        4. has_csc2_match_within_3sigma; csc2_covered = 'Unknown' (no CSC2
           footprint API yet — a non-match is NOT an X-ray non-detection)
Output: outputs/v4_crossmatch_table.csv

This is the FIRST step of the table's life, not the whole of it. Later
scripts add or rewrite columns in place — add_milliquas.py, add_decaps.py,
update_sigma_wan2025.py, racs_pchance/fetch_racs_agn.py, wise_pchance/score_wise.py
— so re-running this alone gives a table without those columns. That is why
writing over the existing table needs --overwrite.

SIMBAD via astroquery>=0.4.8 returns `simbad_otype` as short codes ('BLL',
'Rad', 'EB*'); the existing table holds the old long labels ('BLLac', 'Radio',
'EclBin'). Code that filters on otype must accept both (see
plot_spt_blazar_overlay.py). Checked 2026-09-14: same 59 SIMBAD matches, same
names, separations within 0.45″; all VizieR columns identical to the table.

NOT a classification table: a CSC2 match within 3σ_pos is necessary but not
sufficient for a verified X-ray counterpart.

Converted from notebooks/build_v4_crossmatch_table.ipynb (2026-09-14); the
logic is unchanged, the shared pieces moved to src/.

Usage:
  python scripts/build_v4_crossmatch_table.py --out /tmp/v4_check.csv
  python scripts/build_v4_crossmatch_table.py --overwrite
"""
import argparse
import os
import re
import sys
from collections import Counter

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'src'))

from crossmatch import STRICT_THRESH, crossmatch_one
from galactic import add_galactic_coords, bin_counts
from positions import sigma_pos_wan2025

SOURCE_LIST = os.path.join(REPO, 'data', 'centroids 4_yr.cat')
TABLE = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')

LOOSE_THRESH = 30 * u.arcsec

# ra_col / dec_col explicit because CSC2 uses RAICRS/DEICRS, not RAJ2000.
VIZIER_CATALOGS = {
    # Radio
    'AT20G':    {'id': 'J/MNRAS/384/775',            'ra_col': 'RAJ2000', 'dec_col': 'DEJ2000', 'flux_col': 'S20',   'flux_unit': 'Jy',        'extra_cols': []},
    'NVSS':     {'id': 'VIII/65',                    'ra_col': 'RAJ2000', 'dec_col': 'DEJ2000', 'flux_col': 'S1.4',  'flux_unit': 'mJy',       'extra_cols': []},
    'RACS-mid': {'id': 'J/other/PASA/41.3/sourcesm', 'ra_col': 'RAJ2000', 'dec_col': 'DEJ2000', 'flux_col': 'Ftot',  'flux_unit': 'mJy',       'extra_cols': []},
    'VLASS':    {'id': 'J/ApJ/914/42/table5',        'ra_col': 'RAJ2000', 'dec_col': 'DEJ2000', 'flux_col': 'Flux',  'flux_unit': 'Jy',        'extra_cols': []},
    # MIR — extra cols needed for Stern+2012 AGN selection (W1, W2, quality)
    'AllWISE':  {'id': 'II/328',                     'ra_col': 'RAJ2000', 'dec_col': 'DEJ2000', 'flux_col': 'W1mag', 'flux_unit': 'mag',       'extra_cols': ['W2mag', 'W3mag', 'W4mag', 'e_W2mag', 'ccf', 'ex', 'qph']},
    # X-ray (CSC2.1 master); r0/r1 are 95% error-ellipse radii in arcsec
    'CSC2':     {'id': 'IX/57',                      'ra_col': 'RAICRS',  'dec_col': 'DEICRS',  'flux_col': 'Fluxb', 'flux_unit': 'erg/cm2/s', 'extra_cols': ['r0', 'r1']},
    # Milliquas v8 (Flesch 2023). Type: Q=QSO, A=AGN, B=BL Lac, K=NLQSO,
    # N=NLAGN; suffix R/X = radio/X-ray associated; prefix '2' = radio lobe.
    'Milliquas': {'id': 'VII/294',                   'ra_col': 'RAJ2000', 'dec_col': 'DEJ2000', 'flux_col': 'Rmag',  'flux_unit': 'mag',       'extra_cols': ['Name', 'Type', 'z', 'Comment']},
}
# DECaPS2 is not queried here (VizieR II/367 is VHS DR5, not DECaPS);
# see add_decaps.py.

SIMBAD_KEEP = ['source_id', 'simbad_match', 'simbad_status', 'simbad_error',
               'simbad_name', 'simbad_otype', 'simbad_ra', 'simbad_dec',
               'simbad_sep_arcsec', 'simbad_z', 'simbad_confidence']
VIZIER_KEEP = ['source_id', 'match', 'match_status', 'match_error',
               'sep_arcsec', 'confidence', 'search_radius_arcsec',
               'flux', 'flux_unit', 'cat_ra', 'cat_dec']


def load_sources(path=SOURCE_LIST):
    """The pipe-separated centroid catalog -> DataFrame, brightest |SNR| first."""
    df = pd.read_csv(path, sep=r'\s*\|\s*', engine='python')
    df = df.loc[:, ~df.columns.str.match('^Unnamed')]
    df.columns = df.columns.str.strip()
    for col in df.select_dtypes('object').columns:
        df[col] = df[col].str.strip()
    df['ra_deg'] = df['sky_centroid'].apply(lambda x: float(x.split(',')[0]))
    df['dec_deg'] = df['sky_centroid'].apply(lambda x: float(x.split(',')[1]))
    add_galactic_coords(df)
    return df.sort_values('snr_max', key=abs, ascending=False).reset_index(drop=True)


def simbad_one(simbad, source_id, coord):
    """Nearest SIMBAD object within LOOSE_THRESH, as one entry dict."""
    entry = {'source_id': source_id, 'simbad_match': False,
             'simbad_status': 'no_source', 'simbad_error': None,
             'simbad_name': None, 'simbad_otype': None,
             'simbad_ra': None, 'simbad_dec': None,
             'simbad_sep_arcsec': None, 'simbad_z': None,
             'simbad_confidence': None}
    try:
        result = simbad.query_region(coord, radius=LOOSE_THRESH)
        if result is None or len(result) == 0:   # astroquery>=0.4.8 returns an empty table
            return entry
        cols = result.colnames
        ra_col = next((c for c in cols if c.upper() == 'RA'), None)
        dec_col = next((c for c in cols if c.upper() == 'DEC'), None)
        id_col = next((c for c in cols if 'MAIN_ID' in c.upper()), None)
        ot_col = next((c for c in cols if 'OTYPE' in c.upper()), None)
        z_col = next((c for c in cols if c.upper() in ('Z_VALUE', 'RVZ_REDSHIFT')), None)
        if ra_col is None or dec_col is None:
            entry['simbad_status'] = 'schema_mismatch'
            entry['simbad_error'] = f'missing RA/DEC col; got {list(cols)[:6]}'
            return entry
        # astroquery<0.4.8 (SIMBAD script API) gave sexagesimal RA strings;
        # the TAP-based client gives degrees with a unit on the column.
        if result[ra_col].unit is not None and u.Unit(result[ra_col].unit) == u.deg:
            sim_coords = SkyCoord(ra=np.asarray(result[ra_col], dtype=float) * u.deg,
                                  dec=np.asarray(result[dec_col], dtype=float) * u.deg)
        else:
            sim_coords = SkyCoord(ra=list(result[ra_col]), dec=list(result[dec_col]),
                                  unit=(u.hourangle, u.deg))
        seps = coord.separation(sim_coords)
        i = int(np.argmin(seps))
        entry.update({
            'simbad_match': True,
            'simbad_status': 'matched',
            'simbad_name': result[id_col][i] if id_col else None,
            'simbad_otype': result[ot_col][i] if ot_col else None,
            'simbad_ra': round(sim_coords[i].ra.deg, 6),
            'simbad_dec': round(sim_coords[i].dec.deg, 6),
            'simbad_sep_arcsec': round(seps[i].to(u.arcsec).value, 2),
            'simbad_z': result[z_col][i] if z_col else None,
            'simbad_confidence': 'strict' if seps[i] <= STRICT_THRESH else 'loose',
        })
    except Exception as e:
        entry['simbad_status'] = 'query_error'
        entry['simbad_error'] = f'{type(e).__name__}: {e}'
        print(f'    ERROR {type(e).__name__}: {e}')
    return entry


def build(df):
    from astroquery.simbad import Simbad

    coords = SkyCoord(ra=df['ra_deg'].values * u.deg, dec=df['dec_deg'].values * u.deg)
    df['sigma_pos_arcsec_heuristic'] = df['snr_max'].apply(sigma_pos_wan2025)
    print(f"sigma_pos median {df['sigma_pos_arcsec_heuristic'].median():.2f}″ "
          f"(3σ CSC2 radius {3 * df['sigma_pos_arcsec_heuristic'].min():.1f}–"
          f"{3 * df['sigma_pos_arcsec_heuristic'].max():.1f}″)")

    simbad = Simbad()
    simbad.reset_votable_fields()
    simbad.add_votable_fields('ra', 'dec', 'otype', 'rvz_redshift')
    simbad_rows = []
    for sid, coord in zip(df['id'], coords):
        print(f'  SIMBAD {sid} ...')
        simbad_rows.append(simbad_one(simbad, sid, coord))
    simbad_df = pd.DataFrame(simbad_rows)
    print('SIMBAD status:', simbad_df['simbad_status'].value_counts().to_dict())

    out = df.merge(simbad_df[[c for c in SIMBAD_KEEP if c in simbad_df.columns]],
                   left_on='id', right_on='source_id', how='left'
                   ).drop(columns='source_id', errors='ignore')

    for cat, cfg in VIZIER_CATALOGS.items():
        print(f"\n--- {cat} ({cfg['id']}) ---")
        rows = []
        for row, coord in zip(df.itertuples(), coords):
            radius = (3 * row.sigma_pos_arcsec_heuristic * u.arcsec
                      if cat == 'CSC2' else LOOSE_THRESH)
            e = crossmatch_one(cfg, coord, radius)
            e['source_id'] = row.id
            rows.append(e)
        print(f"  status: {dict(Counter(r['match_status'] for r in rows))}")
        vdf = pd.DataFrame(rows)
        cols = [c for c in VIZIER_KEEP + cfg['extra_cols'] if c in vdf.columns]
        vdf = vdf[cols].rename(columns={c: f'{cat}_{c}' for c in cols if c != 'source_id'})
        out = out.merge(vdf, left_on='id', right_on='source_id', how='left'
                        ).drop(columns='source_id', errors='ignore')

    out['has_csc2_match_within_3sigma'] = out['CSC2_match'].fillna(False).astype(bool)
    out['csc2_covered'] = 'Unknown'
    # Stern+2012 W1-W2 (Vega); needs quality/depth cuts at low |b|
    out['W1_W2'] = out['AllWISE_flux'] - out['AllWISE_W2mag']

    # Fold flux units into column names; drop unit / bin / Comment / empty cols.
    for cat, cfg in VIZIER_CATALOGS.items():
        unit = re.sub(r'[^A-Za-z0-9]+', '_', str(cfg['flux_unit'])).strip('_')
        out = out.rename(columns={f'{cat}_flux': f'{cat}_flux_{unit}'})
        out = out.drop(columns=[f'{cat}_flux_unit'], errors='ignore')
    out = out.drop(columns=['abs_b_bin'], errors='ignore')
    out = out.drop(columns=[c for c in out.columns if 'Comment' in c])
    return out.dropna(axis=1, how='all')


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default=TABLE)
    ap.add_argument('--overwrite', action='store_true',
                    help='allow replacing an existing output file')
    args = ap.parse_args()
    if os.path.exists(args.out) and not args.overwrite:
        raise SystemExit(f'{args.out} exists; later scripts have added columns to it. '
                         'Pass --out elsewhere, or --overwrite if you mean it.')

    df = load_sources()
    print(f'Loaded {len(df)} sources from {os.path.basename(SOURCE_LIST)}')
    print(bin_counts(df).to_string())
    out = build(df)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, index=False)
    n_csc2 = int(out['has_csc2_match_within_3sigma'].sum())
    print(f'\nSaved {len(out)} rows × {len(out.columns)} columns -> {args.out}')
    print(f'has_csc2_match_within_3sigma: {n_csc2} / {len(out)} '
          "(csc2_covered = 'Unknown' for all rows)")


if __name__ == '__main__':
    main()
