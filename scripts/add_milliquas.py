"""
Top-up script: add Milliquas v8 (VII/294) cross-match columns to
outputs/v4_crossmatch_table.csv without re-running scripts/build_v4_crossmatch_table.py.

Dual-anchor logic (the build script itself is now SPT-anchor only):
  Pass 1 — SPT-centroid anchor, 60″ search radius
  Pass 2 — CSC2-anchor, 10″ search radius (only for sources with a CSC2 match)

Idempotent: drops any pre-existing Milliquas* columns before re-merging.
"""
import os
import sys
import warnings
from collections import Counter

import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

warnings.filterwarnings('ignore', message='Warning: converting a masked element to nan.')

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, 'src'))

from crossmatch import crossmatch_one  # noqa: E402

CSV_PATH  = os.path.join(REPO_ROOT, 'outputs', 'v4_crossmatch_table.csv')

SPT_ANCHOR_RADIUS  = 60 * u.arcsec     # widened from 30″ → 1′ (user request)
CSC2_ANCHOR_RADIUS = 10 * u.arcsec

MILLIQUAS_CFG = {
    'id': 'VII/294',
    'ra_col': 'RAJ2000',
    'dec_col': 'DEJ2000',
    'flux_col': 'Rmag',
    'flux_unit': 'mag',
    'extra_cols': ['Name', 'Type', 'z', 'Comment'],
}


# ── Load existing crossmatch table ──
final_df = pd.read_csv(CSV_PATH)
print(f'Loaded {len(final_df)} sources × {len(final_df.columns)} cols ← {CSV_PATH}')

# Drop any pre-existing Milliquas columns so re-runs are idempotent.
existing = [c for c in final_df.columns if c.startswith('Milliquas')]
if existing:
    print(f'  dropping {len(existing)} pre-existing Milliquas* cols')
    final_df = final_df.drop(columns=existing)


# ── Pass 1: SPT-centroid anchor ──
spt_rows = []
print(f'\n--- SPT-anchor Milliquas ({MILLIQUAS_CFG["id"]}, r={SPT_ANCHOR_RADIUS}) ---')
for row in final_df.itertuples():
    coord = SkyCoord(ra=row.ra_deg*u.deg, dec=row.dec_deg*u.deg)
    e = crossmatch_one(MILLIQUAS_CFG, coord, SPT_ANCHOR_RADIUS, strip_strings=True)
    e['source_id'] = row.id
    spt_rows.append(e)
spt_df = pd.DataFrame(spt_rows)
print(f'  status: {dict(Counter(r["match_status"] for r in spt_rows))}')


# ── Pass 2: CSC2-anchor ──
print(f'\n--- CSC2-anchor Milliquas (r={CSC2_ANCHOR_RADIUS}) ---')
anchor_rows = []
n_skipped = 0
for row in final_df.itertuples():
    csc2_match = bool(getattr(row, 'CSC2_match', False))
    csc2_ra    = getattr(row, 'CSC2_cat_ra', None)
    csc2_dec   = getattr(row, 'CSC2_cat_dec', None)
    if not csc2_match or csc2_ra is None or pd.isna(csc2_ra):
        e = {'source_id': row.id, 'match': False,
             'match_status': 'no_csc2_anchor', 'match_error': None,
             'sep_arcsec': None, 'confidence': None,
             'search_radius_arcsec': float(CSC2_ANCHOR_RADIUS.to(u.arcsec).value),
             'flux': None, 'flux_unit': MILLIQUAS_CFG['flux_unit'],
             'cat_ra': None, 'cat_dec': None}
        for ec in MILLIQUAS_CFG['extra_cols']:
            e[ec] = None
        anchor_rows.append(e)
        n_skipped += 1
        continue
    anchor = SkyCoord(ra=float(csc2_ra)*u.deg, dec=float(csc2_dec)*u.deg)
    e = crossmatch_one(MILLIQUAS_CFG, anchor, CSC2_ANCHOR_RADIUS, strip_strings=True)
    e['source_id'] = row.id
    anchor_rows.append(e)
anchor_df = pd.DataFrame(anchor_rows)
sc = dict(Counter(r['match_status'] for r in anchor_rows))
print(f'  status: {sc}  (skipped {n_skipped}/{len(final_df)} for no CSC2 anchor)')


# ── Merge into final_df, mirroring notebook column conventions ──
SPT_KEEP_BASE    = ['source_id', 'match', 'match_status', 'match_error',
                    'sep_arcsec', 'confidence', 'search_radius_arcsec',
                    'flux', 'flux_unit', 'cat_ra', 'cat_dec']
ANCHOR_KEEP_BASE = ['source_id', 'match', 'match_status',
                    'sep_arcsec', 'confidence', 'flux', 'flux_unit',
                    'cat_ra', 'cat_dec']
extras = MILLIQUAS_CFG['extra_cols']

cols_spt = [c for c in (SPT_KEEP_BASE + extras) if c in spt_df.columns]
spt2 = spt_df[cols_spt].rename(
    columns={c: f'Milliquas_{c}' for c in cols_spt if c != 'source_id'})
final_df = final_df.merge(spt2, left_on='id', right_on='source_id',
                          how='left').drop(columns='source_id', errors='ignore')

cols_anc = [c for c in (ANCHOR_KEEP_BASE + extras) if c in anchor_df.columns]
anc2 = anchor_df[cols_anc].rename(
    columns={c: f'Milliquas_csc2anchor_{c}' for c in cols_anc if c != 'source_id'})
final_df = final_df.merge(anc2, left_on='id', right_on='source_id',
                          how='left').drop(columns='source_id', errors='ignore')


# ── Save ──
final_df.to_csv(CSV_PATH, index=False)
print(f'\nSaved {len(final_df)} rows × {len(final_df.columns)} cols → {CSV_PATH}')


# ── Summary ──
n_spt = int(final_df['Milliquas_match'].fillna(False).astype(bool).sum())
n_anc = int(final_df['Milliquas_csc2anchor_match'].fillna(False).astype(bool).sum())
print(f'\n=== Milliquas summary ===')
print(f'  SPT-anchor matches : {n_spt}/{len(final_df)}')
print(f'  CSC2-anchor matches: {n_anc}/{len(final_df)}')

# Blazar / BL Lac candidates — any Milliquas Type code that contains 'B'
# (after stripping leading whitespace). Includes BR, B, 2B, etc.
def is_blazar(t):
    if not isinstance(t, str) or t.strip() == '' or t.strip().lower() == 'nan':
        return False
    return 'B' in t.strip().upper()

bl_mask = final_df['Milliquas_Type'].apply(is_blazar)
if bl_mask.any():
    print(f'\nBL Lac / blazar candidates ({int(bl_mask.sum())}):')
    show = ['id', 'snr_max', 'abs_b_deg', 'Milliquas_sep_arcsec',
            'Milliquas_Name', 'Milliquas_Type', 'Milliquas_z',
            'Milliquas_Comment', 'simbad_otype']
    show = [c for c in show if c in final_df.columns]
    print(final_df[bl_mask][show].to_string(index=False))
else:
    print('\nNo BL Lac / blazar candidates flagged (no Type contains "B").')

if n_spt > 0:
    print(f'\nAll Milliquas matches:')
    show = ['id', 'snr_max', 'Milliquas_sep_arcsec', 'Milliquas_Name',
            'Milliquas_Type', 'Milliquas_z', 'Milliquas_Comment']
    show = [c for c in show if c in final_df.columns]
    print(final_df[final_df['Milliquas_match'].fillna(False).astype(bool)][show].to_string(index=False))


# ── Derive focused Milliquas summary CSV (for upload / sharing) ──
summary_cols = [
    # source context
    'id', 'ra_deg', 'dec_deg', 'snr_max', 'abs_b_deg',
    'sigma_pos_arcsec_heuristic',
    'has_csc2_match_within_3sigma',
    # independent SIMBAD classification (sanity cross-check)
    'simbad_name', 'simbad_otype',
    # SPT-anchor Milliquas
    'Milliquas_match', 'Milliquas_match_status',
    'Milliquas_sep_arcsec', 'Milliquas_search_radius_arcsec',
    'Milliquas_confidence',
    'Milliquas_Name', 'Milliquas_Type', 'Milliquas_z', 'Milliquas_Comment',
    'Milliquas_flux',  # Rmag
    'Milliquas_cat_ra', 'Milliquas_cat_dec',
    # CSC2-anchor Milliquas
    'Milliquas_csc2anchor_match', 'Milliquas_csc2anchor_match_status',
    'Milliquas_csc2anchor_sep_arcsec', 'Milliquas_csc2anchor_confidence',
    'Milliquas_csc2anchor_Name', 'Milliquas_csc2anchor_Type',
    'Milliquas_csc2anchor_z', 'Milliquas_csc2anchor_Comment',
]
summary_cols = [c for c in summary_cols if c in final_df.columns]
summary = final_df[summary_cols].rename(columns={'Milliquas_flux': 'Milliquas_Rmag'}).copy()
# Sort: matched first, then by |SNR| descending
summary['_match'] = summary['Milliquas_match'].fillna(False).astype(bool)
summary = summary.sort_values(
    ['_match', 'snr_max'],
    ascending=[False, False],
    key=lambda c: c.abs() if c.name == 'snr_max' else c
).drop(columns='_match')

summary_path = os.path.join(REPO_ROOT, 'outputs', 'v4_milliquas_summary.csv')
summary.to_csv(summary_path, index=False)
print(f'\nSaved summary {len(summary)} rows × {len(summary.columns)} cols → {summary_path}')
