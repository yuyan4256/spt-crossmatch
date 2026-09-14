"""Nearest-match cone search against one VizieR catalog.

One call = one catalog at one position. The returned entry dict is what
`scripts/build_v4_crossmatch_table.py` and `scripts/add_milliquas.py` turn into
`<catalog>_<field>` columns of outputs/v4_crossmatch_table.csv.

A catalog config is
    {'id': VizieR id, 'ra_col', 'dec_col', 'flux_col', 'flux_unit',
     'extra_cols': [more columns to copy from the best match]}
"""
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

STRICT_THRESH = 10 * u.arcsec     # 'strict' vs 'loose' confidence split


def empty_entry(cfg, search_radius, status='no_source'):
    """The no-match entry: every field present, all values empty."""
    entry = {'match': False, 'match_status': status, 'match_error': None,
             'sep_arcsec': None, 'confidence': None,
             'search_radius_arcsec': float(search_radius.to(u.arcsec).value),
             'flux': None, 'flux_unit': cfg['flux_unit'],
             'cat_ra': None, 'cat_dec': None}
    for ec in cfg.get('extra_cols', []):
        entry[ec] = None
    return entry


def crossmatch_one(cfg, target_coord, search_radius, strip_strings=False):
    """Query one VizieR catalog at one position; keep the nearest source.

    `match_status` is 'matched', 'no_source', 'schema_mismatch' or
    'query_error' — a failed query is recorded, never mistaken for a
    non-detection. Non-numeric extra columns are stored as str
    (whitespace-stripped if `strip_strings`).
    """
    flux_col = cfg['flux_col']
    ra_col, dec_col = cfg['ra_col'], cfg['dec_col']
    extras = cfg.get('extra_cols', [])
    entry = empty_entry(cfg, search_radius)
    try:
        v = Vizier(columns=['*', '+_r'], row_limit=50)
        res = v.query_region(target_coord, radius=search_radius, catalog=cfg['id'])
        if not res:
            return entry
        tbl = next((t for t in res if flux_col in t.colnames), None)
        if tbl is None:
            entry['match_status'] = 'schema_mismatch'
            entry['match_error'] = f'flux_col {flux_col!r} not in any of {len(res)} tables'
            return entry
        if ra_col not in tbl.colnames or dec_col not in tbl.colnames:
            entry['match_status'] = 'schema_mismatch'
            entry['match_error'] = f'declared {ra_col!r}/{dec_col!r} missing'
            return entry
        ra_sample = str(tbl[ra_col][0]).strip()
        if ' ' in ra_sample or 'h' in ra_sample or ':' in ra_sample:
            cat_coords = SkyCoord(ra=list(tbl[ra_col]), dec=list(tbl[dec_col]),
                                  unit=(u.hourangle, u.deg))
        else:
            cat_coords = SkyCoord(ra=np.array(tbl[ra_col], dtype=float) * u.deg,
                                  dec=np.array(tbl[dec_col], dtype=float) * u.deg)
        seps = target_coord.separation(cat_coords)
        best_idx = int(np.argmin(seps))
        best_sep = seps[best_idx]
        if best_sep > search_radius:
            return entry
        entry.update({
            'match': True,
            'match_status': 'matched',
            'sep_arcsec': round(best_sep.to(u.arcsec).value, 2),
            'confidence': 'strict' if best_sep <= STRICT_THRESH else 'loose',
            'flux': float(tbl[flux_col][best_idx]),
            'cat_ra': float(cat_coords[best_idx].ra.deg),
            'cat_dec': float(cat_coords[best_idx].dec.deg),
        })
        for ec in extras:
            if ec in tbl.colnames:
                val = tbl[ec][best_idx]
                try:
                    entry[ec] = float(val)
                except (TypeError, ValueError):
                    entry[ec] = str(val).strip() if strip_strings else str(val)
    except Exception as e:
        entry['match_status'] = 'query_error'
        entry['match_error'] = f'{type(e).__name__}: {e}'
    return entry
