"""
Optical (DECaPS2 = DECam Plane Survey 2) CATALOG cross-match for the 61 SPT
variable sources, via NOIRLab Astro Data Lab TAP (table decaps_dr2.object).

Makes 'optical' a real catalog match (previously DECaPS was only used for image
cutouts). For each source: nearest DECaPS2 object (sep + g/r/i/z/y mags), the
LOCAL source density, and the chance-alignment probability within the SPT
positional error -- because on this Galactic-plane field DECaPS2 is so deep/dense
that a bare 'nearest optical source' is almost always a foreground coincidence.

Needs pyvo. Saves outputs/decaps_crossmatch.csv.
"""
import os, numpy as np, pandas as pd, pyvo
from astropy.coordinates import SkyCoord
import astropy.units as u
import warnings; warnings.filterwarnings('ignore')

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(REPO, 'outputs', 'v4_crossmatch_table.csv')
OUT  = os.path.join(REPO, 'outputs', 'decaps_crossmatch.csv')
TAP  = pyvo.dal.TAPService('https://datalab.noirlab.edu/tap')
R    = 15.0                       # arcsec: query box half-size + density aperture
MAGS = ['mean_mag_g', 'mean_mag_r', 'mean_mag_i', 'mean_mag_z', 'mean_mag_y']

src = pd.read_csv(SRC)
rows = []
for i, s in enumerate(src.itertuples()):
    ra0, dec0 = s.ra_deg, s.dec_deg
    sig = float(s.sigma_pos_arcsec_heuristic)
    dd = R/3600.0; dra = dd/np.cos(np.radians(dec0))
    q = (f"SELECT ra,dec,{','.join(MAGS)} FROM decaps_dr2.object "
         f"WHERE ra BETWEEN {ra0-dra} AND {ra0+dra} AND dec BETWEEN {dec0-dd} AND {dec0+dd}")
    rec = {'id': s.id, 'sigma_pos_arcsec': round(sig, 2)}
    try:
        d = TAP.search(q, maxrec=6000).to_table().to_pandas()
    except Exception as e:
        rec['decaps_error'] = f'{type(e).__name__}'; rows.append(rec)
        print(f'{i+1}/61 {s.id}  ERROR {type(e).__name__}'); continue
    if len(d) == 0:
        rec.update({'decaps_match': False, 'decaps_sep_arcsec': np.nan,
                    'decaps_n_within_R': 0, 'decaps_density_per_as2': 0.0, 'decaps_p_chance': 0.0})
        for m in MAGS: rec[m] = np.nan
    else:
        sep = SkyCoord(ra0*u.deg, dec0*u.deg).separation(
              SkyCoord(d.ra.values*u.deg, d.dec.values*u.deg)).arcsec
        d['sep'] = sep
        inR = d[d.sep <= R]
        n = len(inR); rho = n/(np.pi*R*R)                 # sources / arcsec^2
        pch = 1 - np.exp(-rho*np.pi*sig*sig)              # P(>=1 random src within sigma_pos)
        near = d.loc[d.sep.idxmin()]
        rec.update({'decaps_match': bool(near.sep <= R),
                    'decaps_sep_arcsec': round(float(near.sep), 2),
                    'decaps_n_within_R': int(n),
                    'decaps_density_per_as2': round(float(rho), 4),
                    'decaps_p_chance': round(float(pch), 3)})
        for m in MAGS:
            v = float(near[m]); rec[m] = np.nan if not np.isfinite(v) else round(v, 2)
    rows.append(rec)
    print(f'{i+1}/61 {s.id}  sep={rec.get("decaps_sep_arcsec")}"  n<={R:.0f}\"={rec.get("decaps_n_within_R")}  '
          f'P_chance={rec.get("decaps_p_chance")}')

out = pd.DataFrame(rows); out.to_csv(OUT, index=False)
ok = out[out.decaps_match == True] if 'decaps_match' in out else out.iloc[:0]
print(f'\n=== DECaPS2 optical cross-match (R={R:.0f}") ===')
print(f'match (optical src within {R:.0f}"): {int(out.get("decaps_match", pd.Series(dtype=bool)).sum())}/61')
print(f'nearest sep: median {ok.decaps_sep_arcsec.median():.2f}"  (min {ok.decaps_sep_arcsec.min():.2f}")')
print(f'within sigma_pos: {int((ok.decaps_sep_arcsec<=ok.sigma_pos_arcsec).sum())}/61')
print(f'local density: median {ok.decaps_density_per_as2.median():.3f}/arcsec^2 '
      f'(~{ok.decaps_density_per_as2.median()*3600:.0f}/arcmin^2)')
print(f'chance-alignment P within sigma_pos: median {ok.decaps_p_chance.median():.3f}')
print(f'saved -> {OUT}')
