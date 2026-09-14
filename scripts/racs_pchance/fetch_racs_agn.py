#!/usr/bin/env python
"""Pull RACS-mid1 (Duchesne+2024, VizieR J/other/PASA/41.3/sourcesm) catalog
rows for the AGN-classified sources in data/gp_long_transient_match.xlsx.

Keeps the official catalog columns only: RAJ2000/DEJ2000 and their formal
uncertainties e_RAJ2000/e_DEJ2000 (degrees), plus Ftot / SCode / Flag.
Nearest RACS source within SEARCH_ARCSEC of the SPT centroid wins, matching
the 30" convention already used in outputs/v4_crossmatch_table.csv.

Writes outputs/racs_agn_positions.csv
"""
import io
import os
import time

import numpy as np
import pandas as pd
import requests

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
XLSX = os.path.join(REPO, "data", "gp_long_transient_match.xlsx")
URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
TABLE = '"J/other/PASA/41.3/sourcesm"'
SEARCH_ARCSEC = 30.0
FLOOR_RA_RACS = 0.65   # arcsec, Duchesne+2024 Table 5 (vs ICRF3)
FLOOR_DEC_RACS = 0.97  # arcsec, Duchesne+2024 Table 5 (corrected dec)
COLS = ["RACS-MID1", "RAJ2000", "DEJ2000", "DEdegc",
        "e_RAJ2000", "e_DEJ2000", "e_DEdegc",
        "Ftot", "e_Ftot", "Fpeak", "Noise", "Maj", "Min", "SCode", "NGauss",
        "Flag"]


def select_agn(df):
    """AGN-marked rows (classification contains 'agn'), plus Blazars flagged."""
    cl = df["classification"].fillna("").astype(str)
    is_agn = cl.str.contains("agn", case=False)
    is_blazar = cl.str.contains("blazar", case=False)
    out = df[is_agn | is_blazar].copy()
    out["agn_group"] = np.where(is_blazar[is_agn | is_blazar], "blazar", "agn")
    out = out.sort_values("agn_group").drop_duplicates("name", keep="first")
    return out.sort_values("name").reset_index(drop=True)


def cone(ra, dec, radius_arcsec):
    q = (f"SELECT {','.join(chr(34)+c+chr(34) for c in COLS)} FROM {TABLE} "
         f"WHERE 1=CONTAINS(POINT('ICRS',RAJ2000,DEJ2000),"
         f"CIRCLE('ICRS',{ra:.6f},{dec:.6f},{radius_arcsec/3600.0:.6f}))")
    for attempt in (1, 2, 3, 4):
        try:
            r = requests.post(URL, data=dict(REQUEST="doQuery", LANG="ADQL",
                                             FORMAT="csv", QUERY=q),
                              timeout=(30, 180))
            r.raise_for_status()
            return pd.read_csv(io.StringIO(r.text))
        except Exception as e:
            print(f"  retry {attempt}: {type(e).__name__}: {e}", flush=True)
            time.sleep(5 * attempt)
    raise RuntimeError("VizieR query failed")


def sep_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = map(np.deg2rad, (ra1, dec1, ra2, dec2))
    num = np.hypot(np.cos(d2) * np.sin(r2 - r1),
                   np.cos(d1) * np.sin(d2) - np.sin(d1) * np.cos(d2) * np.cos(r2 - r1))
    den = np.sin(d1) * np.sin(d2) + np.cos(d1) * np.cos(d2) * np.cos(r2 - r1)
    return np.rad2deg(np.arctan2(num, den)) * 3600.0


def main():
    src = select_agn(pd.read_excel(XLSX))
    print(f"{len(src)} AGN-marked sources "
          f"({(src.agn_group=='agn').sum()} agn / "
          f"{(src.agn_group=='blazar').sum()} blazar)", flush=True)

    rows = []
    for _, s in src.iterrows():
        c = cone(s["RA"], s["Dec"], SEARCH_ARCSEC)
        rec = {"name": s["name"], "classification": s["classification"],
               "agn_group": s["agn_group"], "spt_ra": s["RA"],
               "spt_dec": s["Dec"], "l": s["l"], "b": s["b"],
               "n_racs_in_30as": len(c)}
        if len(c):
            c = c.copy()
            c["sep"] = sep_arcsec(s["RA"], s["Dec"], c["RAJ2000"], c["DEJ2000"])
            m = c.loc[c["sep"].idxmin()]
            rec.update({f"racs_{k}": m[k] for k in COLS})
            rec["racs_sep_arcsec"] = m["sep"]
        rows.append(rec)
        print(f"  {s['name']}: {len(c)} RACS in 30\""
              + (f", nearest {rec.get('racs_sep_arcsec', float('nan')):.2f}\""
                 f" e_RA={3600*rec['racs_e_RAJ2000']:.2f}\""
                 f" e_Dec={3600*rec['racs_e_DEdegc']:.2f}\"" if len(c) else ""),
              flush=True)
        time.sleep(0.2)

    out = pd.DataFrame(rows)
    # published RACS-mid astrometric floor: Duchesne+2024 (PASA 41, e003)
    # Table 5, primary catalogue vs ICRF3: da.cos(d) = -0.14 +/- 0.65",
    # dd_corr = -0.03 +/- 0.97".  Used exactly like the Wan+2025 SPT pointing
    # floors in the existing Gaia pipeline: sigma = sqrt(floor^2 + fit^2).
    out["e_ra_fit_arcsec"] = out["racs_e_RAJ2000"] * 3600.0
    out["e_dec_fit_arcsec"] = out["racs_e_DEdegc"] * 3600.0
    out["sigma_ra_arcsec"] = np.hypot(out["e_ra_fit_arcsec"], FLOOR_RA_RACS)
    out["sigma_dec_arcsec"] = np.hypot(out["e_dec_fit_arcsec"], FLOOR_DEC_RACS)
    out["racs_dec_used"] = out["racs_DEdegc"]
    path = os.path.join(REPO, "outputs", "racs_agn_positions.csv")
    out.to_csv(path, index=False)
    print("saved", path)
    print(f"matched {out['racs_RAJ2000'].notna().sum()}/{len(out)}")


if __name__ == "__main__":
    main()
