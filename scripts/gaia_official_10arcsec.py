#!/usr/bin/env python
"""Gaia crossmatch using SPT-3G official get_gaia_prob() method.

For each of the 89 v4-merged sources:
  1. Query Gaia DR3 in a 1° cone around the SPT centroid (via astroquery)
  2. Feed the results into spt3g.sources.gaia.get_gaia_prob(return_all=True)
  3. Filter kept matches to those within 10" of SPT centroid
  4. Compute the SPT-3G association-TS-based p-value for every kept match

TS convention: TS = snr_max**2. For matched-filter Gaussian-noise point-source
detection, likelihood-ratio theory gives TS = SNR^2 exactly; this is the same
identity that reduces the Wan+2025 Eq.6 form theta_beam/sqrt(TS) to
theta_beam/|SNR| we've been using throughout.

Output: outputs/gaia_official_10arcsec.csv
Row = one SPT source × one Gaia candidate within 10".
"""
import os, sys, warnings
warnings.filterwarnings("ignore")

# Point at the arm64-built spt3g install
sys.path.insert(0, os.path.expanduser("~/spt3g_software/build"))

import numpy as np
import pandas as pd
from spt3g import core
from spt3g.sources.gaia import get_gaia_prob
from astroquery.gaia import Gaia

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE  = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
OUT    = os.path.join(REPO, "outputs", "gaia_official_10arcsec.csv")

SEARCH_RADIUS_ARCSEC = 10.0
# Background cone. get_gaia_association's own bg_radius is 1 deg, but its
# association_TS depends on the background only through the brightness
# QUANTILE log(num_total/num_brighter), which is sample-size invariant.
# 0.2 deg keeps the quantile estimate stable while cutting the per-source
# row count ~25x (Galactic-bulge cones at 1 deg exceed 2M rows).
BG_CONE_DEG          = 0.12   # (first 60 sources ran at 0.2; quantile-invariant)


def query_gaia_local(ra_deg, dec_deg, cone_deg=BG_CONE_DEG):
    """ADQL cone query — get every Gaia DR3 source in the cone.

    Uses the async endpoint with an explicit ORDER BY sep_as, so the nearest
    matches survive any truncation, wherever it might come from.

    Checked against the archive 2026-08-24 (do not repeat the old claim of a
    "2000-row sync cap" — it is not what the service does):
      - the TAP capabilities declare outputLimit default = hard = 3,000,000
        rows (https://gea.esac.esa.int/tap-server/tap/capabilities)
      - an anonymous /sync query returns 50,000 rows untruncated
      - COUNT(*) in our 0.12 deg cones matches the local catalog exactly for
        the sparsest / median / densest source (1,232 / 19,301 / 52,431)
    Client-side truncation is still a real trap: astroquery applies its own
    Gaia.ROW_LIMIT to the cone/object helpers, so use launch_job_async with a
    raw ADQL query as done here.
    """
    q = (
        "SELECT source_id, ra, dec, phot_g_mean_mag, parallax, "
        "parallax_over_error, pmra, pmra_error, pmdec, pmdec_error, "
        f"DISTANCE(POINT('ICRS', ra, dec), POINT('ICRS', {ra_deg}, {dec_deg})) "
        "*3600 AS sep_as "
        "FROM gaiadr3.gaia_source "
        f"WHERE CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {ra_deg}, {dec_deg}, {cone_deg})) = 1 "
        "AND phot_g_mean_mag <= 20"    # matches official query_gaia max_g_mag=20
        # NOTE: no ORDER BY — sorting >1M rows server-side is what 500'd the
        # Galactic-bulge-center queries. We filter by sep ourselves, and with
        # no TOP cap an unordered result loses nothing.
    )
    job = Gaia.launch_job_async(q, dump_to_file=False)
    return job.get_results()


def query_gaia_vizier_fallback(ra_deg, dec_deg, cone_deg=BG_CONE_DEG):
    """Fallback for positions ESA's TAP persistently 500s on: query the
    official Gaia DR3 mirror at VizieR (I/355/gaiadr3) and rename columns
    to match the ESA result schema. Same underlying catalog."""
    from astroquery.vizier import Vizier
    from astropy.coordinates import SkyCoord
    from astropy.table import Table
    import astropy.units as u

    v = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag",
                        "Plx", "RPlx", "pmRA", "e_pmRA", "pmDE", "e_pmDE"],
               column_filters={"Gmag": "<=20"}, row_limit=-1)
    res = v.query_region(SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg),
                         radius=cone_deg*u.deg, catalog="I/355/gaiadr3")
    if not res:
        return Table(names=["source_id", "ra", "dec", "phot_g_mean_mag",
                            "parallax", "parallax_over_error",
                            "pmra", "pmra_error", "pmdec", "pmdec_error",
                            "sep_as"])
    t = res[0]
    out = Table()
    out["source_id"]           = np.asarray(t["Source"], dtype=np.int64)
    out["ra"]                  = np.asarray(t["RA_ICRS"], dtype=float)
    out["dec"]                 = np.asarray(t["DE_ICRS"], dtype=float)
    out["phot_g_mean_mag"]     = np.asarray(t["Gmag"], dtype=float)
    out["parallax"]            = np.asarray(t["Plx"], dtype=float)
    out["parallax_over_error"] = np.asarray(t["RPlx"], dtype=float)
    out["pmra"]                = np.asarray(t["pmRA"], dtype=float)
    out["pmra_error"]          = np.asarray(t["e_pmRA"], dtype=float)
    out["pmdec"]               = np.asarray(t["pmDE"], dtype=float)
    out["pmdec_error"]         = np.asarray(t["e_pmDE"], dtype=float)
    c0 = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
    cs = SkyCoord(ra=out["ra"]*u.deg, dec=out["dec"]*u.deg)
    out["sep_as"] = c0.separation(cs).arcsec
    return out


def to_spt3g_stars_frame(t):
    """Convert an astroquery Gaia result (deg) into the pandas DF format that
    get_gaia_prob expects: `ra`, `dec` in G3Units (radians internally),
    `g_mag`, plus source_id."""
    if len(t) == 0:
        return pd.DataFrame(columns=["source_id", "ra", "dec", "g_mag"])
    df = pd.DataFrame({
        "source_id": np.asarray(t["source_id"], dtype=np.int64),
        "ra":        np.asarray(t["ra"],  dtype=float) * core.G3Units.deg,
        "dec":       np.asarray(t["dec"], dtype=float) * core.G3Units.deg,
        "g_mag":     np.asarray(t["phot_g_mean_mag"], dtype=float),
    })
    df.index = df["source_id"].values
    return df


DONE_FILE = os.path.join(REPO, "outputs", "gaia_official_done.txt")


def main():
    import time
    v4 = pd.read_csv(TABLE)
    print(f"Loaded {len(v4)} SPT sources from {os.path.basename(TABLE)}")

    # Resume support: skip sources already queried successfully; append rows
    # incrementally so a crash never loses completed work.
    done = set()
    if os.path.isfile(DONE_FILE):
        done = set(ln.strip() for ln in open(DONE_FILE) if ln.strip())
        print(f"Resume: {len(done)} sources already done, skipping them")

    rows = []
    for i, r in v4.iterrows():
        sid       = r["id"]
        if sid in done:
            continue
        ra, dec   = float(r["ra_deg"]), float(r["dec_deg"])
        snr       = float(r["snr_max"])
        ts        = snr * snr                 # TS = SNR^2 proxy (see docstring)

        # 1) Query Gaia — exponential backoff; ESA TAP 500s under bursty load
        gaia_tbl = None
        for attempt in (1, 2, 3, 4):
            try:
                gaia_tbl = query_gaia_local(ra, dec)
                break
            except Exception as e:
                print(f"[{i+1:2d}/{len(v4)}] {sid}: Gaia query ERROR "
                      f"(attempt {attempt}) {type(e).__name__}: {e}")
                if attempt < 4:
                    time.sleep(20 * attempt)
        if gaia_tbl is None:
            # ESA exhausted → VizieR mirror of the same catalog
            try:
                gaia_tbl = query_gaia_vizier_fallback(ra, dec)
                print(f"[{i+1:2d}/{len(v4)}] {sid}: VizieR fallback OK "
                      f"({len(gaia_tbl)} rows)")
            except Exception as e:
                print(f"[{i+1:2d}/{len(v4)}] {sid}: VizieR fallback ERROR "
                      f"{type(e).__name__}: {e}")
                continue

        stars_df = to_spt3g_stars_frame(gaia_tbl)
        n_in_cone = len(stars_df)
        if n_in_cone == 0:
            print(f"[{i+1:2d}/{len(v4)}] {sid}: 0 Gaia sources in 1° cone")
            continue

        # 2) Call spt3g official function with return_all=True
        # Pointing floors: Wan+2025 (arXiv:2509.08962) Sec 4.4 sigma_sys,abs,
        # calibrated on THIS Galactic-plane field (AT20G vs average map) —
        # overriding the function defaults 3.73/4.35 (Tandoi+2024, winter
        # field). Per-event sigma_sys,rel pending (Wan's flare-window
        # definition doesn't transfer to yr-long variables); fold it into
        # these floors in quadrature once estimated.
        try:
            all_matches = get_gaia_prob(
                stars_df,
                event_ra=ra * core.G3Units.deg,
                event_dec=dec * core.G3Units.deg,
                event_ts=ts,
                pointing_floor_ra=3.05 * core.G3Units.arcsec,
                pointing_floor_dec=4.42 * core.G3Units.arcsec,
                return_all=True,
            )
        except Exception as e:
            print(f"[{i+1:2d}/{len(v4)}] {sid}: get_gaia_prob ERROR "
                  f"{type(e).__name__}: {e}")
            continue

        # 3) Keep candidates within SEARCH_RADIUS_ARCSEC
        if len(all_matches) == 0:
            print(f"[{i+1:2d}/{len(v4)}] {sid}: 0 candidates from spt3g")
            continue
        sep_as = all_matches["dist"] / core.G3Units.arcsec
        keep   = all_matches[sep_as < SEARCH_RADIUS_ARCSEC]
        print(f"[{i+1:2d}/{len(v4)}] {sid}: cone={n_in_cone}, in-10\"={len(keep)}")

        # 4) Persist per-match rows (may be 0, 1, or many)
        src_rows = []
        for gaia_id, m in keep.iterrows():
            # source_id from the INDEX: it stays int64. The row Series `m` is
            # float64-upcast by iterrows(), which rounds 19-digit Gaia ids.
            source_id_val = int(gaia_id)
            # Pull Gaia raw astrometry from the original astroquery table
            hit = gaia_tbl[gaia_tbl["source_id"] == source_id_val]
            if len(hit) == 0:
                continue
            hit = hit[0]

            # Column names follow Gaia DR3 natively where a native name exists.
            src_rows.append({
                "spt_id":              sid,
                "p_chance":            float(m["prob"]),
                "spt_ra_deg":          ra,
                "spt_dec_deg":         dec,
                "spt_snr_max":         snr,
                "source_id":           source_id_val,
                "ra":                  float(m["ra"] / core.G3Units.deg),
                "dec":                 float(m["dec"] / core.G3Units.deg),
                "gaia_sep_arcsec":     float(m["dist"] / core.G3Units.arcsec),
                "phot_g_mean_mag":     float(m["g_mag"]),
                "association_TS":      float(m["association_TS"]),
                "parallax":            float(hit["parallax"])            if hit["parallax"] is not None            else np.nan,
                "parallax_over_error": float(hit["parallax_over_error"]) if hit["parallax_over_error"] is not None else np.nan,
                "pmra":                float(hit["pmra"])                if hit["pmra"] is not None                else np.nan,
                "pmra_error":          float(hit["pmra_error"])          if hit["pmra_error"] is not None          else np.nan,
                "pmdec":               float(hit["pmdec"])               if hit["pmdec"] is not None               else np.nan,
                "pmdec_error":         float(hit["pmdec_error"])         if hit["pmdec_error"] is not None         else np.nan,
            })

        # Incremental flush: append this source's rows + mark it done, so an
        # interrupted run never loses or repeats completed sources.
        if src_rows:
            header = not (os.path.isfile(OUT) and os.path.getsize(OUT) > 10)
            pd.DataFrame(src_rows).to_csv(OUT, mode="a", index=False, header=header)
        with open(DONE_FILE, "a") as f:
            f.write(sid + "\n")
        rows.extend(src_rows)

    out = pd.read_csv(OUT) if os.path.isfile(OUT) and os.path.getsize(OUT) > 10 \
        else pd.DataFrame()
    print(f"\nTotal in {os.path.basename(OUT)}: {len(out)} match rows")
    if len(out):
        print(f"Unique SPT sources with ≥1 Gaia within 10\": {out['spt_id'].nunique()}")


if __name__ == "__main__":
    main()
