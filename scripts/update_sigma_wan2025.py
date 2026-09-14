#!/usr/bin/env python
"""Switch sigma_pos from Tandoi+2024 (winter field) to Wan+2025 (our field).

Wan et al. 2025 (arXiv:2509.08962) Sec 4.4, Eqs 4-6 — calibrated ON the
SPT-3G Galactic plane field:
    sigma_RA^2  = 3.05^2 + (50/sqrt(TS))^2 + sigma_rel_RA^2
    sigma_Dec^2 = 4.42^2 + (49/sqrt(TS))^2 + sigma_rel_Dec^2
sigma_sys,abs (3.05/4.42) is the AT20G-vs-average-map offset scatter for
this field, replacing Tandoi's 3.73/4.35 measured on the winter field.

sigma_sys,rel is NOT included yet: Wan's definition (per-obs maps during
the flare vs average map) applies to short flares; our centroids come from
yearly coadds, where per-obs pointing scatter averages down ~1/sqrt(N_obs).
The appropriate analog — yearly-map vs average-map registration scatter,
estimable per obs_max year from data/coadd 3yr yearly catalogs — is pending
confirmation with Yujie. Hook: SIGMA_REL_RA/DEC in src/positions.py, keyed by year.

Offline update (no re-query needed): the new sigma is strictly smaller
(new^2 = old^2 - 3.73^2 - 4.35^2 + 3.05^2 + 4.42^2 = old^2 - 3.997), so the
CSC2 3*sigma search radius only shrinks — existing matches are re-judged
against the recorded CSC2_sep_arcsec; no new match can appear.

Only CSC2 depends on sigma (search radius = 3*sigma per source). All other
catalogs use fixed 30" radius and a fixed 10" strict/loose threshold.

Output: outputs/v4_crossmatch_table.csv updated in place;
        previous version saved to
        outputs/archive/v4_crossmatch_table_backup_tandoi_sigma.csv
"""
import os
import sys

import pandas as pd

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "src"))

from positions import sigma_pos_wan2025   # noqa: E402

TABLE  = os.path.join(REPO, "outputs", "v4_crossmatch_table.csv")
BACKUP = os.path.join(REPO, "outputs", "archive",
                      "v4_crossmatch_table_backup_tandoi_sigma.csv")


def main():
    df = pd.read_csv(TABLE)
    print(f"Loaded {len(df)} sources")
    if not os.path.isfile(BACKUP):
        df.to_csv(BACKUP, index=False)
        print(f"Backup written: {os.path.basename(BACKUP)}")
    else:
        print(f"Backup already exists, not overwriting: {os.path.basename(BACKUP)}")

    old_sigma = df["sigma_pos_arcsec_heuristic"].copy()
    df["sigma_pos_arcsec_heuristic"] = [
        sigma_pos_wan2025(snr, yr) for snr, yr in zip(df["snr_max"], df["obs_max"])
    ]
    df.drop(columns=["sigma_pos_method"], errors="ignore", inplace=True)

    d = df["sigma_pos_arcsec_heuristic"] - old_sigma
    print(f"\nsigma_pos: old median {old_sigma.median():.2f}\" -> "
          f"new median {df['sigma_pos_arcsec_heuristic'].median():.2f}\" "
          f"(shift {d.min():.2f}\" .. {d.max():.2f}\")")

    # ── Re-judge CSC2 matches against the shrunken 3*sigma radius ──
    new_r = 3.0 * df["sigma_pos_arcsec_heuristic"]
    df["CSC2_search_radius_arcsec"] = new_r.round(2)

    was_match = df["CSC2_match"].fillna(False).astype(bool)
    sep       = pd.to_numeric(df["CSC2_sep_arcsec"], errors="coerce")
    dropped   = was_match & (sep > new_r)

    df.loc[dropped, "CSC2_match"]        = False
    df.loc[dropped, "CSC2_match_status"] = "outside_3sigma_wan2025"
    # keep sep/cat_ra/cat_dec/flux for the record; only the flag flips
    df["has_csc2_match_within_3sigma"] = df["CSC2_match"].fillna(False).astype(bool)

    print(f"\nCSC2: {int(was_match.sum())} matches under Tandoi sigma, "
          f"{int(dropped.sum())} dropped under Wan sigma")
    if dropped.any():
        print(df.loc[dropped, ["id", "snr_max", "CSC2_sep_arcsec",
                               "CSC2_search_radius_arcsec"]].to_string(index=False))

    df.to_csv(TABLE, index=False)
    print(f"\nSaved -> {os.path.basename(TABLE)}")


if __name__ == "__main__":
    main()
