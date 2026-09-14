# References

External papers and references this project relies on.

## Available

- **`1205.0811v1.pdf`** — Stern et al. 2012, *Mid-Infrared Selection of AGN with WISE.
  I. Characterizing WISE-Selected AGN in COSMOS* (ApJ, accepted Apr 24 2012).
  Source of the **W1 − W2 ≥ 0.8 (Vega) AGN selection cut** applied in
  `notebooks/analyze_v4_csc2_candidates.ipynb` and
  `notebooks/analyze_v4_no_csc2_candidates.ipynb`.
  ADS: [2012ApJ...753...30S](https://ui.adsabs.harvard.edu/abs/2012ApJ...753...30S)

- **`2509.08962v1.pdf`** — Y. Wan et al. 2025, *Detection of Millimeter-Wavelength
  Flares from Two Accreting White Dwarf Systems in the SPT-3G Galactic Plane
  Survey* (arXiv 2509.08962).
  Core reference for **figure style** (multi-band cutout + SPT contour
  overlay + X-ray contour) and the two SPT-SV reference sources
  **J174417.2-293942** and **J173508.3-292956**. Also the source of the
  positional-uncertainty model in `src/positions.py` (Sec 4.4, Eqs 4–6).
  ADS: [2025arXiv250908962W](https://ui.adsabs.harvard.edu/abs/2025arXiv250908962W)

- **`2511.10785v1.pdf`** — T. Murphy & D. L. Kaplan 2026, *The Dawes Review 13:
  A New Look at The Dynamic Radio Sky* (PASA, 59 pages).
  Comprehensive review of **radio transient classes, image-domain (slow)
  transient surveys, and discovery / classification methods**. Useful as
  the radio-side analogue of what we're trying to do at 95/150 GHz with
  SPT-3G. Murphy is the RACS PI, so RACS-mid context is here too.
  ADS: [2025arXiv251110785M](https://ui.adsabs.harvard.edu/abs/2025arXiv251110785M)

- **`2311.12369_Duchesne2024_RACS-mid_catalogue.pdf`** — S. W. Duchesne et al. 2024,
  *The Rapid ASKAP Continuum Survey V: cataloguing the sky at 1367.5 MHz and the
  second data release of RACS-mid* (PASA 41, e003; arXiv 2311.12369).
  Source of the **RACS-mid positional uncertainties** used in
  `scripts/racs_pchance/`: Table A1 (p. 28) defines `e_RAJ2000` / `e_DEdegc` as
  *fitting* uncertainties only, Section 4.10 (p. 19) documents the ±1–2″
  beam-to-beam astrometric scatter, and Table 5 (p. 21) gives the vs-ICRF3
  offsets **−0.14 ± 0.65″ (RA)** and **−0.03 ± 0.97″ (corrected Dec)** that we
  adopt as the astrometric floor. Also the source of the `DEdegc`
  declination-corrected column convention.
  ADS: [2024PASA...41....3D](https://ui.adsabs.harvard.edu/abs/2024PASA...41....3D)

- **`2401.13525v2_Tandoi2024_flaring_stars.pdf`** — C. Tandoi et al. 2024,
  *Flaring Stars in a Nontargeted Millimeter-wave Survey with SPT-3G* (ApJ 972, 6; arXiv 2401.13525).
  Origin of the Gaia **association_TS and its random-position null** that
  `scripts/null_resample/` rebuilds for the Galactic field; the p_chance
  histograms follow its Fig. 6 layout.

## Wishlist (not yet downloaded)

- Chandra Source Catalog v2.1 release paper (Evans et al. 2024 or release note)
  — for documenting CSC2 thresholds / hardness ratios if we add X-ray
  filtering for the 3-year sources.

- SPT-3G survey description paper — for citing 3G beam, sensitivity, pixel
  scale, matched-filter convention in the eventual writeup.

## Conventions

- Filename = arXiv ID + version, e.g. `1205.0811v1.pdf`, `2509.08962v1.pdf`.
  Keep it traceable.
- **The PDFs themselves are git-ignored** (all public, 37 MB). This README is
  the committed record — add a bullet whenever you drop a new PDF in here.
- One-line bullet per paper above when adding new ones (what it is + what
  we use it for, not just the title).
