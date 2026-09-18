# SPT-3G Galactic Plane — counterpart identification

Multi-wavelength counterpart identification for long-timescale (~1 yr) variable
sources detected in the SPT-3G Galactic plane survey. The question this repo
answers per source is: **is there a real counterpart, and is it extragalactic
(AGN) or Galactic (star / nova / symbiotic)?**

Scope note: this is counterpart ID, not transient discovery and not a
population analysis.

---

## Data layout

**`data/` and `outputs/` are git-ignored** — collaborators are expected to hold
identical local copies. Paths below are what the code hard-codes, so keep the
names byte-identical.

```
data/
  centroids 4_yr.cat                 THE source list — 73 sources, 4-yr detection
  centroids.cat                      v3 list (32 dual-band sources), older figures
  centroids_merged.cat               deprecated 89-row 3yr+4yr merge, do not use
  gp_long_transient_match.xlsx       per-source classification (AGN / nova / artifact / …)

  coadd/
    galaxy_3yr_pc_gc_v2_90.fits      3-yr coadd, 90 GHz  (HDU0 = T in mK_CMB, HDU1 = weight)
    galaxy_3yr_pc_gc_v2_150.fits     3-yr coadd, 150 GHz — also defines the footprint
  coadd 3yr/                         per-year source catalogs (2023/2024 × 90/150 GHz)
  v3/, v4/                           per-source SPT cutouts + per-coadd catalogs

  gaia_null/
    gaia_galactic_g20.npz            Gaia DR3 G<=20 over the field — 162 065 280 stars
    ts_gal_ts{TS}_n250000.npz        null association-TS samples, one per event-TS
    Gaia_association_to_pval_mappings_galactic.g3   TS-dependent mapping, official g3 format
    tiles_galactic/                  raw VizieR download tiles (resumable)

  wise_null/
    allwise_field.npz                AllWISE over the field — 5 948 365 sources (W1)
    wise_null_n400000_seed9001.npz   null association-TS, 400k positions × 26 sigma pairs
    tiles/                           raw IRSA TAP download strips (resumable)
  spt3g_cluster/                     files copied off the SPT-3G cluster

  fits_cache/                        figure-panel FITS cache (auto-filled, safe to delete)
  unwise/                            unWISE neo6 cutout cache
  milliquas_v8.csv                   Milliquas v8 quasar catalog (99 MB)
```

The Gaia catalog and the null samples are large and slow to rebuild but fully
reproducible — see "Rebuilding derived data" below.

```
outputs/
  v4_crossmatch_table.csv            THE master per-source table (73 rows)
  gaia_official_10arcsec.csv         official get_gaia_prob run, one row per SPT × Gaia pair
  gaia_pchance_galactic *.csv/.xlsx  re-scored with this field's null, one pair per floor
  gaia_pchance_floors.csv            all floors side by side  (+ floors_source_ts.csv)
  racs_agn_pchance*.csv/.xlsx        AGN sources re-scored at their RACS-mid position
  wise_pchance.csv/.xlsx             AllWISE counterparts + p_chance  (+ wise_vs_gaia_pchance.csv)
  3yr_classify.csv, crossmatch_results_v3.csv, present_table.csv   3-yr / v3 era tables
  v4_*_summary.csv, spectral_index.csv, roma_bzcat_footprint.csv

  images/                            every figure (subfolders per figure family)
  mappings/                          published two-column association_TS -> p tables
  logs/                              stdout of the long null / scoring runs
  archive/                           superseded one-off products, kept for provenance
  legacy/, pptx_gen/                 pre-v4 material and the slide generator
```

Nothing in `outputs/` is an input to anything else except
`v4_crossmatch_table.csv` and `gaia_official_10arcsec.csv`, which most scripts
read — do not move or rename those two.

---

## Layout

```
src/            importable modules (scripts/notebooks put src/ on sys.path themselves)
scripts/        runnable analysis + figure scripts
notebooks/      looking at results: tables, a few sources, diagnostics
references/     paper notes (README lists arXiv IDs; the PDFs are git-ignored)
requirements.txt
```

Setup: `pip install -r requirements.txt`. `spt3g_software` is only needed by
`scripts/gaia_official_10arcsec.py` and for writing the `.g3` mapping.

### `src/`

| module | what it does |
|---|---|
| `paths.py` | project-relative paths, so anything runs from any CWD |
| `positions.py` | SPT-3G positional uncertainty, Wan+2025 Eqs 4–6 (`sigma_pos_wan2025`) |
| `crossmatch.py` | nearest-match VizieR cone search, one catalog at one position (`crossmatch_one`) |
| `cutouts.py` | SPT cutout / coadd I/O, reprojection, aperture photometry |
| `fits_cutouts.py` | calibrated-FITS fetchers with disk cache — AllWISE W3/W4, DECaPS r, NVSS, RACS-mid (CASDA) |
| `unwise.py` | unWISE W1/W2 neo6 cutouts |
| `source_figure.py` | **the per-source multi-panel cutout figure** — `plot_source(source, ...)`, SPT contours, `spt_markers` / `csc2_marker` |
| `wide_zoom_figure.py` | 3′ DECaPS r + 30″ DECaPS gri figure (HiPS, morphology only) — `plot_wide_zoom` |
| `plotting.py` | coadd cutout grid |
| `allwise_flags.py` | AllWISE `qph` / `ccf` flag decoding |
| `pchance.py` | null sample → chance probability (`survival_table` / `survival_interp`) |
| `units.py` | DN → physical units, colorbar labels, HiPS survey IDs |
| `galactic.py` | Galactic coordinates and \|b\| binning |

### `scripts/`

`spt.py` is the terminal entry point — alias it once
(`alias spt="$HOME/anaconda3/bin/python <repo>/scripts/spt.py"`) and then
`spt plot <id> --fov 90`, `spt precache`, `spt build-table`, or
`spt run <path>` for anything else, from any directory. `spt list` shows the
commands. Table-rewriting scripts are only reachable through `spt run`.

- **master table** — `build_v4_crossmatch_table.py` builds
  `outputs/v4_crossmatch_table.csv` from `data/centroids 4_yr.cat`; then
  `add_milliquas.py`, `add_decaps.py`, `update_sigma_wan2025.py` add / update
  columns in place (so the build refuses to overwrite without `--overwrite`)
- **other crossmatch** — `check_roma_bzcat.py`, `identify_lum_plot_points.py`
- **Gaia association** — `gaia_official_10arcsec.py` (official `get_gaia_prob`)
- **null calibration** — `null_resample/` (Galactic field)
- **RACS-mid rescoring** — `racs_pchance/` (see below)
- **WISE association** — `wise_pchance/` (AllWISE p_chance, same null machinery)
- **figures** — everything named `plot_*.py`; `plot_source_panels.py` is the
  loop that calls `src.source_figure.plot_source` once per source
- **caching** — `precache_fits.py`, run it before any cutout figure

### `notebooks/`

- `analyze_v4_csc2_candidates.ipynb` — sources with a CSC2 match within 3σ_pos:
  flags, fluxes, cutout panels, wide + zoom optical figures
- `analyze_v4_no_csc2_candidates.ipynb` — the rest: counterpart fingerprints,
  \|b\|, radio spectral index, W1−W2
- `diagnostics/coadd_flux_check.ipynb` — flux level of each source in the 3-yr coadd

Superseded code (winter-field null validation, the home-grown Gaia veto, v3-era
scripts) is kept locally in `outputs/archive/legacy_code/` with a MANIFEST, not
in the repo.

---

## The p_chance machinery

A Gaia counterpart is scored with the SPT-3G association test statistic

```
association_TS = -chi2/2 - ln(num_brighter) + ln(4*pi) + ln(N_cone)
chi2           = (d_RA/sigma_RA)^2 + (d_Dec/sigma_Dec)^2
```

and turned into a chance probability by comparing against a **null built for
this field**: random positions drawn uniformly over the valid coadd footprint,
scored against the same 162M-star Gaia catalog with the same 0.12 deg
background cone. `scripts/null_resample/null_ts.py` is a vectorized
reimplementation of `spt3g.sources.gaia.get_gaia_association`, verified against
the real thing (max |ΔTS| 2e-11; the test lives in the local legacy archive).

Positional sigma always has the same shape — a published survey floor in
quadrature with a per-source statistical term:

| position source | floor (RA / Dec) | statistical term | reference |
|---|---|---|---|
| SPT-3G centroid | 3.05″ / 4.42″ | 50″/√TS, 49″/√TS | Wan+2025 (arXiv:2509.08962) §4.4 |
| RACS-mid | 0.65″ / 0.97″ | catalog `e_RAJ2000` / `e_DEdegc` | Duchesne+2024 (PASA 41, e003) Table 5, vs ICRF3 |

**Two nulls exist and they are not interchangeable.** The winter-field mapping
that ships with `spt3g_software` gives much smaller p-values here — the
Galactic plane is ~2 orders of magnitude denser in Gaia stars. Always use the
Galactic-field mapping for this field.

### Where code goes

`src/` is imported, `scripts/` is executed, notebooks are looked at. One test
decides: **is the output read by other code, or by a person?**

- **`src/`** — any function a second caller could want. If you are about to
  copy a function, it belongs here.
- **`scripts/`** — anything that writes a file other code reads (tables, null
  samples, mappings), runs long, or loops over the whole source list; it
  should be re-runnable with one command.
- **notebooks** — tables, a handful of sources, a parameter to try, a
  diagnostic, with markdown saying what the result means. A notebook calls
  `src/`; it does **not** contain `def`. A definition there cannot be imported,
  so the only way to reuse it is to copy it — which is how `_reproject_spt`
  and `_off_source_rms` once ended up with four copies each.

### Per-source figures

One source's figure is one function call — `src/source_figure.plot_source()`.
It takes the source *record* (id, ra_deg, dec_deg, snr_max, …), not a row index
into a global DataFrame, so the same call works for one source, for a notebook
check, and inside a loop. Running the whole list is just that loop:

```
python scripts/plot_source_panels.py SPT3G_J173508.4-293000.8   # one source
python scripts/plot_source_panels.py --all                      # all 73
python scripts/plot_source_panels.py --all --big DECaPS-r,NVSS  # swap a panel
```

Panels are read from `data/fits_cache/` (run `precache_fits.py` first), so a
full run is offline and takes ~6 s per source. Failures never abort the run:
a panel that cannot be drawn becomes a labelled grey box and is listed in
`outputs/logs/plot_source_panels_errors.csv`.

⚠️ `data/v4/` holds the **3-yr / 61-source** centroid run, so none of the 73
ids match a cutout filename exactly. `find_spt_cutout_id()` falls back to the
nearest cutout within 5″ and the figure title says which file it used and how
far away it was. 41 of 73 resolve; the other 32 get no SPT contours. Dropping
the 4-yr cutouts into `data/v4/` makes all 73 match by id and the note vanishes.

### `scripts/racs_pchance/` — scoring at the radio position

Re-scores the AGN-classified sources at their RACS-mid positions instead of the
SPT centroid, which shrinks sigma from ~5–7″ to ~1″.

```
fetch_racs_agn.py     VizieR J/other/PASA/41.3/sourcesm cone search -> outputs/racs_agn_positions.csv
run_null_racs.py      null for every (sigma_RA, sigma_Dec) pair in one pass  (~12 min for 250k events)
score_racs.py         observed TS + p_chance                                 -> outputs/racs_agn_pchance*.csv
verify_pruning.py     checks the candidate-pruning against a full-cone scan
plot_racs_pchance.py  summary figure
```

Two sigma variants are scored: **A** = catalog fit ⊕ floor (primary,
conservative — the Dec floor partly double-counts the correction-model term
already inside `e_DEdegc`), **B** = published ICRF3 scatter only. They agree on
every conclusion.

⚠️ `null_multi_sigma.py` inherits the 60″/40″ candidate-pruning from
`null_ts.py`. That shortcut assumes roughly isotropic sigma. `verify_pruning.py`
confirms it is bit-exact for the sigma ratios A and B use (~1.8; checked over
2000 footprint positions × 26 sigma pairs, max difference 0), but it **breaks
for strongly anisotropic sigma** — a third variant using catalog fit errors
alone, with sigma_RA down to 0.002″ against sigma_Dec 0.5″, was retired for
exactly this reason plus being unphysical. If you ever add such a variant, run
with `--exact`; the full-cone scan is ~40× slower.

---

## Rebuilding derived data

```bash
# Gaia DR3 G<=20 over the field (hours; resumable, per-tile npz)
python scripts/null_resample/download_gaia_galactic.py
python scripts/null_resample/merge_tiles_galactic.py

# the TS-dependent null mapping family (~80 min, needs a spt3g-software env
# only for writing the .g3 file)
python scripts/null_resample/make_galactic_mappings.py

# figure FITS cache (RACS-mid needs free CASDA/OPAL credentials)
python scripts/precache_fits.py
```

---

## Conventions

- **Never display raw detector units.** Colorbars and axes are μJy / mJy / mag /
  μK. `src/units.py` holds the conversions; unWISE is Vega nanomaggies, not AB.
- **Figure panels come from downloaded calibrated FITS**, not HiPS previews.
  RACS-mid is the one exception (CASDA has no anonymous cutout service), and its
  HiPS was checked against the catalog.
- **No hand-rolled composite statistics.** Use official catalog columns
  (e.g. Gaia `parallax_over_error`) or raw columns.
- **Notebook outputs are stripped before committing** (`nbstripout <file>`).
  Install the pre-commit filter if you would rather not remember:
  `nbstripout --install`.
