# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

MATCHMAKER is a Python package for cross-matching astronomical catalogs (radio/optical surveys) and evaluating chance-alignment probabilities. Built around astropy, pandas, and numpy.

## Setup / Commands

```bash
pip install -e .[dev]        # core + pytest/ruff; add [plot] extra for plotmaker deps (matplotlib, aplpy, ...)
pytest tests/                # run test suite
pytest tests/test_crossmatch.py::TestCrossmatch::test_simple_match_within_separation  # single test
ruff check .                 # lint (config in pyproject.toml; must stay clean, CI enforces)
```

Supported: Python >=3.9, pandas 2 and 3, astropy 5–8. Data files are expected under `../../data/` by default — set the `MATCHMAKER_DATA_PATH` environment variable to override (`DATA_BASE_PATH` in `matchmaker/dataset/__init__.py`).

CI: GitHub Actions (`.github/workflows/ci.yml`) runs ruff + pytest on 3.9/3.11/3.13. `.gitlab-ci.yml` only uploads tagged releases to Zenodo via eossr.

## Architecture

Three top-level packages:

- **`matchmaker/`** — core library.
  - `dataset/__init__.py` defines the central abstractions: `Catalog` (holds a pandas DataFrame in `.df`, column metadata in `.cols`, results in `.matches`), `Column`/`Columns` (map survey-specific column names to labels + astropy units), `Match` (mask, target indices `idx`, separations, chance association), `File`, `Image`.
  - `dataset/*.py` — one `Catalog` subclass per survey (LoTSS, VLASS, NVSS, FIRST, RACS, CLU, Reines, SDSS via astroquery). Each subclass declares its column mappings in `set_survey_specific_columns()` and implements `load_data()`. `lotss.py` is the richest and carries science logic (luminosity/SFR outlier selection, compactness filters, DR2 region boxes).
  - `crossmatch.py` — `crossmatch(obj1, obj2, ...)`: astropy `match_to_catalog_sky` nearest-neighbour match accepted when `sep2d <= sep + n_sigma * sigma` (combined source+target radial positional error), optional ellipse containment filter (`Ellipse`: PA East of North, cos(dec)-scaled tangent plane). `legacy=True` reproduces the historical 3×3 shift-grid method verbatim (quirks documented in `_crossmatch_legacy`) for diffing against old runs. Results are written into both catalogs via `Catalog.add_match()`; nothing is returned.
  - `stats.py` — analytical chance-association probabilities (circular and elliptical).
  - `measures.py`, `fitting.py`, `utils.py` — luminosity distance, power-law scaling, pickle-based state save/load (`states/` folder convention), FITS→DataFrame loaders.
  - `model/` — empirical relations (LOFAR luminosity–SFR power law, WISE colors).
- **`plotmaker/`** — diagnostic plots, image cutout fetching (LoTSS, Pan-STARRS, VO) and plotting (aplpy), spectra.
- **`monte_carlo/`** — `mc.py` runs chance-association Monte Carlo: randomly shifts LoTSS coordinates, re-runs crossmatch, counts spurious matches. Doubles as usage example (`__main__`-guarded script at module bottom).

### Key data-flow convention

Cross-match results live on the Catalog objects: `source.matches[target.name]` gets `mask` (bool over source rows); `target.matches[source.name]` gets `idx` (per-source-row nearest target index), `mask`, `sep2d`. Derived accessors: `mask_idx = np.where(mask)[0]`, `filtered_idx = idx[mask]`. Row access is positional (`df.iloc`), so masks/indices assume DataFrames keep their load-time ordering; filters like `filter_box`/`constrain` reset the index and store the original in an `index_df_*` column.

Units are handled via astropy: every `Column` carries a `.unit`, and `Catalog.prop_to_unit(prop, unit)` converts values on the fly.
