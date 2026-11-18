# Repository Guidelines

## Project Structure & Module Organization
- Source package: `pyfosc/` (core modules: `fosc.py`, `extraction.py`, `masters.py`, utilities in `utils.py`).
- Pipeline scripts: `src/` (callable steps like `makezero_ccdp.py`, `wavecal2m.py`, etc.).
- Data assets: `database/`, `extinction/`, `iraf_data/` (reference FITS/tables). Do not commit large new data unless essential.
- Config: `config/` (templates such as `instrconfig.json.temp`).
- CLI helpers: `pyfosc_init` (setup) and `pyfosc_run.sh` (batch pipeline).

## Build, Test, and Development Commands
- Create env: `python -m venv .venv && source .venv/bin/activate && pip install -U pip`.
- Install in editable mode: `pip install -e .` or with GUI extras: `pip install -e .[gui]`.
- CLI (installed): `pyfosc init` (creates `raw/`, `data/`, `myfosc.json`) and `pyfosc run` (executes standard pipeline).
- Legacy (repo): `bash ./pyfosc_init` and `bash ./pyfosc_run.sh` remain usable from the repo root.
- Build distribution (optional): `python setup.py sdist bdist_wheel`.

## Command-Line Usage
- `pyfosc init --telescope LJT --slit slit2.5 --grism G3 --copy-raw`
- `pyfosc run` (pre-wavecal, then wavelength/flux/telluric calibration; no IRAF needed).
- `pyfosc run --std-file wacrf....fits --star feige34 --overwrite` to choose the standard and redo outputs.
- `pyfosc notebook --name prewavecal` (or `postcal`) to scaffold a runnable notebook.

## Python API (Notebook-Friendly)
- Import: `from pyfosc.pipeline import PreWaveCal`.
- Typical flow:
  - `pw = PreWaveCal('.')`; `pw.discover()`; `pw.build_master_bias()`; `pw.calibrate_bias_and_trim()`;
  - `pw.build_master_flat()`; `pw.normalize_flat()`; `pw.apply_flat_correction()`;
  - `pw.cosmic_ray_clean()`; `pw.extract_1d()`.
  - Optional: `pw.extract_1d(guess=<row_index>)` to seed trace.
- Post-extraction: `from pyfosc.postcal import PostExtraction`; `pe = PostExtraction('.')`;
  `pe.identify_lamps()`; `pe.apply_wavelength()`; `pe.build_sensfunc()`; `pe.flux_calibrate()`; `pe.telluric_correct()`.

## Important Session Notes
- Idempotent pipeline: re-runs skip already-processed files using filename patterns and FITS headers (`SUBTRACT_BIAS`, `TRIMSEC`, `FLAT_CORRECT`).
- Output locations: intermediates under `data/` (bias+trim as `<name>_bc`, flat-correct as `f<name>`, CR-clean as `crf<name>`); master frames under `MasterFrames/`.
- File resolution: avoid double prefixes (e.g., `crfcrf…`) by resolving existing filenames on disk before processing.
- Uncertainty handling: CR-cleaned images are saved with a unitless `StdDevUncertainty`; extraction ensures uncertainty exists if missing.
- Gain/readnoise: prefer values from FITS headers (aliases: `GAIN`, `RDNOISE`/`RON`), falling back to instrument defaults only if absent.
- Extraction inputs: use `FOSCFileCollection` (not `ImageFileCollection`) to enable `check_groups`/`set_parameters`.
- list_backup_data: module-complete; renames, backs up to `raw/`, copies to `data/`, writes `imglist.csv`.
- Post-extraction modules: `specio.py` (IRAF-style 1D/4-band files), `wavecal.py` (template matching, `WaveSolution`, `dispcor`), `fluxcal.py`, `telluric.py`, orchestrated by `postcal.PostExtraction`. Pixel coordinates are 0-based in code; IRAF `idlamp_*` records are 1-based.
- Lamp lines are flat-topped (slit images); centre them with `wavecal.recenter_lines` (half-maximum centroid), not a Gaussian fit.
- XLT Fe-Ar line wavelengths come from the BFOSC atlas list (`iraf_data/newlinelists/fear_bfosc_atlas.txt`), not from the `idlamp_*` records, which disagree with each other and with sky lines by several Angstrom between 4900 and 6700 A.
- Flux zero point follows IRAF `standard` (`fluxcal.FNU_ZERO = 3.68e-20`); telluric fits must stay guarded (airmass ratio is the fallback).
- IRAF is legacy only: `src/*.py` and `pyfosc/steps/{ccdotz,makereflat2m,doapall}.py` still import pyraf; do not add new pyraf imports.
- No new third-party dependencies for calibration code: numpy, scipy, astropy, matplotlib only.
- CLI fixes: `pyfosc-init` JSON building via `json.dumps`; `pyfosc notebook` now accepts `--name` with default.

## Coding Style & Naming Conventions
- Python 3.9+; follow PEP 8 with 4‑space indents and 88–100 char lines.
- Names: modules and functions `snake_case`; classes `CapWords`; constants `UPPER_SNAKE_CASE`.
- Prefer explicit imports; avoid wildcard imports. Add lightweight docstrings for public functions/classes.

## Testing Guidelines
- Run `pytest -q`; tests live under `tests/` as `test_*.py` and use synthetic spectra plus the packaged templates.
- Keep tests fast and data‑light; use small samples or fixtures rather than full FITS.
- Example: `pip install pytest` then `pytest -q`.

## Commit & Pull Request Guidelines
- Commit messages: imperative, concise, and scoped (e.g., "use BoxcarExtract for bright sources").
- PRs must include: purpose/summary, minimal reproducible steps, and notes on data requirements. Link issues when applicable.
- Avoid committing generated outputs and large binaries; respect `.gitignore`.

## Security & Configuration Tips
- Instrument config is rendered to `myfosc.json` from `config/instrconfig.json.temp` by `pyfosc_init`; review before running.
- Do not store credentials in the repo. Keep local paths/config outside version control.

## Agent-Specific Instructions
- When editing, keep changes minimal and focused; do not rename entry scripts or move data directories.
- Prefer adding new functions in `pyfosc/` and orchestrate via `src/` scripts to preserve current workflow.
