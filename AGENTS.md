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
- Run initializer: `bash ./pyfosc_init` (creates `raw/`, `data/`, and `myfosc.json`).
- Run pipeline: `bash ./pyfosc_run.sh` (executes scripts in `src/` in order). Run individual steps with `python src/<step>.py`.
- Build distribution (optional): `python setup.py sdist bdist_wheel`.

## Coding Style & Naming Conventions
- Python 3.9+; follow PEP 8 with 4‑space indents and 88–100 char lines.
- Names: modules and functions `snake_case`; classes `CapWords`; constants `UPPER_SNAKE_CASE`.
- Prefer explicit imports; avoid wildcard imports. Add lightweight docstrings for public functions/classes.

## Testing Guidelines
- No test suite exists yet. Use `pytest` with files under `tests/` named `test_*.py`.
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
