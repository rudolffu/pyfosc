# Installation

PyFOSC is installed from its Git repository. It needs Python 3.9 or later and no IRAF.

## Before you start

- A working Python (conda, miniforge or the system Python).
- `git`.

## Install

1. Create and activate an environment (recommended, so PyFOSC does not disturb other work):

   ```bash
   python -m venv .venv
   ```

   ```bash
   source .venv/bin/activate
   ```

2. Get the code and install it in editable mode:

   ```bash
   git clone https://github.com/rudolffu/pyfosc.git
   ```

   ```bash
   cd pyfosc
   ```

   ```bash
   python -m pip install -e .
   ```

   "Editable" means the `pyfosc` command uses the files in this folder, so `git pull`
   updates your installation. Keep the folder: the reference lamps, standard-star tables
   and extinction curves (`database/`, `iraf_data/`, `extinction/`) are read from it.

3. Check that it worked:

   ```bash
   pyfosc --help
   ```

   You should see the three sub-commands `init`, `run` and `notebook`.

## What gets installed

`numpy`, `scipy`, `matplotlib`, `pandas`, `astropy`, `ccdproc`, `specutils`, `specreduce`,
`colorama` and `packaging`. Nothing else is required.

## Optional extras

- **Notebooks.** `pyfosc notebook` writes a `.ipynb` file; to open it you need Jupyter
  (`python -m pip install notebook`).
- **Qt wavelength tool.** `python -m pip install -e .[gui]` adds `pyqtgraph` and `pyqt5`
  for the older Qt tool in `pyfosc/qtwavecal.py`. The standard pipeline does not use it:
  manual line identification opens an ordinary matplotlib window.
- **Tests.** `python -m pip install pytest`, then `pytest -q` from the repository folder.

## Legacy IRAF scripts

The scripts under `src/` and the modules `pyfosc.steps.ccdotz`, `makereflat2m` and
`doapall` are the earlier IRAF-based version of the reduction. They are kept for reference
and still need [IRAF/PyRAF](https://iraf-community.github.io/). Nothing in `pyfosc run`
uses them, so you can ignore them unless you want to reproduce an old reduction exactly.

## Next

Reduce a night with the [tutorial](tutorial.md).
