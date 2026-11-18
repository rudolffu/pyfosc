# PyFOSC
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10967240.svg)](https://doi.org/10.5281/zenodo.10967240)  

__Please note: This program is under active development. The documentation may be inconsistent with the latest codes.__

__IRAF/PyRAF is no longer required: all steps of the standard pipeline now run on astropy, ccdproc, specreduce, numpy and scipy.__

PyFOSC is a pipeline toolbox for long-slit spectroscopy data reduction written in Python. It can be used for FOSC (Faint Object Spectrograph and Camera) data from Xinglong/Lijiang 2-meter telescopes.  

BFOSC (Beijing-Faint Object Spectrograph and Camera) is an instrument of the 2.16-m Telescope in Xinglong Observatory, National Astronomical Observatories, Chinese Academy of Sciences (NAOC) (IAU code: 327, coordinates: 40°23′39″ N, 117°34′30″ E). For more information about BFOSC, please see http://www.xinglong-naoc.cn/html/en/gcyq/216/detail-18.html. The details of the Xinglong 2.16-m telescope and BFOSC are also reported in [Fan et al. 2016](https://ui.adsabs.harvard.edu/abs/2016PASP..128k5005F/abstract) and [Zhao et al. 2018](https://ui.adsabs.harvard.edu/abs/2018RAA....18..110Z/abstract).

YFOSC (Yunnan-Faint Object Spectrograph and Camera) is an instrument of the 2.4-m Telescope in Lijiang Observatory, Yunnan Observatories, Chinese Academy of Sciences (YNAO) (IAU code: O44, coordinates: 26°42′33.1″ N, 100°1′51.6″ E). For more information about YFOSC, please see http://wiki.gmg.org.cn/pages/viewpage.action?pageId=557106. The details of the Lijiang 2.4-m telescope and YFOSC are also reported in [Wang et al. 2019](https://ui.adsabs.harvard.edu/abs/2019RAA....19..149W/abstract).

## 1. Installation 

### 1.1. Dependencies
```
numpy, scipy, matplotlib, pandas
astropy
ccdproc
specutils
specreduce
```
These are installed automatically with `pip install -e .`.

### 1.2. IRAF/PyRAF (optional, legacy scripts only)

The scripts under `src/` and the modules `pyfosc.steps.{ccdotz,makereflat2m,doapall}` are the
older IRAF-based versions of the reduction. They are kept for reference and still need
[IRAF/PyRAF](https://iraf-community.github.io/) ([installation notes for Apple Silicon](https://yumingfu.space/tech/2024-iraf-mac-installation/)).
Nothing in `pyfosc run` uses them.

### 1.3. Install PyFOSC

Clone and install in editable mode:
```bash
git clone https://github.com/rudolffu/pyfosc.git
cd pyfosc
python -m pip install -e .
```
No environment variable changes are required. The `pyfosc` command is installed automatically.

## 2. Usage

**New documentation under development!**

### 2.1. Quick start (CLI)

From a working directory containing your FITS files:
```bash
# Initialize project layout and config, and optionally copy/rename raws
pyfosc init --telescope XLT --slit slit2.3 --grism G4 --copy-raw

# Run reductions: pre-wavecal, wavelength, flux and telluric calibration
pyfosc run
```
Notes:
- Raw files are backed up under `raw/` and renamed copies placed in `data/`.
- Intermediate products (bias/trim `_bc`, `f*`, `crf*`) are written under `data/`.
- Master frames are saved under `MasterFrames/` at the repository root.
- Wavelength solutions are saved under `data/database/`, and QA plots under `QAplots/`.

### 2.2. Running the pipeline step by step

#### 2.2.1. Optional: Jupyter workflow

Generate a starter notebook for the pre-wavecal steps:
```bash
pyfosc notebook --name prewavecal --cwd ./
```
Run cells to perform bias/trim, flat/normalization, flat correction, CR cleaning, and 1D extraction.

The basic reduction steps include:
- Bias construction and subtraction, and CCD trimming
- Flat field construction and correction
- Cosmic ray removal
- Tracing and extraction of 1d spectra

#### 2.2.2. Wavelength, flux and telluric calibration

These steps are run by `pyfosc run`. To execute them one by one, from the project directory:

```bash
python -m pyfosc.steps.reidentlamp2m    # lamp solutions from the packaged templates
# python -m pyfosc.steps.identlamp2m    # or: mark lamp lines by hand (matplotlib window)
python -m pyfosc.steps.wavecal2m        # dispersion correction + flux calibration
python -m pyfosc.steps.telluric_base2m  # telluric correction
```

or from Python / a notebook (`pyfosc notebook --name postcal`):

```python
from pyfosc.postcal import PostExtraction
pe = PostExtraction('.')
pe.identify_lamps()
pe.apply_wavelength()
pe.build_sensfunc()      # or build_sensfunc(std_file='wacrf....fits', star='feige34')
pe.flux_calibrate()
pe.telluric_correct()
```

What each step does:

- **Lamp identification.** Each extracted lamp (`af*.fits`) is cross-correlated with the
  reference lamps in `database/` for the same telescope; the lines stored in the matching
  `idlamp_*` record are re-centred and fitted with a sigma-clipped Legendre polynomial.
  The solution is written to `data/database/wavesol_<lamp>.json` and plotted in
  `QAplots/wavesol_<lamp>.pdf`. **Check these plots**: the fit RMS should be a few tenths
  of a pixel. If no template matches (new grism or slit), run `identlamp2m` and mark the
  lines by hand: `m` mark a line and type its wavelength, `d` delete, `f` fit,
  `l` add lines from the line list, `q` accept.
  For XLT Fe-Ar lamps the template only locates the lines: the solution is then refitted
  on the lines of the BFOSC lamp atlas (`iraf_data/newlinelists/fear_bfosc_atlas.txt`),
  which also covers 4900-6700 A where the packaged records have few lines
  (`--no-refine` keeps the template's own lines). A warning is logged when a solution
  has no lines over more than a fifth of the detector.
- **Dispersion correction.** `acrf*` -> `wacrf*`, rebinned (flux conserving) to a linear
  wavelength grid using the lamp nearest in time. The four bands (optimal, raw, sky, sigma)
  are kept.
- **Flux calibration.** The standard star is recognised from its `OBJECT` keyword among the
  stars in `iraf_data/onedstds`; override with `--std-file`, `--star`, `--caldir`
  (the tables of one star differ by a few percent between directories, so pass
  `--caldir` to reproduce an older reduction). The
  sensitivity function is saved as `sens.fits` (`QAplots/sensfunc.pdf`) and applied with
  the site extinction curve: `wacrf*` -> `cwacrf*` in erg/cm2/s/A.
- **Telluric correction.** The transmission in the O2 and H2O bands is measured from the
  standard star and scaled to each target: `cwacrf*` -> `<OBJECT>_<tel><n>_<star>_tel.fits`.
  The scaling is the airmass ratio; a fitted shift and O2/H2O scale replace it only when
  they clearly smooth the spectrum, so faint targets are not fitted to noise
  (`--no-fit` always uses the airmass ratio).

Existing outputs are skipped on re-runs; pass `--overwrite` to redo them.

## 3. Credits

This software uses BSD 3-Clause License.  

Copyright (c) 2019-2025, Yuming Fu  
All rights reserved.  

This software contains sources from third-party softwares.  

The `pyfosc$iraf_data/onedstds` directory is from `IRAF`, and it contains standard calibration data for extinction and sensitivity calibration.  

The `removecr_ccdp.py` script uses the [`ccdproc.cosmicray_lacosmic`](https://ccdproc.readthedocs.io/en/latest/api/ccdproc.cosmicray_lacosmic.html#ccdproc.cosmicray_lacosmic) module to remove cosmic rays. The `ccdproc.cosmicray_lacosmic` module is based on the L.A.Cosmic algorithm for Laplacian Cosmic Ray Identification by Pieter G. van Dokkum (Yale) from http://www.astro.yale.edu/dokkum/lacosmic/. L.A.Cosmic detects cosmic rays of arbitrary shapes and sizes, and distinguishes between undersampled point sources and cosmic rays. If you use this program please refer to P. G. van Dokkum, 2001, PASP, 113, 1420.  

Please see `COPYRIGHT` file and `pyfosc$doc/LICENSES` directory for detailed copyright information.  

## 4. How to cite

Yuming Fu. (2024, April 13). PyFOSC: a pipeline toolbox for BFOSC/YFOSC long-slit spectroscopy data reduction (Version v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.10967240

Cite the current version (v1.1.0) with bibtex:  

```
@software{yuming_fu_2024_10967240,
  author       = {Yuming Fu},
  title        = {{PyFOSC: a pipeline toolbox for BFOSC/YFOSC long-slit spectroscopy data reduction}},
  month        = apr,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.1.0},
  doi          = {10.5281/zenodo.10967240},
  url          = {https://doi.org/10.5281/zenodo.10967240}
}
```
