# PyFOSC
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10967240.svg)](https://doi.org/10.5281/zenodo.10967240)  

__Please note: This program is under active development. The documentation may be inconsistent with the latest codes.__

__IRAF/PyRAF installation guide for Mac with Apple Silicon (M1/M2) can be found [here](https://yumingfu.space/tech/2024-iraf-mac-installation/).__

PyFOSC is a pipeline toolbox for long-slit spectroscopy data reduction written in Python. It can be used for FOSC (Faint Object Spectrograph and Camera) data from Xinglong/Lijiang 2-meter telescopes.  

BFOSC (Beijing-Faint Object Spectrograph and Camera) is an instrument of the 2.16-m Telescope in Xinglong Observatory, National Astronomical Observatories, Chinese Academy of Sciences (NAOC) (IAU code: 327, coordinates: 40°23′39″ N, 117°34′30″ E). For more information about BFOSC, please see http://www.xinglong-naoc.cn/html/en/gcyq/216/detail-18.html. The details of the Xinglong 2.16-m telescope and BFOSC are also reported in [Fan et al. 2016](https://ui.adsabs.harvard.edu/abs/2016PASP..128k5005F/abstract) and [Zhao et al. 2018](https://ui.adsabs.harvard.edu/abs/2018RAA....18..110Z/abstract).

YFOSC (Yunnan-Faint Object Spectrograph and Camera) is an instrument of the 2.4-m Telescope in Lijiang Observatory, Yunnan Observatories, Chinese Academy of Sciences (YNAO) (IAU code: O44, coordinates: 26°42′33.1″ N, 100°1′51.6″ E). For more information about YFOSC, please see http://wiki.gmg.org.cn/pages/viewpage.action?pageId=557106. The details of the Lijiang 2.4-m telescope and YFOSC are also reported in [Wang et al. 2019](https://ui.adsabs.harvard.edu/abs/2019RAA....19..149W/abstract).

## 1. Installation 

### 1.1. Dependencies
```
IRAF
PyRAF
pandas
ccdproc
specutils
specreduce
```

`pandas` is already included in the Anaconda distribution. To install `ccdproc` with `conda`, you can use:
```sh
conda install conda-forge::ccdproc
```

`specutils` is a package for representing and manipulating spectroscopic data in Python. You can install it with:
```sh
conda install conda-forge::specutils
```
or 
```sh
pip install specutils
```

`specreduce` is a package for tracing and extracting 1d spectra. You can install it with:
```sh
pip install git+https://github.com/astropy/specreduce.git
```

This package depends on [IRAF](http://iraf.noao.edu/) and [PyRAF](http://www.stsci.edu/institute/software_hardware/pyraf). You can download and install the [IRAF Community Distribution](https://iraf-community.github.io/).

### 1.2. Install the [IRAF/PyRAF Community Distribution](https://iraf-community.github.io/)

For detailed description on how to install IRAF and PyRAF, visit:
- https://iraf-community.github.io/install.html
- https://iraf-community.github.io/x11iraf.html
- https://iraf-community.github.io/pyraf.html

If you are using a Mac with Apple Silicon (M1/M2), you can follow the instructions in [this blog post](https://yumingfu.space/tech/2024-iraf-mac-installation/) to install IRAF/PyRAF.

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

# Run reductions: pre-wavecal (ccdproc/astropy) + wavecal/telluric (IRAF)
pyfosc run
```
Notes:
- Raw files are backed up under `raw/` and renamed copies placed in `data/`.
- Intermediate products (bias/trim `_bc`, `f*`, `crf*`) are written under `data/`.
- Master frames are saved under `MasterFrames/` at the repository root.

### 2.2. Running the pipeline (currently in a python + pyraf hybrid mode)

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

#### 2.2.2. Wave/flux/telluric (IRAF)

These steps are run by `pyfosc run`. To execute manually, run from the `data/` directory:

The remaining steps include:
- Wavelength calibration
- Flux calibration
- Telluric correction

These steps can be done with the following scripts:

```bash
cd data/ # Go to data directory
python -m pyfosc.steps.reidentlamp2m   # Reidentify lamp spectra with stored references
# or: python -m pyfosc.steps.identlamp2m  # Identify lamps manually
python -m pyfosc.steps.wavecal2m       # Wavelength + flux calibration
python -m pyfosc.steps.telluric_base2m # Telluric correction and 1D spectra
```

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
