# PyFOSC
PyFOSC is a pipeline toolbox for spectroscopic data reduction written in Python. It can be used for FOSC data from Xinglong/Lijiang 2-meter telescopes.

## 1. Installation
### 1.1. Dependencies
```
IRAF
PyRAF
pandas
```
This package depends on [IRAF](http://iraf.noao.edu/) and [PyRAF](http://www.stsci.edu/institute/software_hardware/pyraf). Although these two softwares will no longer be supported by their developers, you can still download and install them from [AstroConda](https://astroconda.readthedocs.io/en/latest/) using its 'Legacy Software Stack'.
### 1.2. Install IRAF and PyRAF from [astroconda](https://astroconda.readthedocs.io/en/latest/)

For detailed description on how to install IRAF and PyRAF from astroconda, you can either go to the documentation of AstroConda's [Legacy Software Stack with IRAF](https://astroconda.readthedocs.io/en/latest/installation.html#legacy-software-stack-with-iraf), or my blog post [How to install and initiate IRAF in Ubuntu/Mac OS](https://rudolffu.github.io/tech/iraf-installation/). The minimal installation can be done with:

```Bash
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n iraf27 python=2.7 iraf-all pyraf-all stsci
conda install -n iraf27 pandas
```
Note that every time you want to use IRAF/PyRAF and PyFOSC, you need to activate the env 'iraf27' that you have created:

```Bash
source activate iraf27
```
And to deactivate, run:  
```Bash
source deactivate
```

### 1.3. Download PyFOSC and set environment variable for it.

You can use `git clone` to download this package.  
```bash
git clone https://github.com/rudolffu/pyfosc.git
```
In order to run PyFOSC commands in the terminal, you need to add the path of PyFOSC and its sub-directory `src` to $PATH, by editing `~/.bashrc` (Linux, e.g., Ubuntu) or `~/.bash_profile` (Mac OS). An example of this can be:
```Bash
export PATH=/Your/Path/to/pydbsp:$PATH
export PATH=/Your/Path/to/pydbsp/src:$PATH
```

## 2. Usage

### 2.1 Preparation

#### 2.1.1 Run `pyfosc_init` to begin

First go to the working directory which contains FOSC spectroscopic data, activate `iraf27` env and `mkiraf`.
```
source activate iraf27 (#or conda activate iraf27)
mkiraf
```  
Select `xgterm` when the terminal prompts. Then run `pyfosc_init` from `iraf27` env:
```
pyfosc_init
```

#### 2.1.2 Generate lists of files

As PyFOSC pipeline reduce data grouped by different types (bias, flat, object, etc), make sure you have lists of fits files as follows:
```
zero.list --------------- List of bias files (e.g. bias*.fits).
flat.list --------------- List of flat field files (e.g. flat*.fits).
objall.list ------------- List of 2d spectra images of all objects (science targets and standard stars).
lampall.list ------------ List of all 2d lamp spectra images.
flatnall.list ----------- List of 2d images that need zero correction (flat.list + lampall.list + objall.list).
specall.list ------------ List of all 2d spectra images (objall.list + lampall.list).
```

The `pyfosc$iraf_scripts` directory contains template `cl` scripts that can be executed in IRAF/PyRAF to generate these lists. Users can modify the scripts for specific cases.  

#### 2.1.3 Update headers for BFOSC data obtained with 2.16-m Telescope, Xinglong Station, NAOC (BAO).
Activate `iraf27` env and run `addheaderinfo.py`.
```bash
source activate iraf27
addheaderinfo.py
```

### 2.2 Running the pipeline
You can run `pyfosc_run.sh` from `iraf27` env:
```
pyfosc_run.sh
```
Alternatively, you can run the scripts step by step, following the order as:
```
makezero2m.py     # Combine zero(bias) frames.
ccdotz.py         # Do zero(bias), overscan correction and trimming.
makeflat2m.py     # Combine flat fields.
makereflat2m.py   # Do (illumination) normalization and get perfect flat.
divideflat2m.py   # Do flat correction.
removecr2m.py     # Remove cosmic rays in 2d images. Please go to 3. Credits for more information.
doapall.py        # Extract spectra.
reidentlamp2m.py  # Reidentify lamp spectra with previously stored ones.
#Can use identlamp2m.py instead to identify lamp oneself.
wavecal2m.py      # Do wavelength calibration, flux calibration and extract 1d spectra.
```

Plot all final 1d spectra with:
```
plotonedsp.py
```

Move some intermediate files into `INTMD` directory:
```
mvintmd.sh
```

## 3. Credits

This software contains sources from third-party softwares:  
The `pyfosc$iraf_data/onedstds` directory from `IRAF`, contains standard calibration data for extinction and sensitivity calibration.  
The `pyfosc$iraf_task/lacos_im.cl` script is the `IRAF` version for Laplacian Cosmic Ray Identification by Pieter G. van Dokkum (Yale) from http://www.astro.yale.edu/dokkum/lacosmic/. L.A.Cosmic is an algorithm for robust cosmic ray identification. It detects cosmic rays of arbitrary shapes and sizes, and distinguishes between undersampled point sources and cosmic rays. If you use this program please refer to P. G. van Dokkum, 2001, PASP, 113, 1420.  
Please see `COPYRIGHT` file and `pyfosc$doc/LICENSES` directory for detailed copyright information.  
