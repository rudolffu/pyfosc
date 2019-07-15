# PyFOSC
PyFOSC is a pipeline toolbox for spectroscopic data reduction written in Python. It can be used for FOSC data from Xinglong/Lijiang 2-meter telescopes.

## Installation
This package depends on [IRAF](http://iraf.noao.edu/) and [PyRAF](http://www.stsci.edu/institute/software_hardware/pyraf). Although these two softwares will no longer be supported by their developers, you can still download and install them from [AstroConda](https://astroconda.readthedocs.io/en/latest/) using its 'Legacy Software Stack'.
### 1. Install IRAF and PyRAF from [astroconda](https://astroconda.readthedocs.io/en/latest/)

You can refer to the documentation of AstroConda for instructions on how to install [Legacy Software Stack with IRAF](https://astroconda.readthedocs.io/en/latest/installation.html#legacy-software-stack-with-iraf). The minimal installation can be done with:

```
conda config --add channels http://ssb.stsci.edu/astroconda
conda create -n iraf27 python=2.7 iraf-all pyraf-all
```
Note that every time you want to use IRAF/PyRAF and PyFOSC, you need to activate the env 'iraf27' that you have created:

```
source activate iraf27
```
And to deactivate, run:
```
source deactivate
```

### 2. Download PyFOSC and set environment variable for it.

You can use `git clone` to download this package.
```
git clone https://github.com/rudolffu/pyfosc.git
```
In order to run PyFOSC commands in the terminal, you need to add the path which PyFOSC is located to $PATH, by editting `~/.bashrc` (Linux, e.g., Ubuntu) or `~/.bash_profile`
