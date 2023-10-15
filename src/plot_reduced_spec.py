#!/usr/bin/env python
import matplotlib as mpl
import glob
import matplotlib.pyplot as plt
import multiprocessing as mp
from astropy.io import fits
import numpy as np
import json
from specbox.basemodule import SpecIRAF
import sys
import pathlib

SMALL_SIZE = 14
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

# controls default text sizes
plt.rc('font', size=SMALL_SIZE, family='serif')
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
mpl.rcParams['mathtext.rm'] = 'serif'

# get the filename or expression of the reduced spectra from the user
input_arg = sys.argv[1]
inputlist = glob.glob(input_arg)
if len(inputlist) == 0:
    print("No files found.")
    sys.exit()
elif len(inputlist) == 1:
    print("Only one file found.")
    print("The file will be plotted.")
else:
    print("Multiple files found.")
    print("The files will be plotted.")
    print("The files are:\n" + ", ".join(p for p in inputlist))

pathlib.Path('plots').mkdir(parents=True, exist_ok=True)
for fname in inputlist:
    sp = SpecIRAF(fname)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax = sp.plot(ax=ax)
    ax.set_title(fname)
    plt.savefig('./plots/'+fname.strip('.fits') + '.png', dpi=300) 
    plt.show()   
