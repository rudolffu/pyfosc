#!/usr/bin/env python
import matplotlib as mpl
import glob
import matplotlib.pyplot as plt
import multiprocessing as mp
from astropy.io import fits
import numpy as np

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


class Spiraf():
    def __init__(self, fname):
        hdu = fits.open(fname)
        self.name = fits.getheader(fname)['object']
        CRVAL1 = hdu[0].header['CRVAL1']
        CD1_1 = hdu[0].header['CD1_1']
        CRPIX1 = hdu[0].header['CRPIX1']
        l = len(hdu[0].data)
        self.wave = np.linspace(CRVAL1, CRVAL1 + (l - CRPIX1) * CD1_1, l)
        self.flux = hdu[0].data

    def plot(self, axlim='auto'):
        plt.figure(figsize=(8, 6))
        plt.plot(self.wave, self.flux, lw=1)
        plt.xlim([3800, 9000])
        plt.ylim(bottom=0)
        plt.xlabel(r'$\mathrm{Wavelength(\AA})$')
        plt.ylabel(r'$\mathrm{Flux(erg/s/cm^{2}/\AA)}$')
        plt.title('Spectrum of ' + self.name)
        plt.savefig(self.name + '_spec.pdf', dpi=300)


sfiles = glob.glob('J*fits')
for sp in sfiles:
    spec = Spiraf(sp)
    spec.plot()
