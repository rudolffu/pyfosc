#!/usr/bin/env python
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack,Table
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
import glob
import os
from pathlib import Path
import shutil
import multiprocess as mp
import json
import datetime
from astropy.time import Time
import warnings
from pyraf import iraf

def headertable(filename, telescope=None):
    basename = os.path.basename(filename)
    hdr = fits.getheader(filename)
    arr = np.asarray(hdr.cards).T
    tab = Table(data=arr[1], names=arr[0])
    if telescope=='XLT':
        tab.add_column(basename, name='FILENAME')
        fields = ['OBJECT','RA','DEC','IMAGETYP',
                  'DATE-OBS','EXPTIME','OBSTYPE',
                  'FILTER','FILENAME','TELESCOP',
                  'INSTRUME','JD','AIRMASS','OBSERVAT']
    elif telescope=='LJT':
        tab.add_column(basename, name='FILENAME')
        fields = ['OBJECT','CAT-RA','CAT-DEC',
                  'FILTER1','FILTER3',
                  'DATE-OBS','EXPTIME','OBSTYPE',
                  'FILTER','FILENAME','TELESCOP',
                  'INSTRUME','MJD-OBS','AIRMASS',
                  'SITE','GAIN','RDNOISE']
    elif telescope=='HCT':
        fields = ['OBJECT','RA','DEC','IMAGETYP',
                  'DATE-OBS','EXPTIME','GRISM',
                  'FILTER','FILENAME','TELESCOP',
                  'INSTRUME','OBSERVAT']
    return tab[fields]

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    linelists = 'linelists$fear.dat'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    linelists = 'linelists$henear.dat'
elif teles == "HCT":
    print("Settings for HCT will be used.")
    linelists = 'linelists$fear.dat'
else:
    print("Error detected.")

flist = glob.glob('af*fits')
tablist = []
for item in flist:
    tablist.append(headertable(item, teles))
tb = vstack(tablist)
tb.sort(['DATE-OBS'])
df = tb.to_pandas()
if teles == "LJT":
    idx_fear = df.FILTER.str.contains('lamp_fe_argon')
    idx_hene = df.FILTER.str.contains('lamp_neon_helium')
    print("Fe-Ar lamp spectra:")
    print(df.loc[idx_fear, 'FILENAME'])
    iraf.onedspec()
    iraf.onedspec.identify.unlearn()
    iraf.onedspec.identify.fwidth = 10
    iraf.onedspec.identify.coordli = 'linelists$fear.dat'
    file_str = df.loc[idx_fear, 'FILENAME'].str.cat(sep=",")
    # for item in df.loc[idx_fear, 'FILENAME']:
    iraf.onedspec.identify(images=file_str)
    print("He-Ne lamp spectra:")
    print(df.loc[idx_hene, 'FILENAME'])
    iraf.onedspec()
    iraf.onedspec.identify.unlearn()
    iraf.onedspec.identify.fwidth = 10
    iraf.onedspec.identify.coordli = 'linelists$henear.dat'
    file_str = df.loc[idx_hene, 'FILENAME'].str.cat(sep=",")
    # for item in df.loc[idx_hene, 'FILENAME']:
    iraf.onedspec.identify(images=file_str)
else:
    iraf.onedspec()
    iraf.onedspec.identify.unlearn()
    iraf.onedspec.identify.fwidth = 10
    iraf.onedspec.identify.coordli = linelists
    iraf.onedspec.identify(images='af*fits')
