#!/usr/bin/env python
import os
from astropy.nddata import CCDData
import ccdproc as ccdp
from astropy import units as u
from packaging import version
import platform
import sys
import json
import shutil
import glob
import re
from astropy.io import fits
import numpy as np

if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    gain = 2.2
    readnoise = 7.8
elif teles == "LJT":
    print("Settings for LJT will be used.")
    gain = 1.1
    readnoise = 6.6
elif teles == "HCT":
    print("Settings for HCT will be used.")
    gain = 1.22
    readnoise = 4.8
else:
    print("Error detected.")

objall = []
with open('objall.list') as file:
    for line in file:
        line = 'f'+line.strip('\n')
        objall.append(line)

# Image file collection
keys = ['file','imagetyp','obstype','filter','exptime','object']
files = ccdp.ImageFileCollection('./', keywords=keys,filenames=objall)
co_targets = files.filter()

# Process standard star and targets separately
print('Spectra after flat field correction:\n' + ", ".join(p for p in objall))
co_targets.summary
stdspec = str(raw_input("Enter filename of the standard star spectra \
(separated by ',' for more than 1): "))
stdlist = re.split(',|,\s+', stdspec)
stdexptime = []
for stdspec in stdlist:
    stdspec = stdspec.strip('.fits') + '.fits'
    hdr = fits.getheader(stdspec)
    exptime = float(hdr['EXPTIME'])
    stdexptime.append(exptime)
stdexptime = np.array(stdexptime)
maxexptime = np.max(stdexptime)
idx_tar = (co_targets.summary['exptime'] > maxexptime)
tbex = co_targets.summary[co_targets.summary['exptime'] <= maxexptime]

# For standard stars, done
for i, pathname in enumerate(tbex['file']):
    head, tail = os.path.split(pathname)
    newtail = 'cr'+tail
    newpname = os.path.join(head, newtail)
    print('Object: {}'.format(tbex['object'][i]))
    print('Copying {} to {}'.format(pathname, newpname))
    shutil.copy(pathname, newpname)

# CCD mask
flat1 = CCDData.read('./perFlat.fits', unit='adu')
fmask = ccdp.ccdmask(flat1)

# For targets
tb_tar = co_targets.summary[idx_tar]
try:
    hdr = fits.getheader(tb_tar['file'].tolist()[0])
    gain = hdr['GAIN']
    readnoise = hdr['RDNOISE']
except:
    pass
for file_path in tb_tar['file'].tolist():
    file_name = os.path.basename(file_path)
    ccd = CCDData.read(file_path, unit='adu')
    ccd = ccdp.gain_correct(ccd, gain * u.electron / u.adu)
    ccd.mask = fmask
    new_ccd = ccdp.cosmicray_lacosmic(ccd, readnoise=readnoise, sigclip=7, verbose=True)
    new_ccd.mask = None
    new_ccd.write('cr'+file_name, overwrite=True)
