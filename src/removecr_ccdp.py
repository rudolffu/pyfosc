#!/usr/bin/env python
import os
from astropy.nddata import fits_ccddata_reader,CCDData
import ccdproc as ccdp
from astropy import units as u
from packaging import version
import platform
import sys
import json
import shutil

if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

basepath = os.path.dirname(sys.argv[0])
lacos_im = os.path.join(basepath, '../iraf_tasks/lacos_im.cl')
with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    gain = 1.41
    readnoise = 4.64
elif teles == "LJT":
    print("Settings for LJT will be used.")
    gain = 3.20
    readnoise = 13.5
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
idx1 = (co_targets.summary['exptime'] > 500)
tbex = co_targets.summary[co_targets.summary['exptime'] < 900]

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
tb1 = co_targets.summary[idx1]
for file_path in tb1['file'].tolist():
    file_name = os.path.basename(file_path)
    ccd = CCDData.read(file_path, unit='adu')
    ccd = ccdp.gain_correct(ccd, gain * u.electron / u.adu)
    ccd.mask = fmask
    new_ccd = ccdp.cosmicray_lacosmic(ccd, readnoise=readnoise, sigclip=7, verbose=True)
    cr_mask = new_ccd.mask
    cr_mask[ccd.mask] = False
    new_ccd.mask = None
    new_ccd.write('cr'+file_name, overwrite=True)
