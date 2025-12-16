#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import sys
import os
from pathlib import Path 
import shutil
import re
import platform
from packaging import version
from pyfosc.fits_headers import print_fits_keyword_table
if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    midname = 'XL'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    midname = 'LJ'
else:
    print("Error detected.")


posans = {"", "y", "yes", "true"}
negans = {"n", "no", "false"}
answer = str(raw_input("Do you want to trim the spectrum(a) ([y]/n): "))
if answer.lower() in negans:
    print("Exit...")
    sys.exit(0)


try:
    starname = settings['stdstarname']
except:
    iraf.images()
    iraf.images.imutil()
    print_fits_keyword_table(glob.glob("wacrf*fits"), keyword="OBJECT")
    stdspec = str(raw_input("Enter filename of the standard star spectrum: "))
    starname = str(raw_input("Enter name of the standard star: "))
    settings['stdstarname'] = starname
    with open('myfosc.json', 'w') as f:
        json.dump(settings,f)


Path("./onedbak").mkdir(exist_ok=True)

try:
    onedlist1 = glob.glob('*_*0001.fits')
    for specfile in onedlist1:
        shutil.move(specfile, "./onedbak/")
    print("Old onedspec files backed up.")
except:
    pass


olist = glob.glob('*_*tel*fits')
lower = raw_input("Enter the new lower limit of wavelength in Angstrom (DEFAULT): ")
upper = raw_input("Enter the new upper limit of wavelength in Angstrom (DEFAULT): ")
hdu = fits.open(olist[0])
CRVAL1 = hdu[0].header['CRVAL1']
CD1_1 = hdu[0].header['CD1_1']
CRPIX1 = hdu[0].header['CRPIX1']

try:
    lower = float(lower)
    pix1=int((lower-CRVAL1)/CD1_1+CRPIX1-1)
except:
    pix1=1
    
try:
    upper = float(upper)
    pix2=int((upper-CRVAL1)/CD1_1+CRPIX1-1)
except:
    pix2=max(hdu[0].data.shape)


cpsec = "["+str(pix1)+":"+str(pix2)+"]"

iraf.images.imutil.imcopy.unlearn()
for telspec in olist:
    newname = re.sub('\.fits$', '', telspec) +"trim"
    iraf.images.imutil.imcopy(input=telspec+cpsec, 
                              output=newname)

olist2 = glob.glob('*_*teltrim.fits')

iraf.twodspec()
iraf.twodspec.longslit()
iraf.twodspec.longslit.scopy.unlearn()
for obj in olist2:
    objind = str(olist2.index(obj) + 1)
    objname = fits.getheader(obj)['OBJECT'] + \
        '_' + midname + objind + '_' + starname
    iraf.twodspec.longslit.scopy(input=obj, output=objname, bands=1,
                                 format='onedspec')
print('---DONE---')
