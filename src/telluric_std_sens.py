#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import sys
import os
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
elif teles == "HCT":
    print("Settings for HCT will be used.")
    midname = 'HCT'
else:
    print("Error detected.")

iraf.images()
iraf.images.imutil()
print_fits_keyword_table(glob.glob("wacrf*fits"), keyword="OBJECT")

list_con = glob.glob('*con*fits')
list_tel = glob.glob('*tel*fits')
list_con.extend(list_tel)
if len(list_con) > 0:
    for f in list_con:
        os.remove(f)

stdspec = str(raw_input("Enter filename of the standard star spectrum: "))
starname = str(raw_input("Enter name of the standard star: "))
stdspec1 = stdspec.strip('.fits') + '.fits'
conname = "con"+starname+".fits"
settings['stdstarname'] = starname
with open('myfosc.json', 'w') as f:
    json.dump(settings,f)

iraf.onedspec()
iraf.onedspec.continuum.unlearn()
iraf.onedspec.continuum(input=stdspec1, output=conname)

hdu = fits.open(conname)
CRVAL1 = float(hdu[0].header['CRVAL1'])
CD1_1 = float(hdu[0].header['CD1_1'])
CRPIX1 = float(hdu[0].header['CRPIX1'])

if teles == "LJT":
    pix1=int((6750-CRVAL1)/CD1_1+CRPIX1-1)
    pix2=int((10000-CRVAL1)/CD1_1+CRPIX1-1)
else:
    pix1=int((6750-CRVAL1)/CD1_1+CRPIX1-1)
    pix2=int((8000-CRVAL1)/CD1_1+CRPIX1-1)
if pix2>max(hdu[0].data.shape):
    pix2 = max(hdu[0].data.shape)
cpsec = "["+str(pix1)+":"+str(pix2)+"]"

iraf.images.imutil.imcopy.unlearn()
iraf.images.imutil.imcopy(input=conname+cpsec, output="cal"+conname)

output_name = stdspec.strip('.fits') + '_tellcorr.fits'
iraf.onedspec.telluric.unlearn()
iraf.onedspec.telluric(input=stdspec,
                       output=output_name,
                       cal="cal"+conname,
                       answer="YES")
