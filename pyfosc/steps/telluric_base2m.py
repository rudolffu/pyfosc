#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import sys
import os
import platform
from packaging import version
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
elif teles == "P200":
    print("Settings for P200 will be used.")
    side = settings['mysettings']['side']
    midname = 'P200'+side
else:
    print("Error detected.")

iraf.images()
iraf.images.imutil()
iraf.images.imutil.imheader(images="wacrf*fits")

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

if teles == "LJT" or teles == "P200":
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

olist3 = glob.glob('cwac*fits')
for obj in olist3:
    objind = str(olist3.index(obj) + 1)
    objname = fits.getheader(obj)['OBJECT'] + \
        '_' + midname + objind + '_' + starname + "_tel"
    iraf.onedspec.telluric.unlearn()
    iraf.onedspec.telluric(input=obj,
                           output=objname,
                           cal="cal"+conname,
                           answer="YES")

print('---DONE---')

