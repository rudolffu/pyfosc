#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import sys
import os

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

iraf.images()
iraf.images.imutil()
iraf.images.imutil.imheader(images="wacrf*fits")

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
CRVAL1 = hdu[0].header['CRVAL1']
CD1_1 = hdu[0].header['CD1_1']
CRPIX1 = hdu[0].header['CRPIX1']

pix1=int((6750-CRVAL1)/CD1_1+CRPIX1-1)
pix2=int((8300-CRVAL1)/CD1_1+CRPIX1-1)
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

print('Copy to onedspec...')
olist4 = glob.glob('*_*tel*fits')
print(olist4)

iraf.twodspec()
iraf.twodspec.longslit()
iraf.twodspec.longslit.scopy.unlearn()
for obj in olist4:
    objind = str(olist4.index(obj) + 1)
    objname = fits.getheader(obj)['OBJECT'] + \
        '_' + midname + objind + '_' + starname
    iraf.twodspec.longslit.scopy(input=obj, output=objname, bands=1,
                                 format='onedspec')
print('---DONE---')
