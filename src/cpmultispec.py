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


starname = str(raw_input("Enter name of the standard star: "))

iraf.twodspec()
iraf.twodspec.longslit()
iraf.twodspec.longslit.scopy.unlearn()
olist3 = glob.glob('cwac*fits')
for obj in olist3:
    objind = str(olist3.index(obj) + 1)
    objname = fits.getheader(obj)['OBJECT'] + \
        '_' + midname + objind + '_' + starname + "_4d"
    iraf.twodspec.longslit.scopy(input=obj, output=objname, bands="",
                                 format='multispec')
print('---DONE---')