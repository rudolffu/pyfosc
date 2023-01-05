#!/usr/bin/env python
import glob
from pyraf import iraf
import json
from astropy.io import fits

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

iraf.onedspec()
iraf.onedspec.identify.unlearn()
iraf.onedspec.identify.fwidth = 10
iraf.onedspec.identify.coordli = linelists
iraf.onedspec.identify(images='af*fits')
