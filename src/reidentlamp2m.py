#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import os
import sys
from shutil import copy2
import platform
from packaging import version
if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

CWD = os.getcwd()
dbpath = os.path.join(os.path.dirname(sys.argv[0]), '../database')

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    lamplist = [os.path.basename(x) for x in glob.glob(
        dbpath + '/lamp_XL*fits')]
    linelists = 'linelists$fear.dat'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    lamplist = [os.path.basename(x) for x in glob.glob(
        dbpath + '/lamp_LJ*fits')]
    linelists = 'linelists$henear.dat'
else:
    print("Error detected.")

print('Possible lamp spectrum(a) for references:\n' +
      ", ".join(p for p in lamplist))
refspec = str(
    raw_input("Enter filename of the lamp spectrum you want to use: "))
refspec1 = refspec.strip('.fits') + '.fits'
copy2(dbpath + '/' + refspec1, CWD)
copy2(dbpath + '/id' + refspec.strip('.fits'), CWD + '/database')
iraf.onedspec()
iraf.onedspec.reidentify.unlearn()
# iraf.onedspec.reidentify.fwidth = 10
iraf.onedspec.reidentify.coordli = linelists
iraf.onedspec.reidentify(reference=refspec1,
                         images='af*fits')

iraf.onedspec.identify.unlearn()
iraf.onedspec.identify.fwidth = 10
iraf.onedspec.identify.coordli = linelists
iraf.onedspec.identify(images='af*fits')
print('---DONE---')
