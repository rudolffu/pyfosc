#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import os
import sys
from shutil import copy2

CWD = os.getcwd()
basepath = os.path.dirname(sys.argv[0])
with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    lamplist = [os.path.basename(x) for x in glob.glob(
        basepath + '/database/lamp_XL*fits')]
elif teles == "LJT":
    print("Settings for LJT will be used.")
    lamplist = [os.path.basename(x) for x in glob.glob(
        basepath + '/database/lamp_LJ*fits')]
else:
    print("Error detected.")

print('Possible lamp spectrum(a) for references:\n' +
      ", ".join(p for p in lamplist))
refspec = str(
    raw_input("Enter filename of the lamp spectrum you want to use: "))
refspec1 = refspec.strip('.fits') + '.fits'
copy2(basepath + '/database/' + refspec1, CWD)
copy2(basepath + '/database/id' + refspec.strip('.fits'), CWD + '/database')
iraf.onedspec()
iraf.onedspec.reidentify.unlearn()
# iraf.onedspec.reidentify.fwidth = 10
iraf.onedspec.reidentify.coordli = 'linelists$henear.dat'
iraf.onedspec.reidentify(reference=refspec1,
                         images='af*fits')
print('---DONE---')
