#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import os
from pathlib import Path
from shutil import copy2
import platform
from packaging import version
from pyfosc.paths import find_data_root
if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

CWD = Path.cwd()


DATA_ROOT = find_data_root(module_path=Path(__file__), cwd=CWD)
dbpath = DATA_ROOT / 'database'

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    lamplist = sorted(p.name for p in dbpath.glob('lamp_XL*fits'))
    linelists = 'linelists$fear.dat'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    lamplist = sorted(p.name for p in dbpath.glob('lamp_LJ*fits'))
    linelists = 'linelists$henear.dat'
elif teles == "HCT":
    print("Settings for HCT will be used.")
    lamplist = sorted(p.name for p in dbpath.glob('lamp_HCT*fits'))
    linelists = 'linelists$fear.dat'
else:
    print("Error detected.")

print('Possible lamp spectrum(a) for references:\n' +
      ", ".join(p for p in lamplist))
refspec = str(
    raw_input("Enter filename of the lamp spectrum you want to use: "))
refspec1 = Path(refspec).with_suffix('.fits').name
refspec_path = dbpath / refspec1
if not refspec_path.exists():
    raise FileNotFoundError(f"{refspec1} not found in {dbpath}")
copy2(refspec_path, CWD)
target_db = CWD / 'database'
target_db.mkdir(exist_ok=True)
idfile = dbpath / f"id{Path(refspec).stem}"
copy2(idfile, target_db)
iraf.onedspec()
iraf.onedspec.reidentify.unlearn()
# iraf.onedspec.reidentify.fwidth = 10
iraf.onedspec.reidentify.coordli = linelists
iraf.onedspec.reidentify(reference=refspec1,
                         images='af*.fit*')

iraf.onedspec.identify.unlearn()
iraf.onedspec.identify.fwidth = 10
iraf.onedspec.identify.coordli = linelists
iraf.onedspec.identify(images='af*.fit*')
print('---DONE---')
