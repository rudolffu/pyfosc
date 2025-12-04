#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import sys
import os
import platform
from pathlib import Path
from packaging import version
from pyfosc.paths import find_data_root
if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

DATA_ROOT = find_data_root(module_path=Path(__file__), cwd=Path.cwd())
myonedstds = os.path.abspath(os.path.join(DATA_ROOT, 'iraf_data/onedstds'))
with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    extinct = os.path.join(DATA_ROOT, 'extinction/baoextinct.dat')
    midname = 'XL'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    extinct = os.path.join(DATA_ROOT, 'extinction/LJextinct.dat')
    midname = 'LJ'
elif teles == "HCT":
    print("Settings for HCT will be used.")
    extinct = os.path.join(DATA_ROOT, 'extinction/LJextinct.dat')
    midname = 'HCT'
elif teles == "P200":
    print("Settings for P200 will be used.")
    side = settings['mysettings']['side']
    extinct = 'onedstds$kpnoextinct.dat'
    midname = 'P200' + side
else:
    print("Error detected.")

# check if the following files exist, if true, delete them
list_wa = glob.glob('wa*.fit*')
list_cwa = glob.glob('cwa*.fit*')
list_tel = glob.glob('*tel*.fit*')
list_to_del = list_wa + list_cwa + list_tel
list_to_del.append('std')
list_to_del.extend(glob.glob('sens.fit*'))
list_to_del.extend(glob.glob('con*.fit*'))
list_to_del.extend(glob.glob('calcon*.fit*'))
if len(list_wa) > 0 or len(list_cwa) > 0 or len(list_tel) > 0:
    delete = input("Files that are wave/flux calibrated exist. Delete them? (y/n): ")
    if delete == 'y':
        for f in list_to_del:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
    else:
        print("Please rename or delete the files and run the script again.")
        sys.exit(1)

inputlist = glob.glob('a*.fit*')
iraf.onedspec()
print('Add refspectra to header...')
iraf.onedspec.refspectra.unlearn()
for img in inputlist:
    iraf.onedspec.refspectra(input=img, referen='af*.fit*',
                             sort='', group='', time='no')

print('Dispersion correction...')
iraf.onedspec.dispcor.unlearn()
for img in inputlist:
    iraf.onedspec.dispcor(input=img, output='w' + img)

olist1 = glob.glob('w*.fit*')
print('Spectra after wavelength calibration:\n' + ", ".join(p for p in olist1))
iraf.images()
iraf.images.imutil()
iraf.images.imutil.imheader(images="wacrf*.fit*")
stdspec = str(raw_input("Enter filename of the standard star spectrum: "))
starname = str(raw_input("Enter name of the standard star: "))
starname = starname.lower()
stdspec_name = Path(stdspec).name

stdlist = glob.glob(myonedstds + "/**/" + starname + ".dat")
dirnames = [os.path.basename(os.path.dirname(x)) for x in stdlist]
print("The \"onedstds\" dirs containing your standard star:\n" +
      ", ".join(p for p in dirnames))
stddir = str(raw_input("Enter name of the onedstds directory: "))
# mycaldir = "onedstds$" + stddir + '/'
mycaldir = f"{myonedstds}/{stddir}/"

iraf.twodspec()
iraf.twodspec.longslit()
print('Add std info to standard star spectrum...')
iraf.twodspec.longslit.standard.unlearn()
iraf.twodspec.longslit.standard(input=stdspec, output='std',
                                extinct=extinct,
                                caldir=mycaldir,
                                star_name=starname)
print('Make sensfunc...')
iraf.twodspec.longslit.sensfunc.unlearn()
iraf.twodspec.longslit.sensfunc(
    standard='std', sensitiv='sens', extinct=extinct)

print('Flux calibration...')
olist2 = [Path(p).name for p in glob.glob('wac*.fit*')]
try:
    olist2.remove(stdspec_name)
except ValueError:
    pass
iraf.twodspec.longslit.calibrate.unlearn()
iraf.twodspec.longslit.calibrate.extinction = extinct
for obj in olist2:
    iraf.twodspec.longslit.calibrate(input=obj,
                                     output='c' + obj,
                                     sensitivity='sens')

print('---DONE---')
