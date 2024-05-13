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

basepath = os.path.dirname(sys.argv[0])
myonedstds = os.path.join(basepath, '../iraf_data/onedstds')
myonedstds = os.path.abspath(myonedstds)
with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    extinct = os.path.join(basepath, '../extinction/baoextinct.dat')
    midname = 'XL'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    extinct = os.path.join(basepath, '../extinction/LJextinct.dat')
    midname = 'LJ'
elif teles == "HCT":
    print("Settings for HCT will be used.")
    extinct = os.path.join(basepath, '../extinction/LJextinct.dat')
    midname = 'HCT'
else:
    print("Error detected.")

# check if the following files exist, if true, delete them
list_wa = glob.glob('wa*fits')
list_cwa = glob.glob('cwa*fits')
list_tel = glob.glob('*tel*fits')
list_to_del = list_wa + list_cwa + list_tel
list_to_del.append('std')
list_to_del.append('sens.fits')
list_to_del.append('con*.fits')
list_to_del.append('calcon*fits')
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

inputlist = glob.glob('a*fits')
iraf.onedspec()
print('Add refspectra to header...')
iraf.onedspec.refspectra.unlearn()
for img in inputlist:
    iraf.onedspec.refspectra(input=img, referen='af*fits',
                             sort='', group='', time='no')

print('Dispersion correction...')
iraf.onedspec.dispcor.unlearn()
for img in inputlist:
    iraf.onedspec.dispcor(input=img, output='w' + img)

olist1 = glob.glob('w*fits')
print('Spectra after wavelength calibration:\n' + ", ".join(p for p in olist1))
iraf.images()
iraf.images.imutil()
iraf.images.imutil.imheader(images="wacrf*fits")
stdspec = str(raw_input("Enter filename of the standard star spectrum: "))
starname = str(raw_input("Enter name of the standard star: "))
starname = starname.lower()
stdspec1 = stdspec.strip('.fits') + '.fits'

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
olist2 = glob.glob('wac*fits')
olist2.remove(stdspec1)
iraf.twodspec.longslit.calibrate.unlearn()
iraf.twodspec.longslit.calibrate.extinction = extinct
for obj in olist2:
    iraf.twodspec.longslit.calibrate(input=obj,
                                     output='c' + obj,
                                     sensitivity='sens')

print('---DONE---')
