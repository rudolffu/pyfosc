#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits
import json
import sys
import os

basepath = os.path.dirname(sys.argv[0])
myonedstds = os.path.join(basepath, '../iraf_data/onedstds')
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
else:
    print("Error detected.")

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
stdspec = str(raw_input("Enter filename of the standard star spectrum: "))
starname = str(raw_input("Enter name of the standard star: "))
stdspec1 = stdspec.strip('.fits') + '.fits'

stdlist = glob.glob(myonedstds + "/**/" + starname + ".dat")
dirnames = [os.path.basename(os.path.dirname(x)) for x in stdlist]
print("The \"onedstds\" dirs containing your standard star:\n" +
      ", ".join(p for p in dirnames))
stddir = str(raw_input("Enter name of the onedstds directory: "))
mycaldir = "onedstds$" + stddir + '/'


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

print('Copy to onedspec...')
iraf.twodspec.longslit.scopy.unlearn()
olist3 = glob.glob('cwac*fits')
for obj in olist3:
    objind = str(olist3.index(obj) + 1)
    objname = fits.getheader(obj)['OBJECT'] + \
        '_' + midname + objind + '_' + starname
    iraf.twodspec.longslit.scopy(input=obj, output=objname, bands=1,
                                 format='onedspec')
print('---DONE---')
