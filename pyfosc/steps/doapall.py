#!/usr/bin/env python
import glob
from pyraf import iraf
import json
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
    disp_axis = 1
elif teles == "LJT":
    print("Settings for LJT will be used.")
    disp_axis = 2
elif teles == "HCT":
    print("Settings for HCT will be used.")
    disp_axis = 2
else:
    disp_axis = str(raw_input("Dispersion axis: 1 for line; 2 for column."))

print('Loading IRAF packages ...')
iraf.imred()
iraf.ccdred()
iraf.twodspec()
iraf.apextract()

print('unlearning previous settings...')
iraf.ccdred.unlearn()
iraf.ccdred.ccdproc.unlearn()
iraf.ccdred.combine.unlearn()
iraf.apextract.apall.unlearn()
iraf.apextract.dispaxis = disp_axis
iraf.apextract.verbose = 'no'

print('Extracting object aperture spectrum...')
iraf.apextract.apall.unlearn()
iraf.apextract.apall.readnoise = 'rdnoise'
iraf.apextract.apall.gain = 'gain'
iraf.apextract.apall.format = 'multispec'
iraf.apextract.apall.interac = True
iraf.apextract.apall.find = True
iraf.apextract.apall.recente = True
iraf.apextract.apall.resize = True
iraf.apextract.apall.edit = True
iraf.apextract.apall.trace = True
iraf.apextract.apall.fittrac = True
iraf.apextract.apall.extract = True
iraf.apextract.apall.extras = True
iraf.apextract.apall.review = True
iraf.apextract.apall.background = 'fit'
iraf.apextract.apall.pfit = 'fit2d'
iraf.apextract.apall.weights = 'variance'
iraf.apextract.apall(input='crf//@objall.list', output='acrf//@objall.list')

print('Extracting lamp spectrum...')
iraf.apextract.apall.unlearn()
iraf.apextract.apall.readnoise = 'rdnoise'
iraf.apextract.apall.gain = 'gain'
iraf.apextract.apall.format = 'onedspec'
iraf.apextract.apall.reference = 'crf//@objall.list'
iraf.apextract.apall.interac = False
iraf.apextract.apall.find = False
iraf.apextract.apall.recente = False
iraf.apextract.apall.resize = False
iraf.apextract.apall.edit = False
iraf.apextract.apall.trace = False
iraf.apextract.apall.fittrac = False
iraf.apextract.apall.extract = True
iraf.apextract.apall.extras = False
iraf.apextract.apall.review = True
iraf.apextract.apall.background = 'none'
iraf.apextract.apall.pfit = 'fit1d'
iraf.apextract.apall.weights = 'none'
iraf.apextract.apall(input='f//@lampall.list', output='af//@lampall.list')

print('--- DONE ---')

