#!/usr/bin/env python
# Load Python standard modules
import glob
import os
import shutil
from pyraf import iraf
import json

with open('myfosc.json') as file:
    settings = json.loads(file.read())

teles = settings['mysettings']['telescope']
Grism = settings['mysettings']['Grism']


print('Loading IRAF packages ...')
iraf.imred()
iraf.ccdred()

print('unlearning previous settings...')
iraf.ccdred.unlearn()
iraf.ccdred.ccdproc.unlearn()

print('Applying Bias, Overscan and pre Trimming...')
iraf.ccdred.ccdproc.ccdtype = ''
iraf.ccdred.ccdproc.noproc = False
iraf.ccdred.ccdproc.fixpix = False
iraf.ccdred.ccdproc.darkcor = False
iraf.ccdred.ccdproc.illumcor = False
iraf.ccdred.ccdproc.fringecor = False
iraf.ccdred.ccdproc.readcor = False
iraf.ccdred.ccdproc.scancor = False
iraf.ccdred.ccdproc.readaxis = 'line'
iraf.ccdred.ccdproc.flatcor = False

if teles == "XLT":
    print("Settings for XLT will be used.")
    iraf.ccdred.ccdproc.trimsec = '[51:1750,881:1350]'
    iraf.ccdred.ccdproc(images='@flatnall.list', overscan='no',
                        trim='yes', zerocor='yes', zero='Zero')
elif teles == "LJT":
    print("Settings for LJT will be used.")
    if Grism == "G3":
        iraf.ccdred.ccdproc.biassec = '[10:40,2301:4130]'
        iraf.ccdred.ccdproc.trimsec = '[751:1450,2301:4130]'
        # If you want more trimming on the blue side, use:
        # iraf.ccdred.ccdproc.biassec = '[10:40,2491:4130]'
        # iraf.ccdred.ccdproc.trimsec = '[751:1450,2491:4130]'
    elif Grism == "G10":
        iraf.ccdred.ccdproc.biassec = '[10:40,2591:3340]'
        iraf.ccdred.ccdproc.trimsec = '[751:1450,2591:3340]'
    elif Grism == "G14":
        iraf.ccdred.ccdproc.biassec = '[10:40,1901:4200]'
        iraf.ccdred.ccdproc.trimsec = '[751:1450,1901:4200]'
    iraf.ccdred.ccdproc(images='@flatnall.list', overscan='yes',
                        trim='yes', zerocor='yes', zero='Zero')
else:
    print("Error detected.")

print('--- DONE ---')
