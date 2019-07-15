#!/usr/bin/env python
import glob
# import pyds9
from pyraf import iraf


print('Loading IRAF packages ...')
iraf.imred()
iraf.ccdred()

print('unlearning previous settings...')
iraf.ccdred.unlearn()
iraf.ccdred.ccdproc.unlearn()

print('Dividing perFlat...')
iraf.ccdred.ccdproc.ccdtype = ''
iraf.ccdred.ccdproc.noproc = False
iraf.ccdred.ccdproc.fixpix = False
iraf.ccdred.ccdproc.overscan = False
iraf.ccdred.ccdproc.darkcor = False
iraf.ccdred.ccdproc.illumcor = False
iraf.ccdred.ccdproc.fringecor = False
iraf.ccdred.ccdproc.readcor = False
iraf.ccdred.ccdproc.scancor = False
iraf.ccdred.ccdproc.readaxis = 'line'
iraf.ccdred.ccdproc.zerocor = False
iraf.ccdred.ccdproc.flatcor = True
iraf.ccdred.ccdproc.flat = 'perFlat'
iraf.ccdred.ccdproc.trim = False
iraf.ccdred.ccdproc(images='@specall.list', output='f//@specall.list')

print('--- DONE ---')
