#!/usr/bin/env python
import glob
import os
import shutil
from pyraf import iraf


print('Loading IRAF packages ...')
iraf.imred()
iraf.ccdred()

print('unlearning previous settings...')
iraf.ccdred.unlearn()
iraf.ccdred.ccdproc.unlearn()
iraf.ccdred.combine.unlearn()
iraf.ccdred.flatcombine.unlearn()

# combine flat images
print('Combining flat images ...')
iraf.ccdred.flatcombine.ccdtype = ''
iraf.ccdred.flatcombine.process = 'no'
iraf.ccdred.flatcombine.reject = 'avsigclip'
iraf.ccdred.flatcombine.rdnoise = 'rdnoise'
iraf.ccdred.flatcombine.gain = 'gain'
iraf.ccdred.flatcombine.output = 'Flat'
iraf.ccdred.flatcombine(input='@flat.list')

print('--- DONE ---')
