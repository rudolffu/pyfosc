#!/usr/bin/env python
import glob
import os
import shutil
# import pyds9
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

# for file in flatlist:
#     os.remove(file)
#
# print('openning a ds9 window if not already openned...')
# pyds9.DS9()
#
# # check output flat image
# print('Check output file:')
# iraf.imstatistics('Flat')
# print(' Running "imexamine" task..')
# iraf.imexamine('Flat', 1)

print('--- DONE ---')
