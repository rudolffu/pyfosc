#!/usr/bin/env python
# Load Python standard modules
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
iraf.ccdred.zerocombine.unlearn()

# combine bias images
print('Combining zero level images...')
iraf.ccdred.zerocombine.ccdtype = ''
iraf.ccdred.zerocombine.reject = 'ccdclip'
iraf.ccdred.zerocombine.rdnoise = 'rdnoise'
iraf.ccdred.zerocombine.gain = 'gain'
iraf.ccdred.zerocombine(input='@zero.list')

# for zero in zerolist:
#     os.remove(zero)
#
# print('openning a ds9 window if not already openned...')
# pyds9.DS9()
#
# # check output image
# print('Check output file:')
# iraf.imstatistics('Zero')
#
# print(' Running "imexamine" task..')
# iraf.imexamine('Zero', 1)

print('--- DONE ---')
