#!/usr/bin/env python
# Load Python standard modules
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
iraf.ccdred.zerocombine.unlearn()

# combine bias images
print('Combining zero level images...')
iraf.ccdred.zerocombine.ccdtype = ''
iraf.ccdred.zerocombine.reject = 'ccdclip'
iraf.ccdred.zerocombine.rdnoise = 'rdnoise'
iraf.ccdred.zerocombine.gain = 'gain'
iraf.ccdred.zerocombine(input='@zero.list')

print('--- DONE ---')
