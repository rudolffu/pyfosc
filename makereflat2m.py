#!/usr/bin/env python
import glob
# import pyds9
from pyraf import iraf


print('Loading IRAF packages ...')
iraf.imred()
iraf.ccdred()
iraf.specred()
iraf.twodspec()
iraf.twodspec.longslit()

print('unlearning previous settings...')
iraf.ccdred.unlearn()
iraf.ccdred.ccdproc.unlearn()
iraf.ccdred.combine.unlearn()
iraf.ccdred.zerocombine.unlearn()
iraf.ccdred.flatcombine.unlearn()
iraf.twodspec.longslit.response.unlearn()

# # regular expression of files (e.g bias_00*.fits, flat-2000jan01_?.*)
# theflat = str(raw_input('Enter flat image: '))
#
# iraf.specred.extinction = ''
# iraf.specred.caldir = ''
# iraf.specred.observatory = 'bao'

print('Create a response for flat...')
iraf.twodspec.longslit.response.interactive = True
iraf.twodspec.longslit.response.high_reject = 3.0
iraf.twodspec.longslit.response.low_reject = 3.0
iraf.twodspec.longslit.response(calibration='Flat', normalization='Flat',
                                response='reFlat')

# check output flat image
# print('openning a ds9 window if not already openned...')
# pyds9.DS9()
#
# print('Check output file:')
# iraf.imstatistics('nFlat')
# print('Running "imexamine" task..')
# iraf.imexamine('nFlat', 1)
print('Illumination normalization for flat...')
iraf.twodspec.longslit.illumination.unlearn()
iraf.twodspec.longslit.illumination.interac = 'yes'
iraf.twodspec.longslit.illumination(images='reFlat', illumina='ilFlat', low_rej=2,
                                    high_rej=2)

print('Make perfect Flat ^_^')
iraf.imarith.unlearn()
iraf.imarith(operand1='reFlat', op='/', operand2='ilFlat', result='perFlat')

print('--- DONE ---')
