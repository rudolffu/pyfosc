#!/usr/bin/env python
import glob
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

print('Create a response for flat...')
iraf.twodspec.longslit.response.interactive = True
iraf.twodspec.longslit.response.high_reject = 3.0
iraf.twodspec.longslit.response.low_reject = 3.0
iraf.twodspec.longslit.response(calibration='Flat', normalization='Flat',
                                response='reFlat')

print('Illumination normalization for flat...')
iraf.twodspec.longslit.illumination.unlearn()
iraf.twodspec.longslit.illumination.interac = 'yes'
iraf.twodspec.longslit.illumination(images='reFlat', illumina='ilFlat', low_rej=2,
                                    high_rej=2)

print('Make perfect Flat ^_^')
iraf.imarith.unlearn()
iraf.imarith(operand1='reFlat', op='/', operand2='ilFlat', result='perFlat')

print('--- DONE ---')

