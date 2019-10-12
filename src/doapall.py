#!/usr/bin/env python
import glob
from pyraf import iraf


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
iraf.apextract.dispaxis = 1
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
