#!/usr/bin/env python
import glob
from pyraf import iraf
from astropy.io import fits

iraf.onedspec()
iraf.onedspec.identify.unlearn()
iraf.onedspec.identify.fwidth = 10
iraf.onedspec.identify.coordli = 'linelists$henear.dat'
iraf.onedspec.identify(images='af*fits')
