#!/usr/bin/env python
import glob
import ccdproc as ccdp
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData


# do flat correction on all the images in specall.list
# read in the list of images as a list of strings
with open('specall.list', 'r') as f:
    images = f.readlines()
images = [x.strip() for x in images]
specall = ImageFileCollection('./',filenames=images)
good_flat = CCDData.read('perFlat.fits', unit='adu')

for hdu, fname in specall.hdus(return_fname=True):
    # read in the image
    img = ccdp.CCDData.read(fname, unit='adu')
    img = ccdp.flat_correct(img, good_flat)
    # write the image to disk
    img.write('f'+fname, overwrite=True)
    print('Flat corrected', fname)
