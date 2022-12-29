#!/usr/bin/env python
from sys import argv
import glob
from pyraf import iraf
import multiprocessing as mp
import os
import shutil

imglist = glob.glob('YF*fits')
raw_dir = './raw'
if not os.path.exists(raw_dir):
    os.makedirs(raw_dir)


def cpyfext(img):
    output = img[6:15]
    iraf.imcopy(input=img + '[1]', output=output, mode='hl', verbose='no')
    shutil.move(img, raw_dir + '/' + img)


pool = mp.Pool(mp.cpu_count())
pool.map(cpyfext, imglist)

print('Copying raw files finished.')
