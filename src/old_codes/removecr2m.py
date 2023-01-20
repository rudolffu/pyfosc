#!/usr/bin/env python
import glob
from pyraf import iraf
import multiprocessing as mp
import sys
import os
import json
import shutil
import platform
from packaging import version
if version.parse(platform.python_version()) < version.parse("3.0.0"):
    print("Using raw_input")
else:
    raw_input = input

basepath = os.path.dirname(sys.argv[0])
lacos_im = os.path.join(basepath, '../iraf_tasks/lacos_im.cl')

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    gain = 1.41
    readnoise = 4.64
elif teles == "LJT":
    print("Settings for LJT will be used.")
    gain = 3.20
    readnoise = 13.5
else:
    print("Error detected.")

objall = []
with open('objall.list') as file:
    for line in file:
        line = line.strip('.fits\n')
        objall.append(line)


iraf.task(lacos_im=lacos_im)
print('Loading IRAF packages ...')
iraf.stsdas()
iraf.lacos_im.gain = gain
iraf.lacos_im.readn = readnoise

iraf.images()
iraf.images.imutil()
iraf.images.imutil.imheader(images="@objall.list")
stdspecs = str(raw_input("Enter filenames of all standard star spectra, \n\
separated by comma (','): "))
stdlist = [x.strip('.fits') for x in stdspecs.split(',')]
stdlist = [x.strip() for x in stdlist]
for item in stdlist:
    objall.remove(item)

def rmcrimg(obj):
    iraf.lacos_im(input='f' + obj, output='crf' + obj,
                  outmask='mask' + obj,
                  niter=4,
                  verbose='No')
    print(obj + ' finished.')


# pool = mp.Pool(mp.cpu_count())
# pool = mp.Pool(2)
# pool.map(rmcrimg, objall)
for obj in objall:
    rmcrimg(obj)

if len(stdlist)>0:
    for item in stdlist:
        if teles == "XLT":
            shutil.copy('f'+item+'.fit', 
                        'crf'+item+'.fits')
        elif teles == "LJT":
            shutil.copy('f'+item+'.fits', 
                        'crf'+item+'.fits')