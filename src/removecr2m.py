#!/usr/bin/env python
import glob
from pyraf import iraf
import multiprocessing as mp
import sys
import os
import json

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


def rmcrimg(obj):
    iraf.lacos_im(input='f' + obj, output='crf' + obj,
                  outmask='mask' + obj,
                  verbose='No')
    print(obj + ' finished.')


pool = mp.Pool(mp.cpu_count())
pool.map(rmcrimg, objall)
