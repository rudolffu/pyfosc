#!/usr/bin/env python
# Load Python standard modules
import os
import json
from ccdproc import ImageFileCollection,Combiner
from fnmatch import fnmatch
from astropy import units as u

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']

zero_list=[]
with open('zero.list') as file:
    for line in file:
        line = line.strip('\n')
        zero_list.append(line)

zeros = ImageFileCollection('./',filenames=zero_list)

c = Combiner(zeros.ccds(ccd_kwargs={'unit':u.adu}))
c.clip_extrema(nlow=0, nhigh=1)
avg_combined = c.average_combine()
avg_combined.mask = None
avg_combined.uncertainty = None
if teles == "XLT" :
    avg_combined.write('Zero.fit', overwrite=True)
elif teles in ("LJT","HCT") :
    avg_combined.write('Zero.fits', overwrite=True)

print('--- DONE ---')

