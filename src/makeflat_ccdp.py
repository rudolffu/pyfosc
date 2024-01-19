#!/usr/bin/env python
import os
import json
from ccdproc import ImageFileCollection,Combiner
from fnmatch import fnmatch
from astropy import units as u

with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']

flat_list=[]
with open('flat.list') as file:
    for line in file:
        line = line.strip('\n')
        flat_list.append(line)

flats = ImageFileCollection('./',filenames=flat_list)
c = Combiner(flats.ccds(ccd_kwargs={'unit':u.adu}))
c.sigma_clipping() #dev_func
avg_combined = c.average_combine()
if teles == "XLT" :
    avg_combined.write('Flat.fit', overwrite=True)
elif teles in ("LJT","HCT") :
    avg_combined.write('Flat.fits', overwrite=True)

print('--- DONE ---')
