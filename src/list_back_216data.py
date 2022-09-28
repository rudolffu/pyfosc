#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack,Table
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import glob
import os
from pathlib import Path
import shutil
import multiprocess as mp


def headertable(filename):
    basename = os.path.basename(filename)
    hdr = fits.getheader(filename)
    arr = np.asarray(hdr.cards).T
    tab = Table(data=arr[1], names=arr[0])
    tab.add_column(basename, name='filename')
    # fields = ['OBJECT','RA','DEC','IMAGETYP',
    #           'DATE-OBS','EXPTIME','OBSTYPE',
    #           'FILTER','filename','TELESCOP',
    #           'INSTRUME','JD','AIRMASS',
    #           'OBSERVAT','GAIN','RDNOISE']
    fields = ['OBJECT','RA','DEC','IMAGETYP',
              'DATE-OBS','EXPTIME','OBSTYPE',
              'FILTER','filename','TELESCOP',
              'INSTRUME','JD','AIRMASS','OBSERVAT']
    return tab[fields]


def get_file_ext(filename):
    return os.path.splitext(filename)[1]


def rename_raw(df):
    num = df.filename.str[8:12]
    obstype = df.OBSTYPE.str.lower()
    ext = df.filename.apply(get_file_ext)
    newname = obstype + num + ext
    return newname


def mod_header(filename, field, value):
    with fits.open(filename, 'update') as f:
        for hdu in f:
            hdu.header[field] = value


gain    = 2.2                                                         
ron = 7.8
flist = glob.glob('*fit')
tablist = []
for item in flist:
    mod_header(item, field='gain', value=gain)
    mod_header(item, field='rdnoise', value=ron)
    tablist.append(headertable(item))
tb = vstack(tablist)
tb.sort(['DATE-OBS'])
df = tb.to_pandas()
# df.query('OBSTYPE=="BIAS"')

imglist = df.filename
newname = rename_raw(df)
df['newname'] = newname
raw_dir = './raw'
Path(raw_dir).mkdir(exist_ok=True)

def backup_raw(img, newname):
    shutil.copy(img, newname)
    shutil.move(img, raw_dir + '/' + img)


if __name__ == '__main__':
    try:
        pool = mp.Pool(mp.cpu_count())
        pool.starmap(backup_raw, zip(imglist, newname))
        print('Copying raw files finished.')
    except FileNotFoundError:
        print('File not found. The files might '+
              'have been backed-up, otherwise a wrong location is provided.')


df.to_csv('imglist.csv', index=False)

list_bias = glob.glob('bias*.fit')
list_obj = glob.glob('specltar*.fit')
list_flat = glob.glob('speclflat*.fit')
list_lamp = glob.glob('specllamp*.fit')
list_flatnall = list_flat.copy()
list_flatnall.extend(list_lamp)
list_flatnall.extend(list_obj)
list_specall = list_obj.copy()
list_specall.extend(list_lamp)

with open('zero.list', 'w') as f:
    for line in list_bias:
        f.write(f"{line}\n")

with open('flat.list', 'w') as f:
    for line in list_flat:
        f.write(f"{line}\n")

with open('objall.list', 'w') as f:
    for line in list_obj:
        f.write(f"{line}\n")

with open('lampall.list', 'w') as f:
    for line in list_lamp:
        f.write(f"{line}\n")

with open('flatnall.list', 'w') as f:
    for line in list_flatnall:
        f.write(f"{line}\n")

with open('specall.list', 'w') as f:
    for line in list_specall:
        f.write(f"{line}\n")