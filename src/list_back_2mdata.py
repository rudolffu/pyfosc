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
import json


with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
if teles == "XLT":
    print("Settings for XLT will be used.")
    # midname = 'XL'
elif teles == "LJT":
    print("Settings for LJT will be used.")
    # midname = 'LJ'
else:
    print("Error detected.")


def headertable(filename, telescope=None):
    basename = os.path.basename(filename)
    hdr = fits.getheader(filename)
    arr = np.asarray(hdr.cards).T
    tab = Table(data=arr[1], names=arr[0])
    tab.add_column(basename, name='filename')
    if telescope=='XLT':
        fields = ['OBJECT','RA','DEC','IMAGETYP',
                  'DATE-OBS','EXPTIME','OBSTYPE',
                  'FILTER','filename','TELESCOP',
                  'INSTRUME','JD','AIRMASS','OBSERVAT']
    elif telescope=='LJT':
        fields = ['OBJECT','CAT-RA','CAT-DEC',
                  'FILTER1','FILTER3',
                  'DATE-OBS','EXPTIME','OBSTYPE',
                  'FILTER','filename','TELESCOP',
                  'INSTRUME','MJD-OBS','AIRMASS',
                  'SITE','GAIN','RDNOISE']
    return tab[fields]


def get_file_ext(filename):
    return os.path.splitext(filename)[1]


def rename_raw(df, telescope=None):
    if telescope=='XLT':
        num = df.filename.str[8:12]
    elif telescope=='LJT':
        num = df.filename.str[23:27]
    df.OBJECT = df.OBJECT.str.replace('N/A','NULL')
    df.OBJECT = df.OBJECT.str.replace(' ','')
    objname = '_' + df.OBJECT.str[:10]
    objname = objname.str.replace('_NULL', '')
    obstype = df.OBSTYPE.str.lower()
    ext = df.filename.apply(get_file_ext)
    newname = obstype + num + objname + ext
    return newname


def backup_raw(img, newname):
    shutil.copy(img, newname)
    shutil.move(img, raw_dir + '/' + img)


def mod_header(filename, field, value):
    with fits.open(filename, 'update') as f:
        for hdu in f:
            hdu.header[field] = value

if teles=='XLT':
    gain = 2.2                                                        
    ron = 7.8
    flist = glob.glob('*fit')
    tablist = []
    for item in flist:
        mod_header(item, field='gain', value=gain)
        mod_header(item, field='rdnoise', value=ron)
        tablist.append(headertable(item, teles))
elif teles=='LJT':
    flist = glob.glob('lj*fits')
    tablist = []
    for item in flist:
        tablist.append(headertable(item, teles))

tb = vstack(tablist)
tb.sort(['DATE-OBS'])
df = tb.to_pandas()
# df.query('OBSTYPE=="BIAS"')


if teles=='LJT':
    tb = tb[tb['FILTER3']!='grism14']
    tb = tb[tb['FILTER3']!='grism8']
    tb = tb[tb['FILTER3']!='grism10']
    tb = tb[tb['FILTER1']!='lslit5_05']
    tb = tb[tb['FILTER1']!='lslit1_81']
    df = tb.to_pandas()
    df.loc[df.FILTER.str.contains('lamp_neon_helium'),'OBSTYPE'] = 'CAL'
    df.OBSTYPE = df.OBSTYPE.str.replace('EXPERIMENTAL', 'EXPOSE')
    df.loc[
        (df.FILTER.str.contains('grism')) & (df.FILTER.str.contains('slit'))
        & (~df.FILTER.str.contains('lamp_neon_helium')), 'OBSTYPE'] = 'SCI'
    df.loc[df.FILTER.str.contains('lamp_halogen'), 'OBSTYPE'] = 'LAMPFLAT'


imglist = df.filename
newname = rename_raw(df, teles)
df['newname'] = newname
raw_dir = './raw'
Path(raw_dir).mkdir(exist_ok=True)

df.to_csv('imglist.csv', index=False)

if __name__ == '__main__':
    try:
        pool = mp.Pool(mp.cpu_count())
        pool.starmap(backup_raw, zip(imglist, newname))
        print('Copying raw files finished.')
    except FileNotFoundError:
        print('File not found. The files might '+
              'have been backed-up, otherwise a wrong location is provided.')

if teles=='XLT':
    try:
        list_std = glob.glob('speclfluxref*.fit')
    except:
        list_std = []
    list_bias = glob.glob('bias*.fit')
    list_obj = glob.glob('specltar*.fit')
    list_flat = glob.glob('speclflat*.fit')
    list_lamp = glob.glob('specllamp*.fit')
    list_obj.extend(list_std)
    list_flatnall = list_flat.copy()
    list_flatnall.extend(list_lamp)
    list_flatnall.extend(list_obj)
    list_specall = list_obj.copy()
    list_specall.extend(list_lamp)
elif teles=='LJT':
    list_bias = glob.glob('bias*.fits')
    list_obj = glob.glob('sci*.fits')
    list_flat = glob.glob('lampflat*.fits')
    list_lamp = glob.glob('cal*.fits')
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