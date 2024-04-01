#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack,Table
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
import glob
import os
from pathlib import Path
import shutil
import multiprocessing as mp
import json
import datetime
from astropy.time import Time
import warnings


with open('myfosc.json') as file:
    settings = json.loads(file.read())
teles = settings['mysettings']['telescope']
Grism = settings['mysettings']['Grism']
slit = settings['mysettings']['slit']
if teles == "XLT":
    print("Settings for XLT will be used.")
elif teles == "LJT":
    print("Settings for LJT will be used.")
elif teles == "HCT":
    print("Settings for HCT will be used.")
else:
    print("Error detected.")


def headertable(filename, telescope=None):
    basename = os.path.basename(filename)
    hdr = fits.getheader(filename)
    arr = np.asarray(hdr.cards).T
    tab = Table(data=arr[1], names=arr[0])
    if telescope=='XLT':
        tab.add_column(basename, name='FILENAME')
        fields = ['OBJECT','RA','DEC','IMAGETYP',
                  'DATE-OBS','EXPTIME','OBSTYPE',
                  'FILTER','FILENAME','TELESCOP',
                  'INSTRUME','JD','AIRMASS','OBSERVAT']
    elif telescope=='LJT':
        tab.add_column(basename, name='FILENAME')
        fields = ['OBJECT','CAT-RA','CAT-DEC',
                  'FILTER1','FILTER3',
                  'DATE-OBS','EXPTIME','OBSTYPE',
                  'FILTER','FILENAME','TELESCOP',
                  'INSTRUME','MJD-OBS','AIRMASS',
                  'SITE','GAIN','RDNOISE']
    elif telescope=='HCT':
        fields = ['OBJECT','RA','DEC','IMAGETYP',
                  'DATE-OBS','EXPTIME','GRISM',
                  'FILTER','FILENAME','TELESCOP',
                  'INSTRUME','OBSERVAT']
    return tab[fields]


def get_file_ext(filename):
    return os.path.splitext(filename)[1]


def rename_raw(df, telescope=None):
    if telescope=='XLT':
        num = df.FILENAME.str[8:12]
        ext = df.FILENAME.apply(get_file_ext)
    elif telescope=='LJT':
        num = df.FILENAME.str[23:27]
        ext = df.FILENAME.apply(get_file_ext)
    elif telescope=='HCT':
        num = df.FILENAME.str[3:8]
        ext = '.fits'
    df.OBJECT = df.OBJECT.str.replace('N/A','NULL')
    df.OBJECT = df.OBJECT.str.replace(' ','')
    objname = '_' + df.OBJECT.str[:10]
    objname = objname.str.replace('_NULL', '')
    if telescope=='HCT':
        obstype = df.IMAGETYP.str.lower()
    else:
        obstype = df.OBSTYPE.str.lower()
    newname = obstype + num + objname + ext
    return newname


def backup_raw(img, newname):
    raw_dir = './raw'
    data_dir = './data'
    Path(raw_dir).mkdir(exist_ok=True)
    Path(data_dir).mkdir(exist_ok=True)
    shutil.copy(img, data_dir + '/' + newname)
    shutil.move(img, raw_dir + '/' + img)


def mod_header(filename, field, value):
    with fits.open(filename, 'update') as f:
        for hdu in f:
            hdu.header[field] = value
        f.flush()

def date2isostr_HCT(filename):
    """
    Convert date and time stamps to UTC datetime string in ISO format (YYYY-MM-DDTHH:MM:SS).
    Parameter:
        filename : path-like or file-like 
            FITS file to get header from. 
    Returns:
        isostr : str
            UTC datetime string in ISO format.
    """
    hdr = fits.getheader(filename)
    date = hdr['DATE-OBS']
    if 'T' in hdr['DATE-OBS']:
        warnings.warn('Field DATE-OBS already contains date and time in ISO format. \
The value of DATE-OBS is unchanged.')
        isostr = date
    else:
        hms_str = hdr.comments['TM_START'].split(' ')[0]
        hms_str =  hms_str.replace('/', ':')
        isostr = "{}T{}".format(date, hms_str)
    # tnew = datetime.datetime.fromisoformat(isostr)
    # tnew = Time(tnew) 
    return isostr


def isostr2sidereal(isostr, location):
    """
    Convert UTC datetime string in ISO format to sidereal time in a given location.
    Parameter:
        isostr : str
            UTC datetime string in ISO format.
        location : `~astropy.coordinates.EarthLocation` 
            Location on the Earth.
    Returns:
        sdrtime : float64 
            Value of the sidereal time.
    """    
    utcTime = Time(isostr, scale='utc', location=location)
    sdrtime = utcTime.sidereal_time('apparent').value
    return sdrtime


def isostr2UTChours(isostr, location):
    utcTime = Time(isostr, scale='utc', location=location)
    hangle = Angle(utcTime.isot.split('T')[1] + 'hours')
    utchours = hangle.hour
    return utchours

def getairmass(ra, dec, obstime, location):
    try:
        ra = float(ra)
        dec = float(dec)
        coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg,
                         frame='icrs')
    except:
        ra_str = ra
        dec_str = dec
        coord = SkyCoord('{} {}'.format(ra_str, dec_str),
                         unit=(u.hourangle, u.deg),
                         frame='icrs')
    co_altaz = coord.transform_to(AltAz(obstime=obstime,location=location))
    return co_altaz.secz.value


class HCTfits():
    def __init__(self, filename, site='IAO', fix=True) -> None:
        self.filename = filename
        self.location = EarthLocation.of_site(site)
        if fix==True:
            hdu = fits.open(filename, mode='update')
            hdu[0].verify('fix')
            hdu.flush()
            hdu.close()
        hdr = fits.getheader(self.filename)
        try:
            if ':' in hdr['RA'] and ':' in hdr['DEC']:
                coord = SkyCoord('{} {}'.format(hdr['RA'], hdr['DEC']),
                                 unit=(u.hourangle, u.deg), frame='icrs')
            else:
                ra = float(hdr['RA'])
                dec = float(hdr['DEC'])
                coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
            self.coord = coord
            self.ra = coord.ra.value
            self.dec = coord.dec.value
        except:
            self.coord = None
            self.ra = None
            self.dec = None
        self.hdr = hdr

    def fixhdu(self):
        filename = self.filename
        hdu = fits.open(filename, mode='update')
        hdu[0].verify('fix')
        hdu.flush()
        hdu.close()
    
    def update_field(self, field, value):
        filename = self.filename
        hdu = fits.open(filename, mode='update')
        hdu[0].header[field] = value
        hdu.flush()
        hdu.close()

    def update_dateobs(self):
        hdr = self.hdr
        date = hdr['DATE-OBS']
        if 'T' in hdr['DATE-OBS']:
            warnings.warn('Field DATE-OBS already contains date and time in ISO format. \
The value of DATE-OBS is unchanged.')
            isostr = date
        else:
            hms_str = hdr.comments['TM_START'].split(' ')[0]
            hms_str =  hms_str.replace('/', ':')
            isostr = "{}T{}".format(date, hms_str)
        self.update_field('DATE-OBS', isostr)
        exptime = hdr['EXPTIME']
        t_delta = datetime.timedelta(seconds=exptime/2)
        datemid = datetime.datetime.fromisoformat(isostr) + t_delta
        datemid = Time(datemid, location=self.location)
        self.datemid = datemid
        self.update_field('UTC-MID', datemid.isot)

    def update_airmass(self):
        co_altaz = self.coord.transform_to(AltAz(obstime=self.datemid,location=self.location))
        self.co_altaz = co_altaz
        self.airmass = co_altaz.secz.value
        self.update_field('AIRMASS', self.airmass)
    

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
elif teles=='HCT':
    gain = 1.22                                                        
    ron = 4.8
    print("Now printing files in the current directory......")
    print(os.listdir('./'))
    filename_prefix = str(input("Enter prefix of the HCT spectra, e.g.'af': "))
    flist = glob.glob(filename_prefix+'*')
    tablist = []
    for item in flist:
        hf = HCTfits(item)
        hf.update_field('gain', gain)
        hf.update_field('rdnoise', ron)
        if hf.coord is not None:
            hf.update_dateobs()
            hf.update_airmass()
        else:
            hf.update_field('RA', 0.0)
            hf.update_field('DEC', 0.0)
        if 'IMAGETYP' not in hf.hdr.keys():
            hf.update_field('IMAGETYP', 'undefined')
        tablist.append(headertable(item, teles))
tb = vstack(tablist)
tb.sort(['DATE-OBS'])
df = tb.to_pandas()
# df.query('OBSTYPE=="BIAS"')


grisms_ljt = {'G3': 'grism3', 
              'G8': 'grism8', 
              'G10': 'grism10', 
              'G14': 'grism14'} 
slits_ljt = {'slit1.8': 'lslit1_81', 
             'slit2.5': 'lslit2_51', 
             'slit5.0': 'lslit5_05'} 
if teles=='LJT':
    key_grism = grisms_ljt.pop(Grism)
    key_slit = slits_ljt.pop(slit)
    for grism_item in list(grisms_ljt.values()):
        tb = tb[tb['FILTER3']!=grism_item]
    for slit_item in list(slits_ljt.values()):
        tb = tb[tb['FILTER1']!=slit_item]
    df = tb.to_pandas()
    df.loc[df.FILTER.str.contains('lamp_neon_helium'),'OBSTYPE'] = 'CAL'
    df.loc[df.FILTER.str.contains('lamp_fe_argon'),'OBSTYPE'] = 'CAL'
    df.OBSTYPE = df.OBSTYPE.str.replace('EXPERIMENTAL', 'EXPOSE')
    df.loc[
        (df.FILTER.str.contains('grism')) & (df.FILTER.str.contains('slit'))
        & (~df.FILTER.str.contains('lamp_neon_helium'))
        & (~df.FILTER.str.contains('lamp_fe_argon')), 'OBSTYPE'] = 'SCI'
    df.loc[df.FILTER.str.contains('lamp_halogen'), 'OBSTYPE'] = 'LAMPFLAT'
elif teles=='HCT':
    df['OBJECT_low'] = df.OBJECT.str.lower()
    df.loc[(df.IMAGETYP.str.contains('lamp') 
            & (~df.OBJECT_low.str.contains('halogen'))),'IMAGETYP'] = 'CAL'
    df.loc[df.OBJECT_low.str.contains('halogen'), 'IMAGETYP'] = 'LAMPFLAT'


imglist = df.FILENAME
newname = rename_raw(df, teles)
df['newname'] = newname

df.to_csv('imglist.csv', index=False)

if __name__ == '__main__':
    try:
        # pool = mp.Pool(mp.cpu_count())
        # pool.starmap(backup_raw, zip(imglist, newname))
        for img, newname in zip(imglist, newname):
            backup_raw(img, newname)
        print('Copying raw files finished.')
    except FileNotFoundError:
        print('File not found. The files might '+
              'have been backed-up, otherwise a wrong location is provided.')

# if teles=='XLT':
#     try:
#         list_std = glob.glob('speclfluxref*.fit')
#     except:
#         list_std = []
#     list_bias = glob.glob('bias*.fit')
#     list_obj = glob.glob('specltar*.fit')
#     list_flat = glob.glob('speclflat*.fit')
#     list_lamp = glob.glob('specllamp*.fit')
#     list_obj.extend(list_std)
#     list_flatnall = list_flat.copy()
#     list_flatnall.extend(list_lamp)
#     list_flatnall.extend(list_obj)
#     list_specall = list_obj.copy()
#     list_specall.extend(list_lamp)
# elif teles=='LJT':
#     list_bias = glob.glob('bias*.fits')
#     list_obj = glob.glob('sci*.fits')
#     list_flat = glob.glob('lampflat*.fits')
#     list_lamp = glob.glob('cal*.fits')
#     list_flatnall = list_flat.copy()
#     list_flatnall.extend(list_lamp)
#     list_flatnall.extend(list_obj)
#     list_specall = list_obj.copy()
#     list_specall.extend(list_lamp)
# elif teles=='HCT':
#     list_bias = glob.glob('bias*.fits')
#     list_obj = glob.glob('object*.fits')
#     list_flat = glob.glob('lampflat*.fits')
#     list_lamp = glob.glob('cal*.fits')
#     list_flatnall = list_flat.copy()
#     list_flatnall.extend(list_lamp)
#     list_flatnall.extend(list_obj)
#     list_specall = list_obj.copy()
#     list_specall.extend(list_lamp)


# with open('zero.list', 'w') as f:
#     for line in list_bias:
#         f.write(f"{line}\n")

# with open('flat.list', 'w') as f:
#     for line in list_flat:
#         f.write(f"{line}\n")

# with open('objall.list', 'w') as f:
#     for line in list_obj:
#         f.write(f"{line}\n")

# with open('lampall.list', 'w') as f:
#     for line in list_lamp:
#         f.write(f"{line}\n")

# with open('flatnall.list', 'w') as f:
#     for line in list_flatnall:
#         f.write(f"{line}\n")

# with open('specall.list', 'w') as f:
#     for line in list_specall:
#         f.write(f"{line}\n")