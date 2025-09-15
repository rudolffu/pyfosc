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
side = settings['mysettings']['side']
if teles == "XLT":
    print("Settings for XLT will be used.")
elif teles == "LJT":
    print("Settings for LJT will be used.")
elif teles == "HCT":
    print("Settings for HCT will be used.")
elif teles == "P200":
    print("Settings for P200 will be used.")
else:
    print("Error detected.")


def headertable(filename, telescope=None):
    basename = os.path.basename(filename)
    hdr = fits.getheader(filename)
    if 'COMMENT' in hdr:
        del hdr['COMMENT']
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
    elif telescope=='P200':
        tab.add_column(basename, name='FILENAME')
        fields = ['OBJECT','RA','DEC','IMGTYPE',
                  'UTSHUT','EXPTIME',
                  'FILENAME', 'GRATING', 
                  'GAIN','RON']
    return tab[fields]


def get_file_ext(filename):
    return os.path.splitext(filename)[1]


def rename_raw(df, telescope=None, side=side):
    if telescope=='XLT':
        num = df.FILENAME.str[8:12]
        ext = df.FILENAME.apply(get_file_ext)
    elif telescope=='LJT':
        num = df.FILENAME.str[23:27]
        ext = df.FILENAME.apply(get_file_ext)
    elif telescope=='HCT':
        num = df.FILENAME.str[3:8]
        ext = '.fits'
    elif telescope=='P200':
        ext = '.fits'
        num = df.FILENAME.str[-9:-5]
    df.OBJECT = df.OBJECT.str.replace('N/A','NULL')
    df.OBJECT = df.OBJECT.str.replace(' ','')
    objname = '_' + df.OBJECT.str[:10]
    objname = objname.str.replace('_NULL', '')
    if telescope=='HCT':
        obstype = df.IMAGETYP.str.lower()
    elif telescope=='P200':
        obstype = side[0]+df.IMGTYPE.str.lower()
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

