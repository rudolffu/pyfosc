#!/usr/bin/env python
import glob
import re
from pyraf import iraf
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.coordinates import Angle
from astropy.time import Time
import pandas as pd
import astropy.units as u


def MyFilen(string):
    return '0' + string[8:11] + '.fit'


def NewUtc(fitsn):
    timestr = fits.getheader(fitsn)['DATE-OBS']
    return timestr.replace('T', ' ')


def NewSt(utctime):
    t = Time(utctime, scale='utc', location=xinglong)
    t.delta_ut1_utc = 0.
    return t.sidereal_time('apparent').value


def UtcHours(utctime):
    hangle = Angle(utctime.split(' ')[1] + 'hours')
    utchours = hangle.hour
    return utchours


xinglong = EarthLocation.of_site('beijing xinglong observatory')

gain = 1.41
ron = 4.64
txtlist = glob.glob('*.txt')
print('Possible log file in this folder:\n' + ", ".join(p for p in txtlist))
logfile = str(raw_input("Enter filename of the log file: "))
olog = pd.read_fwf(logfile, sep='\s+', skiprows=10)
olog['filen'] = map(MyFilen, olog['Files'])
narows = olog[olog.isna().any(axis=1)]
print(narows[['Obj', 'filen', 'Exp.(s)']])
lampspec = str(raw_input("Enter filename(s) of cal lamp\
spectrum(a) with longer exposures (separated by ',' for more than 1): "))
lamplist = re.split(',|,\s+', lampspec)
for img in lamplist:
    iraf.hedit.add = True
    iraf.hedit.addonly = True
    iraf.hedit.delete = False
    iraf.hedit.update = True
    iraf.hedit.verify = False
    iraf.hedit(img, fields='object', value='fear')

otbl = olog.dropna(subset=['R.A.', 'Dec.']).reset_index(drop=True)
otbl['hmsdms'] = otbl['R.A.'] + otbl['Dec.']
coord = SkyCoord(otbl['hmsdms'].values, frame='icrs',
                 unit=(u.hourangle, u.deg))
otbl['rah'] = coord.ra.to(u.hourangle).value
otbl['decdeg'] = coord.dec.to(u.degree).value
otbl['utctime'] = map(NewUtc, otbl['filen'])
otbl['st'] = map(NewSt, otbl['utctime'])
otbl['ut'] = map(UtcHours, otbl['utctime'])

iraf.twodspec()
iraf.longslit()
allimg = glob.glob('*.fit*')
for img in allimg:
    iraf.hedit.add = True
    iraf.hedit.addonly = True
    iraf.hedit.delete = False
    iraf.hedit.update = True
    iraf.hedit.verify = False
    iraf.hedit(img, fields='gain', value=gain)
    iraf.hedit(img, fields='rdnoise', value=ron)
    iraf.hedit(img, fields='epoch', value=2000.0)

for i in range(len(otbl)):
    iraf.hedit.add = True
    iraf.hedit.addonly = True
    iraf.hedit.delete = False
    iraf.hedit.update = True
    iraf.hedit.verify = False
    iraf.hedit(images=otbl['filen'][i], fields='ra', value=otbl['rah'][i])
    iraf.hedit(images=otbl['filen'][i], fields='dec', value=otbl['decdeg'][i])
    iraf.hedit(images=otbl['filen'][i], fields='object', value=otbl['Obj'][i])
    iraf.hedit(images=otbl['filen'][i], fields='st', value=otbl['st'][i])
    iraf.hedit(images=otbl['filen'][i], fields='ut', value=otbl['ut'][i])
    iraf.longslit.setairmass(images=otbl['filen'][i],
                             observatory='bao',
                             ra='ra',
                             dec='dec',
                             equinox='epoch',
                             st='st',
                             ut='ut',
                             date='date-obs')
