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
    """
    Added try-except due to changes in fits headers of BFOSC data.
    """
    try:
        timestr1 = fits.getheader(fitsn)['DATE-OBS']
        timestr2 = fits.getheader(fitsn)['TIME-OBS']
        timestr = timestr1 + " "+ timestr2 
    except:
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
fitslist = glob.glob('*.fit*')
fitslist = pd.DataFrame(fitslist, columns=['filen'])
print('Possible log file in this folder:\n' + ", ".join(p for p in txtlist))
logfile = str(raw_input("Enter filename of the log file: "))
olog = pd.read_fwf(logfile, sep='\s+', skiprows=10)
olog['filen'] = map(MyFilen, olog['Files'])
narows = olog[olog.isna().any(axis=1)]
print(narows[['Obj', 'filen', 'Exp.(s)']])
lampspec = str(raw_input("Enter filename(s) of cal lamp \
spectrum(a) with longer exposures (separated by ',' for more than 1): "))
lamplist = re.split(',|,\s+', lampspec)
for img in lamplist:
    try:
        iraf.hedit.add = True
        iraf.hedit.addonly = True
        iraf.hedit.delete = False
        iraf.hedit.update = True
        iraf.hedit.verify = False
        iraf.hedit(img, fields='OBJECT', value='fear')
    except:
        pass

otbl = olog.dropna(subset=['R.A.', 'Dec.']).reset_index(drop=True)
otbl['hmsdms'] = otbl['R.A.'] + otbl['Dec.']
coord = SkyCoord(otbl['hmsdms'].values, frame='icrs',
                 unit=(u.hourangle, u.deg))
otbl['rah'] = coord.ra.to(u.hourangle).value
otbl['radeg'] = coord.ra.to(u.degree).value
otbl['decdeg'] = coord.dec.to(u.degree).value
notbl = otbl.merge(fitslist, on="filen")
notbl['utctime'] = map(NewUtc, notbl['filen'])
notbl['st'] = map(NewSt, notbl['utctime'])
notbl['ut'] = map(UtcHours, notbl['utctime'])
notbl['utdatetime'] = notbl.utctime.str.replace(' ', 'T')


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

for i in range(len(notbl)):
    try:
        iraf.hedit.add = True
        iraf.hedit.addonly = True
        iraf.hedit.delete = False
        iraf.hedit.update = True
        iraf.hedit.verify = False
        iraf.hedit(images=notbl['filen'][i],
                   fields='object', value=notbl['Obj'][i])
        iraf.hedit(images=notbl['filen'][i], fields='st', value=notbl['st'][i])
        iraf.hedit(images=notbl['filen'][i], fields='ut', value=notbl['ut'][i])
        iraf.hedit(images=notbl['filen'][i], 
                   fields='OBSERVATORY', 
                   value="BAO")
        iraf.hedit(images=notbl['filen'][i], fields='UTDTIME', value=notbl['utdatetime'][i])
        iraf.longslit.setairmass(images=notbl['filen'][i],
                                 observatory='BAO',
                                 ra='ra',
                                 dec='dec',
                                 equinox='epoch',
                                 st='st',
                                 ut='ut',
                                 date='UTDTIME')
        iraf.hedit(images=notbl['filen'][i],
                   fields='ra', value=notbl['radeg'][i])
        iraf.hedit(images=notbl['filen'][i],
                   fields='dec', value=notbl['decdeg'][i])
    except:
        pass
