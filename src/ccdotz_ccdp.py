#!/usr/bin/env python
# Load Python standard modules
from astropy.io import fits
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection,subtract_overscan,trim_image,subtract_bias
import json
from glob import glob

with open('myfosc.json') as file:
    settings = json.loads(file.read())

teles = settings['mysettings']['telescope']
Grism = settings['mysettings']['Grism']

if teles == "XLT":
    print("Settings for XLT will be used.")
    overscan = False
    if Grism == "G4":
        trimsec = '[51:1750,681:1350]'
    elif Grism == "G7":
        trimsec = '[51:2000,681:1350]'

elif teles == "LJT":
    print("Settings for LJT will be used.")
    overscan = True
    flist = glob("*fits")
    for f in flist:
        hdu = fits.open(f,mode='update')
        hdr = hdu[0].header
        try:
            del hdr['CCDSEC']
        except:
            pass
        hdu.flush()
        hdu.close()
    if Grism == "G3":
        biassec = '[2100:2148,2301:4130]'
        trimsec = '[651:1350,2301:4130]'
        # If you want more trimming on the blue side, use:
        # biassec = '[10:40,2491:4130]'
        # trimsec = '[751:1450,2491:4130]'
    elif Grism == "G8":
        biassec = '[2100:2148,1136:4250]'
        trimsec = '[751:1450,1136:4250]'
    elif Grism == "G10":
        biassec = '[2100:2148,2591:3340]'
        trimsec = '[751:1450,2591:3340]'
    elif Grism == "G14":
        biassec = '[2100:2148,1901:4200]'
        trimsec = '[751:1450,1901:4200]'

elif teles == "HCT":
    print("Settings for HCT will be used.")
    overscan = False
    trimsec = '[26:250,101:2800]'

else:
    print("Error detected.")

flatnall_list=[]
with open('flatnall.list') as file:
    for line in file:
        line = line.strip('\n')
        flatnall_list.append(line)

flatnalls = ImageFileCollection('./',filenames=flatnall_list)
tb_flatnalls = flatnalls.summary

zero = CCDData.read('./Zero.fit', unit='adu')
for file_name in tb_flatnalls['file'].tolist():
    ccd = CCDData.read(file_name, unit='adu')
    if overscan:
        ccd = subtract_overscan(ccd,fits_section=biassec,overscan_axis=1)
    ccd = subtract_bias(ccd,zero)
    ccd = trim_image(ccd,fits_section=trimsec)
    ccd.write(file_name,overwrite=True)

print('--- DONE ---')
