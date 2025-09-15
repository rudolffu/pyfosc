#!/usr/bin/env python
# coding: utf-8
from astropy.io import fits
import pandas as pd
import numpy as np
import astropy.units as u
import glob
from pathlib import Path
import json


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
elif teles=='HCT':
    list_bias = glob.glob('bias*.fits')
    list_obj = glob.glob('object*.fits')
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

print('Wrote zero.list, flat.list, objall.list, lampall.list, flatnall.list, specall.list')

