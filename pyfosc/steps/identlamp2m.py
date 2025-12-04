#!/usr/bin/env python
import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table, vstack
from pyraf import iraf


def headertable(filename: str, telescope: str | None = None) -> Table:
    """Extract key header fields into a table for sorting."""
    basename = os.path.basename(filename)
    hdr = fits.getheader(filename)
    arr = np.asarray(hdr.cards).T
    tab = Table(data=arr[1], names=arr[0])
    tab.add_column(basename, name="FILENAME")
    if telescope == "XLT":
        fields = [
            "OBJECT",
            "RA",
            "DEC",
            "IMAGETYP",
            "DATE-OBS",
            "EXPTIME",
            "OBSTYPE",
            "FILTER",
            "FILENAME",
            "TELESCOP",
            "INSTRUME",
            "JD",
            "AIRMASS",
            "OBSERVAT",
        ]
    elif telescope == "LJT":
        fields = [
            "OBJECT",
            "CAT-RA",
            "CAT-DEC",
            "FILTER1",
            "FILTER3",
            "DATE-OBS",
            "EXPTIME",
            "OBSTYPE",
            "FILTER",
            "FILENAME",
            "TELESCOP",
            "INSTRUME",
            "MJD-OBS",
            "AIRMASS",
            "SITE",
            "GAIN",
            "RDNOISE",
        ]
    elif telescope == "HCT":
        fields = [
            "OBJECT",
            "RA",
            "DEC",
            "IMAGETYP",
            "DATE-OBS",
            "EXPTIME",
            "GRISM",
            "FILTER",
            "FILENAME",
            "TELESCOP",
            "INSTRUME",
            "OBSERVAT",
        ]
    else:
        fields = tab.colnames
    return tab[fields]


with open("myfosc.json") as file:
    settings = json.loads(file.read())
teles = settings["mysettings"]["telescope"]
if teles == "XLT":
    print("Settings for XLT will be used.")
    linelists = "linelists$fear.dat"
elif teles == "LJT":
    print("Settings for LJT will be used.")
    linelists = "linelists$henear.dat"
elif teles == "HCT":
    print("Settings for HCT will be used.")
    linelists = "linelists$fear.dat"
else:
    raise RuntimeError("Unknown telescope; please set it in myfosc.json.")

flist = glob.glob("af*.fit*")
if not flist:
    raise RuntimeError("No lamp files matched pattern 'af*.fit*'.")

tablist = [headertable(item, teles) for item in flist]
tb = vstack(tablist)
tb.sort(["DATE-OBS"])
df = tb.to_pandas()

if teles == "LJT":
    idx_fear = df.FILTER.str.contains("lamp_fe_argon")
    idx_hene = df.FILTER.str.contains("lamp_neon_helium")
    print("Fe-Ar lamp spectra:")
    print(df.loc[idx_fear, "FILENAME"])
    iraf.onedspec()
    iraf.onedspec.identify.unlearn()
    iraf.onedspec.identify.fwidth = 10
    iraf.onedspec.identify.coordli = "linelists$fear.dat"
    file_str = df.loc[idx_fear, "FILENAME"].str.cat(sep=",")
    iraf.onedspec.identify(images=file_str)
    print("He-Ne lamp spectra:")
    print(df.loc[idx_hene, "FILENAME"])
    iraf.onedspec()
    iraf.onedspec.identify.unlearn()
    iraf.onedspec.identify.fwidth = 10
    iraf.onedspec.identify.coordli = "linelists$henear.dat"
    file_str = df.loc[idx_hene, "FILENAME"].str.cat(sep=",")
    iraf.onedspec.identify(images=file_str)
else:
    iraf.onedspec()
    iraf.onedspec.identify.unlearn()
    iraf.onedspec.identify.fwidth = 10
    iraf.onedspec.identify.coordli = linelists
    iraf.onedspec.identify(images="af*.fit*")
