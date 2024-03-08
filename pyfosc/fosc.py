#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.table import vstack,Table
from astropy.time import Time
from astropy.stats import sigma_clip
import astropy.units as u
import pandas as pd
import numpy as np
import os
import shutil
import multiprocessing as mp
from pathlib import Path
import json
import datetime
import warnings
import glob
import ccdproc as ccdp
from ccdproc import ImageFileCollection
from matplotlib.colors import LogNorm


class BasicCCDMixin:
    def plot_image(self, 
                   log_scale=True, 
                   use_wcs=False,
                   cmap='gray',
                   vmin=None, vmax=None):
        data = self.data
        if vmin is None and vmax is None:
            vmin = np.percentile(data, 1)
            vmax = np.percentile(data, 99)
        if log_scale is True:
            norm = LogNorm(vmin, vmax)
        else:
            norm = None
        if use_wcs is True and self.wcs is not None:
            projection = self.wcs
            
        else:
            projection = None
        ax = plt.subplot(projection=projection)
        im = ax.imshow(data, cmap=cmap, norm=norm, 
                       origin='lower')
        plt.colorbar(im, ax=ax, label=self.unit)
    
    

class SpecImage(BasicCCDMixin, CCDData):
    def __init__(self, ccddata, disp_axis=None,
                 *args, **kwd):
        super().__init__(ccddata, *args, **kwd)
        self._disp_axis = disp_axis
        
    @classmethod
    def read(cls, filename, hdu=0, unit=None,
             hdu_uncertainty="UNCERT",
             hdu_mask="MASK",
             hdu_flags=None,
             key_uncertainty_type="UTYPE",
             hdu_psf="PSFIMAGE",
             **kwd):
        ccddata = super().read(filename, 
                                hdu, unit, 
                                hdu_uncertainty, 
                                hdu_mask, hdu_flags, 
                                key_uncertainty_type, 
                                hdu_psf, **kwd)
        return cls(ccddata)
        
    @property
    def disp_axis(self):
        return self._disp_axis
    
    @disp_axis.setter
    def disp_axis(self, value):
        # if self.shape[0] != self.shape[1]:
        #     long_axis = np.argmax(self.shape)
        # if value is None:
        #     try:
        #         value = long_axis
        #     except:
        #         value = 0
        #         print('No dispersion axis specified, defaulting to 0')
        self._disp_axis = int(value)
