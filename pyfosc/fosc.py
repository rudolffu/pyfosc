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


class BasicCCDMixin():
    """A mixin class for CCDData objects"""
    def plot_image(self, log_scale=True, 
                   use_wcs=False, cmap='gray',
                   vmin=None, vmax=None):
        """
        Plot the image data with optional log scale 
            and WCS projection.
        Parameters
        ----------
        log_scale : bool, optional
            If True, use log scale for the image.
        use_wcs : bool, optional
            If True, use the WCS projection for the image.
        cmap : str, optional
            The colormap for the image.
        vmin : float, optional
            The minimum value for the image.
        vmax : float, optional
            The maximum value for the image.
        
        Returns
        -------
        None
            
        """    
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
    """A class for 2d spectral images"""
    def __init__(self, ccddata, disp_axis=None,
                 *args, **kwd):
        """
        Parameters
        ----------
        ccddata : CCDData
            The CCDData object to be used.
        disp_axis : int, optional
            The dispersion axis of the image.
        *args : list
            Additional positional arguments.
        **kwd : dict
            Additional keyword arguments.
        """
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
        """
        Read a spectral image from a file and 
            return a SpecImage object.
        Parameters
        ----------
        filename : str
            The filename of the image.
        hdu : int, optional
            The HDU to read from the file.
        unit : astropy.unit or str, optional
            The unit of the image data.
        """
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
