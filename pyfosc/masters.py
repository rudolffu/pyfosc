#!/usr/bin/env python
# coding: utf-8
import numpy as np
import ccdproc as ccdp
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection
from pathlib import Path
import os
from .fosc import SpecImage




class MasterBias:
    """Class to create and handle master bias frames."""
    def __init__(self, 
                 work_dir=None,
                 list_bias=None, 
                 master_file=None):
        """
        Parameters
        ----------
        work_dir : str, optional
            The working directory where the data is located. The default is None.
        list_bias : list, optional
            List of bias frames. The default is None.
        master_file : str, optional
            The name of the master bias file. The default is None.        
        """
        if master_file is not None:
            self.master_bias = SpecImage.read(master_file)
            return 
        else:
            self.master_bias = None
        if work_dir is None:
            location = os.getcwd()
        else:
            location = Path(work_dir) / 'data'
        if list_bias is not None:
            bias_collection = ImageFileCollection(
                location=location,
                filenames=list_bias)
        self.work_dir = work_dir
        self.list_bias = list_bias
        self.bias_collection = bias_collection
        self.master_file = master_file
        
    def build(self, save=False, output_file='master_bias.fits'):
        """
        Create the master bias frame.
        
        Parameters
        ----------
        save : bool, optional
            Save the master bias frame. The default is False.
        output_file : str, optional
            The name of the output file. The default is 'master_bias.fits'.
            
        Returns
        -------
        master_bias : SpecImage(ccdproc.CCDData)
            The master bias frame.
        """
        if self.master_bias is not None:
            print('Master bias already exists. Returning the master bias.')
            return self.master_bias
        if self.bias_collection is None:
            print('No bias frames found. Returning None')
            return None
        master_bias = ccdp.combine(
            self.bias_collection.ccds(ccd_kwargs={'unit': 'adu'}), 
            output_file=None,
            method='average', sigma_clip=True, 
            sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, 
            sigma_clip_func=np.ma.median)
        master_bias.header['HISTORY'] = 'Master bias frame created'
        master_bias.mask = None
        master_bias.data = master_bias.data.astype(np.float32)
        master_bias.uncertainty.array = master_bias.uncertainty.array.astype(np.float32)
        master_bias = SpecImage(master_bias)
        if save==True:
            output_dir = Path(self.work_dir) / 'MasterFrames'
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = Path(output_dir) / output_file
            self.master_file = output_file
            master_bias.write(output_file, overwrite=True)
        self.master_bias = master_bias
        return self.master_bias
