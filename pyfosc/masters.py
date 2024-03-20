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




class MasterBuilder:
    """Base class to create and handle master frames."""
    def __init__(self, 
                 work_dir=None,
                 list_frames=None, 
                 master_file=None,
                 frametype=None):
        """
        Parameters
        ----------
        work_dir : str, optional
            The working directory where the data is located. The default is None.
        list_frames : list, optional
            List of frames. The default is None.
        master_file : str, optional
            The name of the master frame file. The default is None.   
        frametype : str, optional
            The type of frame. The default is None.
            exapmle: 'bias', 'flat'     
        """
        if master_file is not None:
            self.master_frame = SpecImage.read(master_file)
            return 
        else:
            self.master_frame = None
        if work_dir is None:
            location = os.getcwd()
        else:
            location = Path(work_dir) / 'data'
        if list_frames is not None:
            frame_collection = ImageFileCollection(
                location=location,
                filenames=list_frames)
        self.work_dir = work_dir
        self.list_frames = list_frames
        self.frame_collection = frame_collection
        self.master_file = master_file
        if frametype is None:
            print('No frametype specified. Setting to `None`.')
            frametype = 'None'
        self.frametype = frametype
        
        
    def build(self, save=False, output_file=None, method='average'):
        """
        Create the master frame.
        
        Parameters
        ----------
        save : bool, optional
            Save the master frame. The default is False.
        output_file : str, optional
            The name of the output file. The default is 'master_frame.fits'.
            
        Returns
        -------
        master_frame : SpecImage(ccdproc.CCDData)
            The master frame.
        """
        if output_file is None:
            output_file = f'master_{self.frametype}.fits'
        if self.master_frame is not None:
            print('Master frame already exists. Returning the master frame.')
            return self.master_frame
        if self.frame_collection is None:
            print('No frames found. Returning None')
            return None
        master_frame = ccdp.combine(
            self.frame_collection.ccds(ccd_kwargs={'unit': 'adu'}), 
            output_file=None,
            method=method, sigma_clip=True, 
            sigma_clip_low_thresh=5, sigma_clip_high_thresh=5, 
            sigma_clip_func=np.ma.median)
        master_frame.header['HISTORY'] = f'Master {self.frametype} created'
        master_frame.mask = None
        master_frame.data = master_frame.data.astype(np.float32)
        master_frame.uncertainty.array = master_frame.uncertainty.array.astype(np.float32)
        master_frame = SpecImage(master_frame, 
                                 framename=output_file, 
                                 frametype=self.frametype)
        if save==True:
            output_dir = Path(self.work_dir) / 'MasterFrames'
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = Path(output_dir) / output_file
            self.master_file = output_file
            master_frame.write(output_file, overwrite=True)
        self.master_frame = master_frame
        return self.master_frame


class MasterBias(MasterBuilder):
    """Class to create and handle master bias frames."""
    def __init__(self, 
                 work_dir=None,
                 list_bias=None, 
                 master_file=None):
        super().__init__(work_dir, list_bias,
                         master_file, frametype='bias')


class MasterFlat(MasterBuilder):
    """Class to create and handle master flat frames."""
    def __init__(self, 
                 work_dir=None,
                 list_flat=None, 
                 master_file=None):
        super().__init__(work_dir, list_flat,
                         master_file, frametype='flat')
        
    def make_bad_pixel_mask(self):
        """
        Create a bad pixel mask from a single flat field frame.
        """
        raw_flat_path = os.path.join(
            self.frame_collection.location, 
            self.frame_collection.files[0])
        flat = SpecImage.read(raw_flat_path, unit='adu')
        mask = ccdp.ccdmask(flat)
        self.bad_pixel_mask = mask
        fig, axes = plt.subplots(1, 2, figsize=(15, 10))
        flat.plot_image(ax=axes[0])
        axes[0].set_title('Single raw flat')
        axes[1].imshow(mask, origin='lower', cmap='binary')
        axes[1].set_title('Derived mask')
        plt.show()
