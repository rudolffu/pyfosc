#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.ndimage import convolve, gaussian_filter
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import CCDData
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
import astropy.units as u
import ccdproc as ccdp
from ccdproc import ImageFileCollection
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
        if save is True:
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


class FlatNormalizer:
    """Class to normalize flat field frames."""
    def __init__(self, master_flat, disp_axis,
                 work_dir=None):
        """
        Parameters
        ----------
        master_flat : SpecImage(ccdproc.CCDData)
            The master flat field frame.
        disp_axis : int
            The dispersion axis of the spectrograph. 
            0 for vertical, 1 for horizontal (python indexing).
        work_dir : str, optional
            The working directory. The default is None.
        """
        self.master_flat = master_flat
        self.disp_axis = disp_axis
        self.normalized_flat = None
        if work_dir is None:
            work_dir = os.getcwd()
        self.work_dir = work_dir
        
    def response_correct(self, degree=25):
        """
        Correct for the (near blackbody) spectrum response of the halogen lamp.
        
        Parameters
        ----------
        degree : int, optional
            The degree of the Legendre1D polynomial fit. The default is 25. 
        """
        # check if disp_axis is 1 (horizontal). if 0, transpose the data
        if hasattr(self, 'resp_corr_flat'):
            if self.resp_corr_flat is not None:
                print('Response corrected flat already exists. Returning the flat.')
                self.resp_corr_flat.plot_image()
                return self.resp_corr_flat
        master_flat = self.master_flat.copy()
        if self.disp_axis == 0:
            master_flat.data = master_flat.data.T
            master_flat.uncertainty.array = master_flat.uncertainty.array.T
        length_spatial = master_flat.data.shape[0]
        length_disp = master_flat.data.shape[1]
        # get the central 1/3 rows of the flat
        spatial_start = length_spatial // 3
        spatial_end = 2 * length_spatial // 3
        cent_section = master_flat.data[spatial_start:spatial_end, :]
        cent_median = np.median(cent_section, axis=0)
        fitter = fitting.FittingWithOutlierRemoval(
            fitting.LinearLSQFitter(), sigma_clip, 
            niter=5,
            sigma=2)
        model_init = models.Legendre1D(degree=degree)
        model_fit, mask = fitter(model_init, np.arange(length_disp), cent_median)
        # plot the fit and the residuals
        x = np.arange(length_disp)
        y = model_fit(x)
        residuals = cent_median - y
        fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True,
                                 gridspec_kw={'height_ratios': [2, 1, 6],
                                              'hspace': 0})
        axes[0].plot(cent_median, color='black', label='data')
        axes[0].set_title(
            f'Median profile of row {spatial_start} '
            f'to {spatial_end} along the dispersion axis')
        axes[0].plot(y, color='red', label='fit')
        axes[0].legend()
        axes[1].plot(residuals, color='black', label='residuals')
        axes[1].legend()
        master_flat.plot_image(ax=axes[2])
        axes[2].set_aspect('auto')
        axes[2].set_title(None)
        plt.show()
        # divide the 2d master_flat by the 1d fit
        resp_corr_flat = master_flat.copy()
        resp_corr_flat.data /= y
        resp_corr_flat.uncertainty.array /= y
        resp_corr_flat.unit = u.dimensionless_unscaled
        resp_corr_flat.uncertainty.unit = u.dimensionless_unscaled
        if self.disp_axis == 0:
            resp_corr_flat.data = resp_corr_flat.data.T
            resp_corr_flat.uncertainty.array = resp_corr_flat.uncertainty.array.T
        self.resp_corr_flat = resp_corr_flat
        resp_corr_flat.plot_image()
        
    def illumination_correct(self, gauss_sigma=10, save=False, output_file=None):
        """
        Correct for spatial variations caused by illumination
        in the flat field with a gaussian filter.
        
        Parameters
        ----------
        gauss_sigma : int, optional
            The sigma value for the gaussian filter. The default is 10.
        save : bool, optional
            Save the normalized flat. The default is False.
        output_file : str, optional
            The name of the output file. The default is 'normalized_flat.fits'.
        """
        # check if self.resp_corr_flat exists, if not, run dispersion_correct
        if not hasattr(self, 'resp_corr_flat'):
            self.dispersion_correct()
        if output_file is None:
            output_file = 'normalized_flat.fits'
        resp_corr_flat = self.resp_corr_flat
        smooth_dc_flat = gaussian_filter(
            resp_corr_flat.data, sigma=gauss_sigma)
        normalized_flat = resp_corr_flat.copy()
        normalized_flat.data /= smooth_dc_flat
        normalized_flat.uncertainty.array /= smooth_dc_flat
        median_of_flat = np.median(normalized_flat.data)
        normalized_flat.data /= median_of_flat
        normalized_flat.uncertainty.array /= median_of_flat
        normalized_flat.unit = u.dimensionless_unscaled
        normalized_flat.framename = 'normalized_flat.fits'
        normalized_flat.header['HISTORY'] = 'Flat normalized'
        normalized_flat.header['BUNIT'] = ''
        normalized_flat.plot_image()
        self.normalized_flat = normalized_flat
        if save is True:
            output_dir = Path(self.work_dir) / 'MasterFrames'
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = Path(output_dir) / output_file
            normalized_flat.write(output_file, overwrite=True)
            print(f'Normalized flat saved as {output_file}')