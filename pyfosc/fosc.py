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
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle



class BasicCCDMixin():
    """A mixin class for CCDData objects"""
    def plot_image(self, log_scale=False, use_wcs=False, cmap='viridis',
                   vmin=None, vmax=None, percentile_range=(1, 99), 
                   normalization=None):
        """
        Plot the image data with options for scaling, WCS projection, and colormap.
        
        Parameters
        ----------
        log_scale : bool, optional
            If True, use log scale for the image. Ignored if normalization is provided.
        use_wcs : bool, optional
            If True, use the WCS projection for the image.
        cmap : str, optional
            The colormap for the image.
        vmin, vmax : float, optional
            The minimum and maximum value for the image scaling. 
            If None, determined from percentile_range.
        percentile_range : tuple, optional
            The lower and upper percentiles to use for automatic vmin/vmax calculation.
        normalization : matplotlib.colors.Normalize or subclass, optional
            Custom normalization. Overrides log_scale if provided.
        
        Returns
        -------
        None
        """
        data = self.data
        if vmin is None or vmax is None:
            vmin, vmax = np.percentile(data, percentile_range)
        
        if normalization:
            norm = normalization(vmin=vmin, vmax=vmax)
        elif log_scale==True:
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)
        
        if use_wcs and self.wcs is not None:
            projection = self.wcs
        else:
            projection = None
            
        fig, ax = plt.subplots(subplot_kw={'projection': projection})
        im = ax.imshow(data, cmap=cmap, norm=norm, origin='lower')
        plt.colorbar(im, ax=ax, label=self.unit)
        plt.show()
    
    

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
        self._disp_axis = int(value)


    def plot_sampled_rows_columns(self):
        """
        Plot each sampled row and column at 1/4, 1/2, and 3/4 of the image size.
        """
        data = self.data
        nrows, ncols = data.shape
        
        # Calculate indices for rows and columns at 1/4, 1/2, and 3/4 positions
        rows_indices = [nrows // 4, nrows // 2, 3 * nrows // 4]
        cols_indices = [ncols // 4, ncols // 2, 3 * ncols // 4]
        
        # Determine the layout of subplots
        nrows_plot = 2  # Number of rows of subplots
        ncols_plot = 3  # Number of columns of subplots
        
        # Create figure with determined layout
        fig, axes = plt.subplots(nrows_plot, ncols_plot, figsize=(12, 8))
        
        # Plot sampled rows
        for i, row_index in enumerate(rows_indices):
            axes[0, i].plot(data[row_index, :], label=f'Row {row_index}')
            axes[0, i].set_title(f'Sampled Row at {row_index}')
            axes[0, i].set_xlabel('Column Index')
            axes[0, i].set_ylabel(f'Intensity ({self.unit})')
            axes[0, i].legend()
        
        # Plot sampled columns
        for j, col_index in enumerate(cols_indices):
            axes[1, j].plot(data[:, col_index], label=f'Column {col_index}')
            axes[1, j].set_title(f'Sampled Column at {col_index}')
            axes[1, j].set_xlabel('Row Index')
            axes[1, j].set_ylabel(f'Intensity ({self.unit})')
            axes[1, j].legend()

        # Adjust layout
        plt.tight_layout()
        plt.show()

class SlitImage(BasicCCDMixin, CCDData):
    def plot_with_zoom(self, zoom_factor=2, percentile_range=(0.05, 99.95), cmap='viridis'):
        """
        Plots full image and zoomed view of central part, with box indicating zoomed
        area on full image, and labels zoomed region using full image coordinates.
        
        Parameters
        ----------
        zoom_factor : int or float, optional
            Factor by which to zoom into central part of image.
        percentile_range : tuple, optional
            Percentile range for scaling zoomed view to reveal bright sources.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
        
        vmin, vmax = np.percentile(self.data, (1, 99))  # Full image scaling
        ax1.imshow(self.data, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
        ax1.set_title('Full Image')
        
        nrow, ncol = self.data.shape
        center_row, center_col = nrow // 2, ncol // 2  # Zoom center row/column
        zext_row, zext_col = nrow // zoom_factor // 2, ncol // zoom_factor // 2  # Zoom extents
        
        zoom_data = self.data[center_row-zext_row:center_row+zext_row, center_col-zext_col:center_col+zext_col]
        vzmin, vzmax = np.percentile(zoom_data, percentile_range)  # Zoomed view scaling
        
        ax2.imshow(zoom_data, cmap=cmap, origin='lower', norm=LogNorm(vmin=vzmin, vmax=vzmax))
        ax2.set_title('Zoomed Central Part')
        
        # Box on full image for zoomed region
        rect = Rectangle((center_row-zext_col, center_col-zext_row), 
                         2*zext_col, 2*zext_row, lw=1, edgecolor='r', facecolor='none')
        ax1.add_patch(rect)
        
        # Adjust zoomed image's ticks to match full image's coordinates
        xticks = np.linspace(0, 2*zext_col, 10)
        yticks = np.linspace(0, 2*zext_row, 10)
        ax2.set_xticks(xticks)
        ax2.set_yticks(yticks)
        ax2.minorticks_on()
        ax2.set_xticklabels((center_row-zext_col + xticks).astype(int))
        ax2.set_yticklabels((center_col-zext_row + yticks).astype(int))
        ax2.grid(True)
        
        plt.tight_layout()
        plt.show()
