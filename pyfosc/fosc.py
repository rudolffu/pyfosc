#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
from astropy.io import fits
from astropy.io import registry
from astropy.wcs import WCS
from astropy.nddata import CCDData, NDDataArray
from astropy.nddata.mixins.ndio import NDDataRead
from astropy.table import vstack,Table
from astropy.time import Time
from astropy.stats import sigma_clip
import astropy.units as u
import pandas as pd
import numpy as np
import os
import re
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
import yaml




class BasicCCDMixin():
    """A mixin class for CCDData objects"""
    
    read = registry.UnifiedReadWriteMethod(NDDataRead)
    
    def plot_image(self, log_scale=False, use_wcs=False, cmap='viridis',
                   vmin=None, vmax=None, percentile_range=(1, 99), 
                   normalization=None, ax=None):
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
        ax : matplotlib.axes.Axes, optional
            The existing axes to plot in. If None, plot in a new figure.
        
        Returns
        -------
        None
        """
        data = self.data
        framename = self.framename
        frametype = self.frametype
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
            
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 12),
                                   subplot_kw={'projection': projection})
        else:
            fig = ax.figure
        im = ax.imshow(data, cmap=cmap, norm=norm, origin='lower')
        if framename is not None or frametype is not None:
            ax.set_title(f'{frametype} {framename}')
        plt.colorbar(im, ax=ax, label=self.unit, 
                     orientation='horizontal', shrink=0.6, pad=0.1)
        return ax
        
        # im = ax.imshow(data, cmap=cmap, norm=norm, origin='lower')
        # plt.colorbar(im, ax=ax, label=self.unit, 
        #              orientation='horizontal', shrink=0.6, pad=0.1)
        # plt.show()
        # return ax
    

class SpecImage(BasicCCDMixin, CCDData):
    """A class for 2d spectral images"""
    def __init__(self, *args, disp_axis=None, framename=None, frametype=None, **kwd):
        """
        Enhanced constructor to handle additional properties specific to spectral images.
        """
        if len(args) == 1 and (isinstance(args[0], str) or isinstance(args[0], Path)):
            unit = kwd.get('unit', u.dimensionless_unscaled)
            ccddata = super().read(args[0], unit=unit, **kwd)
            super().__init__(ccddata, **kwd)
            framename = os.path.basename(args[0])
        else:
            super().__init__(*args, **kwd)
        self._disp_axis = disp_axis
        self.framename = framename or (os.path.basename(kwd.get('filename')) if 'filename' in kwd else None)
        self.frametype = frametype or self.meta.get('OBSTYPE', None)
        
    # @classmethod
    # def read(cls, filename, fix=False, **kwd):
    #     """
    #     Reads a spectral image from a file and returns a SpecImage object.
    #     """
    #     if fix is True:
    #         hdu = fits.open(filename, mode='update')
    #         hdu[0].verify('fix')
    #         hdu.flush()
    #         hdu.close()
    #     ccddata = super().read(filename, **kwd)
    #     framename = os.path.basename(filename)
    #     frametype = ccddata.meta.get('OBSTYPE', None)
    #     return cls(ccddata, framename=framename, frametype=frametype, **kwd)
        
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


class FOSCFileCollection(ImageFileCollection):
    """A class for managing FOSC data files based on ImageFileCollection"""
    def __init__(self, location=None, keywords=None, 
                 find_fits_by_reading=False, filenames=None, 
                 glob_include=None, glob_exclude=None, ext=0):
        """
        Parameters
        ----------
        location : str or Path
            The location of the data files.
        keywords : dict, optional
            The keywords to use for filtering the data files.
        find_fits_by_reading : bool, optional
            If True, search for FITS files by reading the headers.
        filenames : list, optional
            The list of filenames to use.
        glob_include : str, optional
            The glob pattern to include files.
        glob_exclude : str, optional    
            The glob pattern to exclude files.
        ext : int, optional
            The extension of the FITS file to use.
        """
        super().__init__(location, keywords, 
                         find_fits_by_reading, 
                         filenames, glob_include, 
                         glob_exclude, ext)
        self.table = self.summary
        # self.parameters = {}
        
    def check_groups(self, grism=None, slit=None, telescope=None):
        """
        Check the groups of files based on grism, slit, and telescope.
        
        Parameters
        ----------
        grism : str, optional
            The grism name.
        slit : str, optional
            The slit name.
        telescope : str, optional
            The telescope name.
        """
        tbs = self.table
        if telescope is None:
            if 'telescop' in tbs.columns:
                telescope_str = tbs['telescop'][0]
                if 'Xinglong' in telescope_str:
                    telescope = 'XLT'
                elif '2m4' in telescope_str:
                    telescope = 'LJT'
                elif 'HCT' in telescope_str:
                    telescope = 'HCT'
            elif 'telid' in tbs.columns:
                telescope_str = tbs['telid'][0]
                if telescope_str == '200':
                    telescope = 'P200'
                    self.grism = tbs['fpa'][0]
        self.telescope = telescope
        if self.telescope != 'P200':
            filter_arr = np.unique(tbs['filter'].value.data)
        if self.telescope == 'XLT':
            filter_components = [re.split(r'_|\*', s) for s in filter_arr]
            filter_components = np.unique(filter_components) 
            self.slit = self.find_filter_comp(
                filter_components,'slit', comp=slit)
            self.grism = self.find_filter_comp(
                filter_components,'grism',comp=grism, comp_key='G')
            self.list_bias = self.files_filtered(obstype='BIAS').tolist()
            self.list_flat = self.files_filtered(obstype='SPECLFLAT').tolist()
            self.list_cal = self.files_filtered(obstype='SPECLLAMP').tolist()
            self.list_sci = self.files_filtered(
                regex_match=True,
                obstype='SPECLTARGET|SPECLFLUXREF').tolist()
            self.list_allbutbias = self.files_filtered(
                regex_match=True, 
                obstype='SPECLFLAT|SPECLLAMP|SPECLTARGET|SPECLFLUXREF').tolist()
            self.list_sci_cal = self.files_filtered(
                regex_match=True, 
                obstype='SPECLTARGET|SPECLLAMP|SPECLFLUXREF').tolist()
            self.list_slitimg = self.files_filtered(obstype='SLITTARGET').tolist()
        elif self.telescope == 'LJT':
            filter_components = [re.split(r'\*', s) for s in filter_arr]
            filter_components = np.unique(filter_components) 
            self.slit = self.find_filter_comp(
                filter_components,'slit', comp=slit)
            self.grism = self.find_filter_comp(
                filter_components,'grism',comp=grism)
            self.list_sci = list(self.files_filtered(file='sci*', regex_match=True))
            naxis1 = fits.getval(
                os.path.join(self.location, self.list_sci[0]), 'naxis1')
            naxis2 = fits.getval(
                os.path.join(self.location, self.list_sci[0]), 'naxis2')
            self.list_bias = list(self.files_filtered(file='bias*', naxis1=naxis1, naxis2=naxis2, regex_match=True))
            self.list_flat = list(self.files_filtered(file='lampflat*', naxis1=naxis1, naxis2=naxis2, regex_match=True))
            self.list_cal = list(self.files_filtered(file='cal*', regex_match=True))
            self.list_allbutbias = self.list_flat + self.list_cal + self.list_sci
            self.list_sci_cal = self.list_cal + self.list_sci
        elif self.telescope == 'HCT':
            filter_arr = np.unique(tbs['grism'].value.data)
            self.grism = self.find_filter_comp(
                filter_arr, 'Grism', comp=grism)
            self.slit = None
            self.list_bias = list(self.files_filtered(file='bias*', regex_match=True))
            if len(self.list_bias) > 0:
                naxis1 = fits.getval(
                    os.path.join(self.location, self.list_bias[0]), 'naxis1')
                naxis2 = fits.getval(
                    os.path.join(self.location, self.list_bias[0]), 'naxis2')
                self.list_sci = list(self.files_filtered(file='object*', 
                                                         naxis1=naxis1, 
                                                         naxis2=naxis2, 
                                                         regex_match=True))
            else:
                self.list_sci = list(self.files_filtered(file='object*', regex_match=True))
            self.list_flat = list(self.files_filtered(file='lampflat*', regex_match=True))
            self.list_cal = list(self.files_filtered(file='cal*', regex_match=True))
            self.list_allbutbias = self.list_flat + self.list_cal + self.list_sci
            self.list_sci_cal = self.list_cal + self.list_sci
        elif self.telescope == 'P200':
            self.slit = None
            self.list_bias = list(self.files_filtered(file='(b|r)bias*', regex_match=True))
            self.list_sci = list(self.files_filtered(file='(b|r)object*', regex_match=True))
            self.list_flat = list(self.files_filtered(file='(b|r)flat*', regex_match=True))
            self.list_cal = list(self.files_filtered(file='(b|r)cal*', regex_match=True))
            self.list_allbutbias = self.list_flat + self.list_cal + self.list_sci
            self.list_sci_cal = self.list_cal + self.list_sci
    
    def find_filter_comp(self, filter_components, comp_name, comp=None, comp_key=None):
        """
        Find a filter component (slit or grism) based on the filter strings.
        
        Parameters
        ----------
        filter_components : list
            The list of filter components.
        comp_name : str
            The name of the filter component.
        comp : str, optional
            The filter component to search for.
        comp_key : str, optional
            The matching pattern to use for searching the filter component.
        """
        if comp_key is None:
            comp_key = comp_name
        comp_list = [s for s in filter_components 
                     if ((comp_key in s) or (comp_name in s))]
        if len(comp_list) == 1:
            comp_guessed = comp_list[0]
            print(f'{comp_name} found:', comp_guessed)
            if comp is not None:
                if comp_guessed == comp:
                    print(f'{comp_name} matches the input ({comp_name}={comp})')
                else:
                    comp_list = [comp_guessed, comp]
                    print(f'{comp_name} does not match the input '
                          f'(1. {comp_name} found={comp_guessed}, '
                          f'2. {comp_name} input={comp})')
                    comp_index = int(
                        input(f'Enter the index of the {comp_name} '
                              f'(enter 1 for {comp_guessed} and 2 for {comp}.): '))
                    comp = comp_list[comp_index-1]
            else:
                comp = comp_guessed       
        elif len(comp_list) > 1:
            print(f'Multiple {comp_name}s found:')
            for i, comp in enumerate(comp_list):
                print(f'Index = {i+1}, name = {comp}')
            comp_index = int(input(f'Enter the index of the {comp_name}: '))
            comp = comp_list[comp_index-1]
        else:
            print(f'No {comp_name} found. Please check the {comp_name} names.')
        return comp      

    def set_parameters(self, grism=None, telescope=None, config_file=None):
        """
        Set parameters based on grism and telescope, optionally overridden by a YAML configuration file.
        
        Parameters
        ----------
        grism : str
            The grism name.
        telescope : str
            The telescope name.
        config_file : str, optional
            Path to a YAML configuration file containing user-defined parameters.
        """
        # Default parameters
        if telescope is None:
            telescope = self.telescope
        if grism is None:
            grism = self.grism
        grisms_ljt_naming = {
            'grism3': 'G3', 
            'grism8': 'G8', 
            'grism10': 'G10', 
            'grism14': 'G14'} 
        grisms_hct_naming = {
            '4 Grism 7': 'Gr7'
        }
        if telescope == 'LJT' and grism in grisms_ljt_naming.keys():
            grism = grisms_ljt_naming[grism]
        if telescope == 'HCT' and grism in grisms_hct_naming.keys():
            grism = grisms_hct_naming[grism]
        default_params = {
            'XLT_NewG4': {
                'trimsec': '[51:1750,681:1350]', 
                'disp_axis': 1, 
                'overscan': False
                },
            'XLT_G4': {
                'trimsec': '[51:1750,681:1350]', 
                'disp_axis': 1, 
                'overscan': False
                },
            'LJT_G3': {
                # 'trimsec': '[651:1350,2301:4130]',
                'trimsec': '[201:1400,2301:4130]',
                'biassec': '[10:40,1:4612]',
                'disp_axis': 0,
                'overscan': True
                },
            'LJT_G14': {
                'trimsec': '[751:1450,2280:4200]',
                'biassec': '[10:40,1:4612]',
                'disp_axis': 0,
                'overscan': True
                },
            'LJT_G8': {
                'trimsec': '[751:1450,1136:4150]',
                'biassec': '[10:40,1:4612]',
                'disp_axis': 0,
                'overscan': True
                },
            'HCT_Gr7': {
                'trimsec': '[26:250,165:2800]',
                'disp_axis': 0,
                'overscan': False
                },
            'P200_DBSP_BLUE': {
                'trimsec': '[101:380,240:2585]',
                'biassec': '[420:460,1:2835]',
                'disp_axis': 0,
                'overscan': True
                },
            'P200_DBSP_RED2': {
                'trimsec': '[751:3700,51:400]',
                'disp_axis': 1,
                'overscan': False
            },
        }
        
        # Fetch the default parameters for the given telescope-grism pair
        key = f"{telescope}_{grism}"
        params = default_params.get(key, {})
        
        # Load user-defined parameters from YAML if provided
        if config_file:
            with open(config_file, 'r') as f:
                user_params = yaml.safe_load(f)
                
                # Check if the user_params contains the specific telescope-grism key
                if key in user_params:
                    # Update the default parameters with any user overrides
                    params.update(user_params[key])
        
        self.parameters = params
        print(f"Parameters set for {telescope}-{grism}: {self.parameters}")

    def write_parameters_to_file(self, output_file, grism=None, telescope=None):
        """
        Write the current parameters to a YAML file, including grism and telescope keys.
        
        Parameters
        ----------
        output_file : str
            The path to the output YAML file where parameters will be written.
        grism : str
            The grism name.
        telescope : str
            The telescope name.
        """
        if not self.parameters:
            print("No parameters to write.")
            return
        if telescope is None:
            telescope = self.telescope
        if grism is None:
            grism = self.grism
        # Creating a nested dictionary with grism and telescope as keys
        params_to_write = {f"{telescope}_{grism}": self.parameters}
        
        with open(output_file, 'w') as f:
            yaml.dump(params_to_write, f, default_flow_style=False)
            
        print(f"Parameters written to {output_file}.")
