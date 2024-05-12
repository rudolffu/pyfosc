#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import vstack,Table
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, EarthLocation, AltAz
import glob
import os
import re
from pathlib import Path
import shutil
import multiprocessing as mp
import json
from astropy.time import Time
from ccdproc import ImageFileCollection
import ccdproc as ccdp
from astropy.nddata import CCDData, nduncertainty
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.stats import mad_std
# import specreduce
from specreduce.tracing import FitTrace, ArrayTrace
from specreduce.background import Background
from specreduce.extract import HorneExtract, BoxcarExtract
# from specutils import Spectrum1D, SpectrumCollection
from astropy.stats import sigma_clip
# from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler
from PyAstronomy import pyasl
from astropy.utils.exceptions import AstropyWarning
import warnings
warnings.simplefilter('ignore', category=AstropyWarning)
warnings.simplefilter('ignore', category=UserWarning)
import logging
from .fosc import FOSCFileCollection,SpecImage
from . import newlog as log
logger = log.getLogger('pyfosc.extraction', level=logging.INFO)
logger.propagate = False

class Extract1dSpec:
    def __init__(self, ic_2dspec=None, filenames=None, data_dir=None, disp_axis=None,
                 trace_list=None, qa_dir=None):
        if ic_2dspec is not None:
            if not isinstance(ic_2dspec, (ImageFileCollection, FOSCFileCollection)):
                raise ValueError(
                    'ic_2dspec must be an instance of ImageFileCollection or FOSCFileCollection')
        elif filenames is not None:
            if isinstance(filenames, str):
                filenames = [filenames]
            if data_dir is not None:
                ic_2dspec = FOSCFileCollection(location=data_dir, filenames=filenames)
            else:
                ic_2dspec = FOSCFileCollection(location='./', filenames=filenames)
        self.ic_2dspec = ic_2dspec
        ic_2dspec.check_groups()
        ic_2dspec.set_parameters()
        self.telescope = ic_2dspec.telescope
        self.filenames = ic_2dspec.files
        self.data_dir = ic_2dspec.location
        if disp_axis is None:
            disp_axis = ic_2dspec.parameters['disp_axis']
        self.disp_axis = disp_axis
        self.spatial_axis = 1 - disp_axis
        if qa_dir is None:
            qa_dir = Path(self.data_dir)/'../QAplots'
            qa_dir.mkdir(exist_ok=True)
        self.trace_list = trace_list
        self.qa_dir = qa_dir
        peak_list = [np.argmax(np.sum(im.data, 
                                      axis=self.disp_axis)) for im in self.ic_2dspec.ccds()]
        median_peak = np.median(peak_list)
        self.median_peak = median_peak
        
    def trace_and_extract(self, bins=None, peak_method='max', window=50, guess=None,
                          sep_one_side=2, width_one_side=20):        
        trace_list = []    
        for i, (im, fname) in enumerate(self.ic_2dspec.ccds(return_fname=True)):
            im = SpecImage(im)
            if self.disp_axis == 0:
                im.data = im.data.T
                im.uncertainty.array = im.uncertainty.array.T
            sum_pro = np.sum(im.data, axis=1)
            if guess is None:
                guess = self.median_peak
            if abs(sum_pro.argmax()-guess)<=50:
                guess = sum_pro.argmax()   
            trace = FitTrace(
                im, 
                bins=bins, 
                peak_method=peak_method, 
                window=window,
                guess=guess)
            g_cent = np.mean(trace.trace.data)
            trace_upper = int(np.max(trace.trace.data)) + 5
            trace_lower = int(np.min(trace.trace.data)) - 5
            if isinstance(im.uncertainty, nduncertainty.StdDevUncertainty):
                trace_uncertainty = im.uncertainty.array[trace_lower:trace_upper, :]
            elif isinstance(im.uncertainty, nduncertainty.VarianceUncertainty):
                im_stddev = np.sqrt(im.uncertainty.array)
                trace_uncertainty = im_stddev[trace_lower:trace_upper, :]
            elif isinstance(im.uncertainty, nduncertainty.InverseVariance):
                im_stddev = np.sqrt(1/im.uncertainty.array)
                trace_uncertainty = im_stddev[trace_lower:trace_upper, :]
            trace_uncertainty_1d = np.sqrt(np.sum(trace_uncertainty**2, axis=0))
            trace_list.append(trace.trace.data)
            if self.telescope == "LJT":
                cent_pro = sum_pro[int(g_cent)-50:int(g_cent)+50] 
                # fit a gaussian to the summed profile at the guessed location
                g_init = models.Gaussian1D(
                    amplitude=np.max(cent_pro), mean=50, stddev=3., bounds={'stddev': (1.5, 10)})
                cbg_init = models.Const1D(amplitude=np.min(cent_pro))
                g_fitter = fitting.LMLSQFitter()
                gc_init = g_init + cbg_init
                gc_fitted = g_fitter(gc_init, np.arange(len(cent_pro)), cent_pro)
                std_pro = gc_fitted.stddev_0.value
                # g_fitted = g_fitter(g_init, np.arange(len(cent_pro)), cent_pro)
                # std_pro = g_fitted.stddev.value
                logger.info(
                    f'Fitting Gaussian to the central 100 pixels of the summed profile of {fname}: \n\t\t'
                    f'center={g_cent:.2f}, std={std_pro:.2f}')
                bg = Background.two_sided(im, trace, separation=10*std_pro, width=4*std_pro)
                extract = HorneExtract(im-bg, trace)
                sp = extract.spectrum
                if np.min(sp.flux.value) < np.median(sp.flux.value) - 3 * mad_std(sp.flux.value):
                    logger.warning(
                        f'{fname} has a minimum flux value: {np.min(sp.flux.value)} < median - 3 * mad_std.\n\t\t'
                        f'Redoing background subtraction with larger separation and width.')
                    bg = Background.two_sided(im, trace, separation=12*std_pro, width=4*std_pro)
                    extract = HorneExtract(im-bg, trace)
                    sp = extract.spectrum
            elif self.telescope == "XLT":
                bg = Background.one_sided(im, trace, separation=sep_one_side, width=width_one_side)
                extract = HorneExtract(im-bg, trace)
                sp = extract.spectrum
            multispecdata = fake_multispec_data(
                (sp.flux.value, sp.flux.value, 
                 bg.bkg_spectrum().flux.value, trace_uncertainty_1d))
            sp_hdu = fits.PrimaryHDU(data=multispecdata)
            hdrcopy = im.header.copy(strip = True)
            sp_hdu.header.extend(hdrcopy, strip=True, update=True,
                                 update_first=False, useblanks=True, bottom=False)
            sp_hdr = sp_hdu.header
            sp_hdr['NAXIS'] = 3
            sp_hdr['NAXIS1'] = len(sp.flux)
            sp_hdr['NAXIS2'] = 1
            sp_hdr['NAXIS3'] = 4
            sp_hdr['WCSDIM'] = 3
            sp_hdr['WAT0_001'] = 'system=equispec'
            sp_hdr['WAT1_001'] = 'wtype=linear label=Pixel'
            sp_hdr['WAT2_001'] = 'wtype=linear'
            sp_hdr['CRVAL1'] = 1
            sp_hdr['CRPIX1'] = 1
            sp_hdr['CD1_1'] = 1
            sp_hdr['CD2_2'] = 1
            sp_hdr['CD3_3'] = 1
            sp_hdr['LTM1_1'] = 1
            sp_hdr['LTM2_2'] = 1
            sp_hdr['LTM3_3'] = 1
            sp_hdr['WAT3_001'] = 'wtype=linear'
            sp_hdr['CTYPE1'] = 'PIXEL'
            sp_hdr['CTYPE2'] = 'LINEAR'
            sp_hdr['CTYPE3'] = 'LINEAR'
            sp_hdr['BANDID1'] = 'spectrum - background fit, weights variance, clean no'               
            sp_hdr['BANDID2'] = 'raw - background fit, weights none, clean no'                        
            sp_hdr['BANDID3'] = 'background - background fit'                                         
            sp_hdr['BANDID4'] = 'sigma - background fit, weights variance, clean no'  
            sp_hdr['APNUM1'] = f'1 1 {trace.trace.data[0]:.2f} {trace.trace.data[-1]:.2f}'
            sp_hdu.writeto(f'{self.data_dir}/a{fname}', overwrite=True)
            fig, axes = plt.subplots(nrows=2, figsize=(8, 8))
            im.plot_image(ax=axes[0])
            # axes[0].imshow((im-bg).data, origin='lower', aspect='auto')
            axes[0].plot(trace.trace.data, 'r-')
            axes[0].set_title(f'trace of {fname}')
            ax1 = axes[1].get_subplotspec()
            axes[1].remove()
            axes[1] = fig.add_subplot(ax1, projection=sp.wcs)
            sp.plot(axes=axes[1])
            axes[1].set_title(f'Extracted 1d spectrum of {fname}')
            plt.savefig(f'{self.qa_dir}/trace_{fname}.pdf')
        self.trace_list = trace_list
        
    def extract_lamp(self, ic_lamp=None, filenames=None,
                     mean_trace=None):
        if ic_lamp is None:
            if filenames is not None:
                if isinstance(filenames, str):
                    filenames = [filenames]
                ic_lamp = FOSCFileCollection(location=self.data_dir, filenames=filenames)
            else:
                raise ValueError('ic_lamp or list of lamp files must be provided')
        if mean_trace is None:
            mean_trace = np.mean(self.trace_list, axis=0)
        for i, (im, fname) in enumerate(ic_lamp.ccds(return_fname=True)):
            im = SpecImage(im)
            if self.disp_axis == 0:
                im.data = im.data.T
                im.uncertainty.array = im.uncertainty.array.T
            trace = ArrayTrace(im, mean_trace)
            extract = BoxcarExtract(im, trace)
            sp = extract.spectrum
            sp_hdr = im.header
            sp_hdr['CTYPE1'] = 'PIXEL'
            sp_hdr['CRVAL1'] = 1
            sp_hdr['CRPIX1'] = 1
            sp_hdr['CD1_1'] = 1
            pyasl.write1dFitsSpec(
                f'{self.data_dir}/a{fname}', 
                flux=sp.flux.value, 
                wvl=sp.spectral_axis.value, 
                header=sp_hdr,
                clobber=True)
            fig, ax = plt.subplots(figsize=(8, 4),
                                subplot_kw={'projection': sp.wcs})
            sp.plot(axes=ax)
            ax.set_title(f'Extracted 1d spectrum of {fname}')
            plt.savefig(f'{self.qa_dir}/lamp1d_{fname}.pdf')
            
            
def fake_multispec_data(arrlist):
# https://github.com/jrthorstensen/opextract/blob/master/opextract.py#L337
   # takes a list of 1-d numpy arrays, which are
   # to be the 'bands' of a multispec, and stacks them
   # into the format expected for a multispec.  As of now
   # there can only be a single 'aperture'.

   return np.expand_dims(np.array(arrlist), 1)