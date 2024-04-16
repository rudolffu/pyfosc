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
from astropy.nddata import CCDData
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.stats import mad_std
import specreduce
from specreduce.tracing import FitTrace, ArrayTrace
from specreduce.background import Background
from specreduce.extract import HorneExtract, BoxcarExtract
from specutils import Spectrum1D, SpectrumCollection
from astropy.stats import sigma_clip
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler
from PyAstronomy import pyasl
from astropy.utils.exceptions import AstropyWarning
from astropy.io.fits.verify import VerifyWarning
from astropy.wcs import FITSFixedWarning
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
        
    def trace_and_extract(self, bins=200, peak_method='max', window=50, guess=None):
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
            trace_list.append(trace.trace.data)
            bg_init = Background.two_sided(im, trace, separation=20, width=5)
            im_bgs = im - bg_init
            new_sum_pro = np.sum(im_bgs.data[:, int(im.data.shape[1]/3):2*int(im.data.shape[1]/3)], axis=1)
            cent_pro = new_sum_pro[int(guess)-50:int(guess)+50] 
            # fit a gaussian to the summed profile at the guessed location
            g_init = models.Gaussian1D(amplitude=np.max(cent_pro), 
                                       mean=50, stddev=3., bounds={'stddev': (1.5, 15)})
            g_fitter = fitting.LMLSQFitter()
            g_fitted = g_fitter(g_init, np.arange(len(cent_pro)), cent_pro)
            std_pro = g_fitted.stddev.value
            logger.info(
                f'Fitting Gaussian to the central 100 pixels of the summed profile of {fname}: \n\t\t'
                f'guess={guess:.2f}, std={std_pro:.2f}')
            bg = Background.two_sided(im, trace, separation=1.4*std_pro, width=0.5*std_pro)
            extract = HorneExtract(im-bg, trace)
            sp = extract.spectrum
            if np.min(sp.flux.value) < -1 * np.median(sp.flux.value) - 3*mad_std(sp.flux.value) and self.telescope == 'LJT':
                logger.warning(
                    f'{fname} has a minimum flux value: {np.min(sp.flux.value)} < median - 3 * mad_std.\n\t\t'
                    f'Redoing background subtraction with larger separation and width.')
                bg = Background.two_sided(im, trace, separation=30, width=12)
                extract = HorneExtract(im-bg, trace)
                sp = extract.spectrum
            sp_hdr = im.header
            sp_hdr['CTYPE1'] = 'PIXEL'
            sp_hdr['CTYPE2'] = 'LINEAR'
            sp_hdr['CTYPE3'] = 'LINEAR'
            sp_hdr['CRVAL1'] = 1
            sp_hdr['CRPIX1'] = 1
            sp_hdr['CD1_1'] = 1
            sp_hdr['WAT3_001'] = 'wtype=linear'
            pyasl.write1dFitsSpec(
                f'{self.data_dir}/a{fname}', 
                flux=sp.flux.value, 
                wvl=sp.spectral_axis.value, 
                header=sp_hdr,
                clobber=True)
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