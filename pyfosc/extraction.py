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
import specreduce
from specreduce.tracing import FitTrace, ArrayTrace
from specreduce.background import Background
from specreduce.extract import HorneExtract, BoxcarExtract
from specutils import Spectrum1D, SpectrumCollection
from astropy.stats import sigma_clip
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler
from PyAstronomy import pyasl

from .fosc import FOSCFileCollection,SpecImage


class Extract1dSpec:
    def __init__(self, ic_2dspec=None, filenames=None, data_dir=None, trace_list=None, qa_dir=None):
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
        self.filenames = ic_2dspec.files
        self.data_dir = ic_2dspec.location
        if qa_dir is None:
            qa_dir = Path(self.data_dir)/'../QAplots'
            qa_dir.mkdir(exist_ok=True)
        self.trace_list = trace_list
        self.qa_dir = qa_dir
        peak_list = [np.argmax(np.sum(im.data, axis=1)) for im in self.ic_2dspec.ccds()]
        median_peak = np.median(peak_list)
        self.median_peak = median_peak
        
    def trace_and_extract(self, bins=200, peak_method='max', window=50, guess=None):
        trace_list = []    
        for i, (im, fname) in enumerate(self.ic_2dspec.ccds(return_fname=True)):
            im = SpecImage(im)
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
            bg = Background.one_sided(im, trace, separation=2, width=20)
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