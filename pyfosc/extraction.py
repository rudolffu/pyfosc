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
import specreduce
from specreduce.tracing import FitTrace, ArrayTrace
from specreduce.background import Background
from specreduce.extract import BoxcarExtract
# from specutils import Spectrum1D, SpectrumCollection
from astropy.stats import sigma_clip
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler
from specutils import Spectrum1D
from scipy.interpolate import interp1d
from PyAstronomy import pyasl
from astropy.utils.exceptions import AstropyWarning
from specreduce.core import SpecreduceOperation
from specreduce.tracing import Trace, FlatTrace
from scipy.integrate import trapezoid
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
        
    def trace_and_extract(self, bins=None, peak_method='max', window=50, guess=None):        
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
            trace_list.append(trace.trace.data)
            cent_pro = sum_pro[int(g_cent)-50:int(g_cent)+50] 
            # fit a gaussian to the summed profile at the guessed location
            g_init = models.Gaussian1D(
                amplitude=np.max(cent_pro), mean=50, stddev=3., bounds={'stddev': (1.5, 10)})
            cbg_init = models.Const1D(amplitude=np.min(cent_pro))
            g_fitter = fitting.LMLSQFitter()
            gc_init = g_init + cbg_init
            gc_fitted = g_fitter(gc_init, np.arange(len(cent_pro)), cent_pro)
            std_pro = gc_fitted.stddev_0.value
            logger.info(
                f'Fitting Gaussian to the central 100 pixels of the summed profile of {fname}: \n\t\t'
                f'center={g_cent:.2f}, std={std_pro:.2f}')
            bg = Background.two_sided(im, trace, separation=10*std_pro, width=4*std_pro)
            extract = HorneExtract(im-bg, trace)
            sp = extract.spectrum
            sp_uncertainty = extract.sp_uncertainty
            if np.min(sp.flux.value) < np.median(sp.flux.value) - 3 * mad_std(sp.flux.value):
                logger.warning(
                    f'{fname} has a minimum flux value: {np.min(sp.flux.value)} < median - 3 * mad_std.\n\t\t'
                    f'Redoing background subtraction with larger separation and width.')
                bg = Background.two_sided(im, trace, separation=12*std_pro, width=4*std_pro)
                extract = HorneExtract(im-bg, trace)
                sp = extract.spectrum
                sp_uncertainty = extract.sp_uncertainty
            # check the residual of the background subtraction
            im_res = im-bg
            im_res_med = np.median(im_res.data)
            logger.info(f'{fname} median residual of background subtraction: {im_res_med:.2f}')
            if im_res_med < -2.0:
                logger.warning(
                    f'{fname} background over-subtracted! \n\t\t'
                    f'Using BoxcarExtract for {fname}.')
                extract = BoxcarExtract(im-bg, trace)
                sp = extract.spectrum
                sp_uncertainty = sp.flux.value * 0.1
            if np.isnan(sp.flux.value).any():
                # interpolate over NaNs
                logger.warning(f'Interpolating over NaNs in {fname}')
                flux_interp = interp1d(sp.spectral_axis.value[~np.isnan(sp.flux.value)],
                                        sp.flux.value[~np.isnan(sp.flux.value)],
                                        kind='linear', fill_value='extrapolate')
                sp = Spectrum1D(flux=flux_interp(sp.spectral_axis.value) * sp.flux.unit,
                                spectral_axis=sp.spectral_axis)
            multispecdata = fake_multispec_data(
                (sp.flux.value, sp.flux.value, 
                 bg.bkg_spectrum().flux.value, sp_uncertainty))
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
            # imres = axes[0].imshow((im-bg).data, origin='lower', aspect='auto')
            axes[0].plot(trace.trace.data, 'r-')
            axes[0].set_title(f'trace of {fname}')
            ax1 = axes[1].get_subplotspec()
            axes[1].remove()
            axes[1] = fig.add_subplot(ax1, projection=sp.wcs)
            sp.plot(axes=axes[1])
            axes[1].set_title(f'Extracted 1d spectrum of {fname}')
            # # add a colorbar to the right of the image
            # cbar = fig.colorbar(imres, ax=axes[0], orientation='vertical')
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


class HorneExtract(specreduce.extract.HorneExtract):
    def __call__(self, image=None, trace_object=None,
                 disp_axis=None, crossdisp_axis=None,
                 bkgrd_prof=None, spatial_profile='interpolated_profile',
                 n_bins_interpolated_profile=None,
                 interp_degree_interpolated_profile=None,
                 variance=None, mask=None, unit=None):
        """
        Run the Horne calculation on a region of an image and extract a
        1D spectrum.

        Parameters
        ----------

        image : `~astropy.nddata.NDData`-like or array-like, required
            The input 2D spectrum from which to extract a source. An
            NDData object must specify uncertainty and a mask. An array
            requires use of the ``variance``, ``mask``, & ``unit`` arguments.

        trace_object : `~specreduce.tracing.Trace`, required
            The associated 1D trace object created for the 2D image.

        disp_axis : int, optional
            The index of the image's dispersion axis.

        crossdisp_axis : int, optional
            The index of the image's cross-dispersion axis.

        bkgrd_prof : `~astropy.modeling.Model`, optional
            A model for the image's background flux.

        spatial_profile : str or dict, optional
            The shape of the object profile. The first option is 'gaussian' to fit
            a uniform 1D gaussian to the average of pixels in the cross-dispersion
            direction. The other option is 'interpolated_profile'  - when this
            option is used, the profile is sampled in bins and these samples are
            interpolated between to construct a continuously varying, empirical
            spatial profile for extraction. For this option, if passed in as a
            string (i.e spatial_profile='interpolated_profile') the default values
            for the number of bins used (10) and degree of interpolation
            (linear in x and y, by default) will be used. To set these parameters,
            pass in a dictionary with the keys 'n_bins_interpolated_profile' (which
            accepts an integer number of bins) and 'interp_degree' (which accepts an
            int, or tuple of ints for x and y degree, respectively).
            [default: gaussian]

        variance : `~numpy.ndarray`, optional
            (Only used if ``image`` is not an NDData object.)
            The associated variances for each pixel in the image. Must
            have the same dimensions as ``image``. If all zeros, the variance
            will be ignored and treated as all ones.  If any zeros, those
            elements will be excluded via masking.  If any negative values,
            an error will be raised.

        mask : `~numpy.ndarray`, optional
            (Only used if ``image`` is not an NDData object.)
            Whether to mask each pixel in the image. Must have the same
            dimensions as ``image``. If blank, all non-NaN pixels are
            unmasked.

        unit : `~astropy.units.Unit` or str, optional
            (Only used if ``image`` is not an NDData object.)
            The associated unit for the data in ``image``. If blank,
            fluxes are interpreted in DN.


        Returns
        -------
        spec_1d : `~specutils.Spectrum1D`
            The final, Horne extracted 1D spectrum.
        """
        image = image if image is not None else self.image
        trace_object = trace_object if trace_object is not None else self.trace_object
        disp_axis = disp_axis if disp_axis is not None else self.disp_axis
        crossdisp_axis = crossdisp_axis if crossdisp_axis is not None else self.crossdisp_axis
        bkgrd_prof = bkgrd_prof if bkgrd_prof is not None else self.bkgrd_prof
        spatial_profile = (spatial_profile if spatial_profile is not None else
                           self.spatial_profile)
        variance = variance if variance is not None else self.variance
        mask = mask if mask is not None else self.mask
        unit = unit if unit is not None else self.unit

        # figure out what 'spatial_profile' was provided
        # put this parsing into another method at some point, its a lot..
        interp_degree_interpolated_profile = None
        n_bins_interpolated_profile = None
        spatial_profile_choices = ('gaussian', 'interpolated_profile')

        if isinstance(spatial_profile, str):
            spatial_profile = spatial_profile.lower()
            if spatial_profile not in spatial_profile_choices:
                raise ValueError("spatial_profile must be one of"
                                 f"{', '.join(spatial_profile_choices)}")
            if spatial_profile == 'interpolated_profile':  # use defaults
                bkgrd_prof = None
                n_bins_interpolated_profile = 10
                interp_degree_interpolated_profile = 1
        elif isinstance(spatial_profile, dict):
            # first, figure out what type of profile is indicated
            # right now, the only type that should use a dictionary is 'interpolated_profile'
            # but this may be extended in the future hence the 'name' key. also,
            # the gaussian option could be supplied as a single-key dict

            # will raise key error if not present, and also fail on .lower if not
            # a string - think this is informative enough
            spatial_profile_type = spatial_profile['name'].lower()

            if spatial_profile_type not in spatial_profile_choices:
                raise ValueError("spatial_profile must be one of"
                                 f"{', '.join(spatial_profile_choices)}")
            if spatial_profile_type == 'gaussian':
                spatial_profile = 'gaussian'
            else:
                if 'n_bins_interpolated_profile' in spatial_profile.keys():
                    n_bins_interpolated_profile = \
                        spatial_profile['n_bins_interpolated_profile']
                else:  # use default
                    n_bins_interpolated_profile = 10
                if 'interp_degree_interpolated_profile' in spatial_profile.keys():
                    interp_degree_interpolated_profile = \
                        spatial_profile['interp_degree_interpolated_profile']
                else:  # use default
                    interp_degree_interpolated_profile = 1
            spatial_profile = spatial_profile_type
            bkgrd_prof = None
        else:
            raise ValueError('``spatial_profile`` must either be string or dictionary.')

        # parse image and replace optional arguments with updated values
        self.image = self._parse_image(image, variance, mask, unit, disp_axis)
        variance = self.image.uncertainty.array
        mask = self.image.mask
        unit = self.image.unit

        img = np.ma.masked_array(self.image.data, mask=mask)

        # create separate mask including any previously uncaught non-finite
        # values for purposes of calculating fit
        or_mask = np.logical_or(img.mask,
                                ~np.isfinite(self.image.data))

        # If the trace is not flat, shift the rows in each column
        # so the image is aligned along the trace:
        if not isinstance(trace_object, FlatTrace):
            img = _align_along_trace(
                img,
                trace_object.trace,
                disp_axis=disp_axis,
                crossdisp_axis=crossdisp_axis)

        if self.spatial_profile == 'gaussian':

            # fit profile to average (mean) profile along crossdisp axis
            fit_ext_kernel = self._fit_gaussian_spatial_profile(img,
                                                                disp_axis,
                                                                crossdisp_axis,
                                                                or_mask,
                                                                bkgrd_prof)

            # this is just creating an array of the trace to shift the mean
            # when iterating over each wavelength. this needs to be fixed in the
            # future to actually account for the trace shape in a non-flat trace
            # (or possibly omitted all togehter as it might be redundant if
            # _align_along_trace is correcting this already)
            if isinstance(trace_object, FlatTrace):
                mean_init_guess = trace_object.trace
            else:
                mean_init_guess = np.broadcast_to(
                    img.shape[crossdisp_axis] // 2, img.shape[disp_axis]
                )

        else:  # interpolated_profile
            # for now, bkgrd_prof must be None because a compound model can't
            # be created with a interpolator + model. i think there is a way
            # around this, but will follow up later
            if bkgrd_prof is not None:
                raise ValueError('When `spatial_profile`is `interpolated_profile`,'
                                 '`bkgrd_prof` must be None. Background should'
                                 ' be fit and subtracted from `img` beforehand.')
            # make sure n_bins doesnt exceed the number of (for now) finite
            # columns. update this when masking is fixed.
            n_finite_cols = np.logical_or.reduce(or_mask, axis=crossdisp_axis)
            n_finite_cols = np.count_nonzero(n_finite_cols.astype(int) == 0)

            # determine interpolation degree from input and make tuple if int
            # this can also be moved to another method to parse the input
            # 'spatial_profile' arg, eventually
            if isinstance(interp_degree_interpolated_profile, int):
                kx = ky = interp_degree_interpolated_profile
            else:  # if input is tuple of ints

                if not isinstance(interp_degree_interpolated_profile, tuple):
                    raise ValueError("``interp_degree_interpolated_profile`` must be ",
                                     "an integer or tuple of integers.")
                if not all(isinstance(x, int) for x in interp_degree_interpolated_profile):
                    raise ValueError("``interp_degree_interpolated_profile`` must be ",
                                     "an integer or tuple of integers.")

                kx, ky = interp_degree_interpolated_profile

            if n_bins_interpolated_profile >= n_finite_cols:

                raise ValueError(f'`n_bins_interpolated_profile` ({n_bins_interpolated_profile}) '
                                 'must be less than the number of fully-finite '
                                 f'wavelength columns ({n_finite_cols}).')

            interp_spatial_prof = self._fit_self_spatial_profile(img, disp_axis,
                                                                 crossdisp_axis,
                                                                 or_mask,
                                                                 n_bins_interpolated_profile,
                                                                 kx, ky)

            # add private attribute to save fit profile. should this be public?
            self._interp_spatial_prof = interp_spatial_prof

        col_mask = np.logical_or.reduce(or_mask, axis=crossdisp_axis)
        nonf_col = [np.nan] * img.shape[crossdisp_axis]

        # array of 'x' values for each wavelength for extraction
        nrows = img.shape[crossdisp_axis]
        xd_pixels = np.arange(nrows)

        kernel_vals = []
        norms = []
        for col_pix in range(img.shape[disp_axis]):

            # for now, skip columns with any non-finite values
            # NOTE: fit and other kernel operations should support masking again
            # once a fix is in for renormalizing columns with non-finite values
            if col_mask[col_pix]:
                kernel_vals.append(nonf_col)
                norms.append(np.nan)
                continue

            if self.spatial_profile == 'gaussian':

                # set compound model's mean to column's matching trace value
                # again, this is probably not necessary
                if bkgrd_prof is not None:  # attr names will be diff. if not compound
                    fit_ext_kernel.mean_0 = mean_init_guess[col_pix]
                else:
                    fit_ext_kernel.mean = mean_init_guess[col_pix]

                # evaluate fit model (with shifted mean, based on trace)
                fitted_col = fit_ext_kernel(xd_pixels)

                # save result and normalization
                # this doesn't need to be in this loop, address later
                kernel_vals.append(fitted_col)

                norms.append(fit_ext_kernel.amplitude_0
                             * fit_ext_kernel.stddev_0 * np.sqrt(2*np.pi))

            else:  # interpolated_profile
                fitted_col = interp_spatial_prof(col_pix, xd_pixels)
                kernel_vals.append(fitted_col)
                norms.append(trapezoid(fitted_col, dx=1)[0])

        # transform fit-specific information
        kernel_vals = np.vstack(kernel_vals).T
        norms = np.array(norms)

        # calculate kernel normalization
        g_x = np.sum(kernel_vals**2 / variance, axis=crossdisp_axis)

        # sum by column weights
        weighted_img = np.divide(img * kernel_vals, variance)
        result = np.sum(weighted_img, axis=crossdisp_axis) / g_x

        # multiply kernel normalization into the extracted signal
        extraction = result * norms
        extracted_variance = np.sum(np.abs(kernel_vals), axis=crossdisp_axis) / g_x * norms**2
        extracted_stddev = np.sqrt(extracted_variance)
        self.sp_uncertainty = extracted_stddev

        # convert the extraction to a Spectrum1D object
        return Spectrum1D(extraction * unit,
                          spectral_axis=self.image.spectral_axis)
        
        
def _align_along_trace(img, trace_array, disp_axis=1, crossdisp_axis=0):
    """
    Given an arbitrary trace ``trace_array`` (an np.ndarray), roll
    all columns of ``nddata`` to shift the NDData's pixels nearest
    to the trace to the center of the spatial dimension of the
    NDData.
    """
    # TODO: this workflow does not support extraction for >2D spectra
    if not (disp_axis == 1 and crossdisp_axis == 0):
        # take the transpose to ensure the rows are the cross-disp axis:
        img = img.T

    n_rows, n_cols = img.shape

    # indices of all columns, in their original order
    rows = np.broadcast_to(np.arange(n_rows)[:, None], img.shape)
    cols = np.broadcast_to(np.arange(n_cols), img.shape)

    # we want to "roll" each column so that the trace sits in
    # the central row of the final image
    shifts = trace_array.astype(int) - n_rows // 2

    # we wrap the indices so we don't index out of bounds
    shifted_rows = np.mod(rows + shifts[None, :], n_rows)

    return img[shifted_rows, cols]