import pandas as pd
import sys
import numpy as np
from gwcs import coordinate_frames as cf
from gwcs import wcs
from astropy import units as u
# from specbox.basemodule import SpecIRAF
from specutils import Spectrum1D
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from specutils.manipulation import FluxConservingResampler, LinearInterpolatedResampler
from .utils import linear_fit_poly1D


class WaveCalibrator:
    """
    Wavelength calibration class for a given set of pixel-wavelength pairs.
    """
    def __init__(self, line_pixels, line_wavelengths, wave_unit='Angstrom',
                 model='Legendre1D', degree=7, sigma=2, plot=False):
        """
        Parameters
        ----------
        line_pixels : array-like
            Pixel positions of the spectral lines.
        line_wavelengths : array-like
            Wavelengths of the spectral lines.
        wave_unit : str, optional
            Unit of the wavelength, by default 'Angstrom'.
        model : str, optional
            Model to fit the pixel-wavelength pairs, by default 'Legendre1D'.
        degree : int, optional
            Degree of the polynomial model, by default 7.
        sigma : float, optional
            Sigma value for sigma clipping, by default 2.
        plot : bool, optional
            Plot the fitting results, by default False.
        """
        self.line_pixels = line_pixels
        self.line_wavelengths = line_wavelengths
        self.wave_unit = getattr(u, wave_unit)
        self.model = model
        self.degree = degree
        self.sigma = sigma
        self.plot = plot
        self.calibrate()
        
    def calibrate(self):
        """
        Calibrate the pixel-wavelength pairs using the specified model and
        construct a wcs object.
        """
        self.model_fit = linear_fit_poly1D(
            self.line_pixels, self.line_wavelengths, model=self.model,
            degree=self.degree, sigma=self.sigma, plot=self.plot)
        # Create a WCS object
        pixel_frame = cf.CoordinateFrame(1, "SPECTRAL", [0,], axes_names=["x",], unit=[u.pix,])
        spectral_frame = cf.SpectralFrame(axes_names=["wavelength",], unit=[self.wave_unit,])
        pipeline = [(pixel_frame, self.model_fit), (spectral_frame, None)]
        wcsobj = wcs.WCS(pipeline)
        self.wcs = wcsobj
        
    def apply_to_spectrum(self, spectrum, resample=False):
        """
        Apply the wavelength calibration to a given spectrum.
        
        Parameters
        ----------
        spectrum : Spectrum1D
            Spectrum to be calibrated.

        Returns
        -------
        updated_spectrum : Spectrum1D
            Calibrated spectrum.
        """
        updated_spectrum = Spectrum1D(spectrum.flux, wcs=self.wcs, mask=spectrum.mask,
                                      uncertainty=spectrum.uncertainty)
        if updated_spectrum.spectral_axis[0] > updated_spectrum.spectral_axis[-1]:
            new_wave = updated_spectrum.spectral_axis[::-1]
            new_flux = updated_spectrum.flux[::-1]
            if updated_spectrum.uncertainty is not None:
                updated_spectrum.uncertainty.array = updated_spectrum.uncertainty.array[::-1]
        else:
            new_wave = updated_spectrum.spectral_axis
            new_flux = updated_spectrum.flux    
        new_spec = Spectrum1D(
            flux=new_flux, 
            spectral_axis=new_wave, 
            uncertainty=updated_spectrum.uncertainty)
        return new_spec
        

def read_identify_feature_list(filename):
    """
    Read the feature list from the output file of IRAF identify task.
    
    Parameters
    ----------
    filename : str
        Name of the IRAF identified feature list file.
        
    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing the feature list.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Find all start indices of the feature lists
    start_indices = []
    for i, line in enumerate(lines):
        if 'features' in line:
            start_indices.append(i + 1)
    if not start_indices:
        raise ValueError("Feature list not found in the file.")
    # Use the last start index to parse the latest feature list
    start_index = start_indices[-1]
    feature_lines = []
    for line in lines[start_index:]:
        feature_lines.append(line.strip())
    # Convert feature lines to DataFrame
    columns = ['Position', 'ObservedWavelength', 
               'ReferenceWavelength', 'Width', 
               'colx1', 'colx2', 'LineName']
    data = []
    for line in feature_lines:
        parts = line.split()
        if len(parts) == 7:  # Ensure line has the correct number of parts
            data.append(parts)
        else:
            print(f"Skipped line: {line}")
    df = pd.DataFrame(data, columns=columns)
    numeric_columns = ['Position', 'ObservedWavelength', 
                       'ReferenceWavelength']
    df[numeric_columns] = df[numeric_columns].astype(float)
    df = df.loc[:, ['Position', 'ObservedWavelength', 
                    'ReferenceWavelength', 'LineName']]  
    return df

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean)**2 / (2 * sigma**2))

def subpix_peaks(lamp_spectrum_flux, height=4000):
    peaks, _ = find_peaks(lamp_spectrum_flux, height=height)  # Define an appropriate threshold
    refined_pixel_positions = []
    for peak in peaks:
        window_size = 5  # Adjust based on your data
        x_data = np.arange(max(0, peak - window_size), 
                           min(len(lamp_spectrum_flux), 
                               peak + window_size + 1))
        y_data = lamp_spectrum_flux[x_data]
    # Initial guesses: amplitude, mean (in pixel units), standard deviation
        p0 = [lamp_spectrum_flux[peak], peak, 1]
        try:
            popt, _ = curve_fit(gaussian, x_data, y_data, p0=p0)
            refined_peak_pixel_position = popt[1]  # The refined mean position in pixel units
            refined_pixel_positions.append(refined_peak_pixel_position)
        except RuntimeError:
            print(f"Could not fit a Gaussian to the peak at index {peak}")
    return refined_pixel_positions
