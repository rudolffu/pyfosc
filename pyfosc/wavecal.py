import pandas as pd
import sys
import numpy as np
from gwcs import coordinate_frames as cf
from gwcs import wcs
import pyqtgraph as pg
from astropy import units as u
from PyQt5.QtGui import QCursor, QFont 
from PyQt5.QtCore import Qt, QThread
from PyQt5.QtWidgets import QApplication, QFrame, QWidget, QInputDialog
from pyqtgraph.dockarea import DockArea, Dock
from specbox.basemodule import SpecIRAF
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

class InteractiveSpectrumPlot(pg.PlotWidget):
    def __init__(self, spectrum_x, spectrum_y, 
                 xinwave=False, fit_peaks=True,
                 df=None,
                 *args, **kwargs):
        super().__init__()
        self.spectrum_x = spectrum_x
        self.spectrum_y = spectrum_y
        self.xinwave = xinwave
        self.fit_peaks = fit_peaks
        if self.fit_peaks==True:
            refined_pixel_positions = subpix_peaks(
                self.spectrum_y, height=5000)
            df = pd.DataFrame({'Position': refined_pixel_positions})
            df['UserWavelength'] = np.nan
        self.df = df
        self.setBackground('w')
        self.plotSpectrum()
        self.showGrid(x=True, y=True)
        self.setMouseEnabled(x=True, y=True)
        self.setLogMode(x=False, y=False)
        self.setAspectLocked(False)
        # Show auto-range button
        self.enableAutoRange()
        self.vb = self.getViewBox()
        # Enable Mouse selection for zooming
        self.vb.setMouseMode(self.vb.RectMode)
        # Connect the mouse click event to the onClick method
        self.scene().sigMouseClicked.connect(self.onClick)
    
    def plotSpectrum(self):
        self.plot(self.spectrum_x, self.spectrum_y, 
                  pen=pg.mkPen('b', width=3),)
        if self.df is not None:
            if 'ReferenceWavelength' in self.df.columns:
                for peak in self.df['ReferenceWavelength']:
                    self.plot([peak, peak], [0, self.spectrum_y.max()], pen=pg.mkPen('m', width=1))
            elif 'Position' in self.df.columns:
                for peak in self.df['Position']:
                    self.plot([peak, peak], [0, self.spectrum_y.max()], pen=pg.mkPen('r', width=1))

    def onClick(self, event):
        # Map the click position to the plot coordinates
        pos = event.scenePos()
        if self.sceneBoundingRect().contains(pos):
            mappedPos = self.plotItem.vb.mapSceneToView(pos)
            x = mappedPos.x()
            y = mappedPos.y()
            print(f"Clicked position: x={x:.2f}, y={y:.2f}")
            # Here you can implement the logic to mark and label the spectral line
    
    def keyPressEvent(self, event):
        if event.modifiers() & Qt.ControlModifier:
            if event.key() == Qt.Key_F: # Cmd + F, flip x-axis
                self.invertX()
            elif event.key() == Qt.Key_R: # Cmd + R, reset view
                self.clear()
                self.plotSpectrum()
            elif event.key() == Qt.Key_L: # Cmd + L, toggle log scale
                self.setLogMode(x=False, y=True)
            elif event.key() == Qt.Key_S:
                self.df.to_csv('ref_wave.csv', index=False)
                print("Saved to ref_wave.csv")
        if event.key() == Qt.Key_Space:
            mouse_pos = self.mapFromGlobal(QCursor.pos())
            self.vb = self.getViewBox()
            x_now = self.vb.mapSceneToView(mouse_pos).x()
            if self.xinwave==True:
                idx = np.abs(self.df['ReferenceWavelength'].values - x_now).argmin()
                wave = self.df['ReferenceWavelength'].values[idx]
                print(f"Nearest line: x_wave={wave:.4f}, line={self.df['LineName'].values[idx]}")
            else:
                idx = np.abs(self.df['Position'].values - x_now).argmin()
                x_pix = self.df['Position'].values[idx]
                print(f"Nearest peak: x_pix={x_pix:.3f}, ref_wave={self.df['UserWavelength'].values[idx]:.4f}")
        if event.key() == Qt.Key_I:
            mouse_pos = self.mapFromGlobal(QCursor.pos())
            self.vb = self.getViewBox()
            x_now = self.vb.mapSceneToView(mouse_pos).x()
            idx = np.abs(self.df['Position'].values - x_now).argmin()
            x_pix = self.df['Position'].values[idx]
            # take user input as float
            user_wave, ok = QInputDialog.getDouble(self.parent(), "User Input", "Enter the wavelength:", decimals=4)
            if ok:
                self.df.loc[idx, 'UserWavelength'] = user_wave
                print(f"User input: x_pix={x_pix:.3f}, ref_wave={user_wave:.4f}")
        

class SpectrumViewerApp(QApplication):
    def __init__(self, spectrum_x, spectrum_y, 
                 ref_x=None, ref_y=None,
                 ref_df=None,
                 *args, **kwargs):
        super().__init__(sys.argv)
        self.spectrum_x = spectrum_x
        self.spectrum_y = spectrum_y
        self.ref_x = ref_x
        self.ref_y = ref_y
        self.ref_df = ref_df
        self.plot = InteractiveSpectrumPlot(
            self.spectrum_x, self.spectrum_y, fit_peaks=True,
            xinwave=False, *args, **kwargs)
        self.make_layout()
        self.exec_()
        self.exit() 
        sys.exit()  
    
    def make_layout(self):
        layout = pg.LayoutWidget()
        layout.resize(1200, 800)
        layout.setWindowTitle("PyFOSC Spectrum Viewer")
        
        # Create a dock area
        dock_area = DockArea()
        
        # Create a dock for the main spectrum plot
        main_dock = Dock("Main Spectrum Plot", size=(800, 600))
        main_dock.addWidget(self.plot)
        dock_area.addDock(main_dock)
        if self.ref_x is not None and self.ref_y is not None:
            ref_spectrum_x = self.ref_x
            ref_spectrum_y = self.ref_y
        # Create a dock for the reference spectrum plot
        ref_dock = Dock("Reference Spectrum Plot", size=(800, 400))
        ref_plot = InteractiveSpectrumPlot(ref_spectrum_x, ref_spectrum_y, fit_peaks=False,
                                           xinwave=True, df=self.ref_df)
        for idx, row in self.ref_df.iterrows():
                ref_plot.plot([row['ReferenceWavelength'], row['ReferenceWavelength']], 
                              [0, ref_spectrum_y.max()], pen=pg.mkPen('m', width=1))
        ref_dock.addWidget(ref_plot)
        dock_area.addDock(ref_dock, position="top", relativeTo=main_dock)
        
        layout.addWidget(dock_area)
        self.layout = layout
        self.layout.show()


class SpectrumViewerThread(QThread):
    def __init__(self, spectrum_x, spectrum_y, *args, **kwargs):
        super().__init__()
        self.spectrum_x = spectrum_x
        self.spectrum_y = spectrum_y
        self.app = SpectrumViewerApp(self.spectrum_x, self.spectrum_y, *args, **kwargs)
        self.app.aboutToQuit.connect(self.exit)

    def run(self):
        self.app.exec_()
        self.exit()
        self.finished.emit()
        sys.exit()
        