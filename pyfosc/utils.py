import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip

def linear_fit_poly1D(x_obs, y_obs, model='Legendre1D', degree=7,
                      sigma=2, plot=True):
    # check if x_obs and y_obs are numpy arrays
    if not isinstance(x_obs, np.ndarray):
        x_obs = np.array(x_obs)
    if not isinstance(y_obs, np.ndarray):
        y_obs = np.array(y_obs)
    if len(x_obs) != len(y_obs):
        raise ValueError('x_obs and y_obs must have the same length')
    if len(x_obs) < 2:
        raise ValueError('x_obs and y_obs must have at least 2 elements')
    # check if model is a valid model
    if model not in ['Legendre1D', 'Chebyshev1D', 'Polynomial1D']:
        raise ValueError('model must be one of "Legendre1D", "Chebyshev1D", "Polynomial1D"')
    fitter = fitting.FittingWithOutlierRemoval(
        fitting.LinearLSQFitter(), sigma_clip, 
        niter=3, sigma=sigma)
    model_init = getattr(models, model)(degree)
    model_fit, mask = fitter(model_init, x=x_obs, y=y_obs)
    x_fit = np.linspace(x_obs.min(), x_obs.max(), 100)
    y_fit = model_fit(x_fit)
    x_ma = np.ma.masked_array(x_obs, mask=mask)
    y_ma = np.ma.masked_array(y_obs, mask=mask)
    y_pred = model_fit(x_ma)
    residuals = y_ma - y_pred
    rms = np.sqrt(np.mean(residuals**2))
    print(
        f"Model: {model} of degree {degree}.\n"
        f"Number of points used in the fitting: {x_ma.count()}.\n"
        f"RMS of the fitting: {rms:.4f}.")
    if plot:
        plot_linear_fit(x_obs, y_obs, x_fit, y_fit, residuals)
    return model_fit
    
def plot_linear_fit(x_obs, y_obs, x_fit, y_fit, residuals):
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                         gridspec_kw={'height_ratios': [2, 1],
                                      'hspace': 0})
    axes[0].plot(x_obs, y_obs, 'o', label='Observed Data')
    axes[0].plot(x_fit, y_fit, color='red', label='fit')
    axes[0].legend()
    axes[1].scatter(x_obs, residuals, color='black', label='residuals')
    axes[1].axhline(y=0, color='r', linestyle='--')
    axes[1].legend()
    plt.show()