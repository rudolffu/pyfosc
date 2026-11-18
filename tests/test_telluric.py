import numpy as np

from pyfosc import telluric as tl
from pyfosc.specio import Spec1D


def _absorption(wave, depth_o2=1.0, depth_h2o=1.0):
    """Optical depth of toy O2 A/B and H2O bands."""
    tau = depth_o2 * (0.9 * np.exp(-0.5 * ((wave - 7630.0) / 18.0) ** 2)
                      + 0.35 * np.exp(-0.5 * ((wave - 6890.0) / 14.0) ** 2))
    tau += depth_h2o * 0.25 * np.exp(-0.5 * ((wave - 7250.0) / 40.0) ** 2)
    return tau


def _spec(wave, flux):
    return Spec1D(bands=[flux, flux, np.ones_like(flux), 0.01 * flux],
                  crval=wave[0], cdelt=wave[1] - wave[0], calibrated=True)


WAVE = 5500.0 + 2.5 * np.arange(1000)   # 5500-8000 A
STD_CONT = 3e-13 * (WAVE / 6000.0) ** -3.5


def test_continuum_ignores_bands():
    flux = STD_CONT * np.exp(-_absorption(WAVE))
    cont = tl.fit_continuum(WAVE, flux, exclude=tl._band_mask(WAVE))
    assert np.max(np.abs(cont / STD_CONT - 1)) < 0.01


def test_transmission_is_one_outside_bands():
    std = _spec(WAVE, STD_CONT * np.exp(-_absorption(WAVE)))
    trans = tl.transmission(std, telescope="XLT")
    truth = np.exp(-_absorption(WAVE))
    outside = ~tl._band_mask(WAVE) | (WAVE < tl.BLUE_LIMIT)
    assert np.all(trans.flux[outside] == 1.0)
    a_band = (WAVE > 7600) & (WAVE < 7680)
    assert np.max(np.abs(trans.flux[a_band] - truth[a_band])) < 0.02


def test_correction_at_other_airmass_and_shift():
    rng = np.random.default_rng(3)
    std = _spec(WAVE, STD_CONT * np.exp(-_absorption(WAVE)))
    trans = tl.transmission(std, telescope="XLT")
    # target: different continuum, 1.6x the O2 column, 0.7x the water, shifted 2 A
    obj_cont = 1e-15 * (1 + 0.3 * np.sin(WAVE / 400.0))
    obj_flux = obj_cont * np.exp(-_absorption(WAVE - 2.0, depth_o2=1.6, depth_h2o=0.7))
    obj_flux = obj_flux * (1 + rng.normal(0, 0.003, len(WAVE)))
    corrected, fit = tl.correct(_spec(WAVE, obj_flux), trans, airmass_ratio=1.2)

    assert abs(fit.scales["O2"] - 1.6) < 0.15
    assert abs(fit.scales["H2O"] - 0.7) < 0.15
    assert abs(fit.shift - 2.0) < 0.7
    for lo, hi in ((7590, 7700), (6860, 6960), (7150, 7350)):
        sel = (WAVE > lo) & (WAVE < hi)
        before = np.max(np.abs(obj_flux[sel] / obj_cont[sel] - 1))
        after = np.max(np.abs(corrected.flux[sel] / obj_cont[sel] - 1))
        assert after < 0.06 and after < 0.3 * before
    # untouched away from the bands, and the sky band is never modified
    clean = WAVE < 6800
    np.testing.assert_allclose(corrected.flux[clean], obj_flux[clean])
    np.testing.assert_allclose(corrected.bands[2], 1.0)


def test_no_fit_uses_airmass_ratio():
    std = _spec(WAVE, STD_CONT * np.exp(-_absorption(WAVE)))
    trans = tl.transmission(std, telescope="XLT")
    _, fit = tl.correct(std, trans, airmass_ratio=1.0, fit=False)
    assert fit.shift == 0.0 and fit.scales == {"O2": 1.0, "H2O": 1.0}
