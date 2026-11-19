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


def _template():
    std = _spec(WAVE, STD_CONT * np.exp(-_absorption(WAVE)))
    return tl.transmission(std, telescope="XLT")


def test_fitted_scale_comes_with_a_realistic_error():
    trans = _template()
    pulls = []
    for seed in range(12):
        rng = np.random.default_rng(seed)
        flux = 1e-16 * np.exp(-_absorption(WAVE, depth_o2=1.3, depth_h2o=1.3))
        flux = flux * (1 + rng.normal(0, 0.03, len(WAVE)))
        _, fit = tl.correct(_spec(WAVE, flux), trans, airmass_ratio=1.0)
        assert fit.source["O2"] == "fit"
        pulls.append((fit.scales["O2"] - 1.3) / fit.errors["O2"])
    pulls = np.array(pulls)
    assert abs(np.mean(pulls)) < 1.0 and 0.4 < np.std(pulls) < 2.0


def test_moderate_noise_still_measures_the_deep_band():
    # the airmass ratio (1.1) is wrong here; the A band is deep enough to tell
    rng = np.random.default_rng(11)
    flux = 1e-16 * np.exp(-_absorption(WAVE)) * (1 + rng.normal(0, 0.15, len(WAVE)))
    _, fit = tl.correct(_spec(WAVE, flux), _template(), airmass_ratio=1.1)
    assert fit.source["O2"] == "fit"
    assert abs(fit.scales["O2"] - 1.0) < 3 * fit.errors["O2"]
    # the shallow water band is not measurable and keeps the airmass ratio
    assert fit.source["H2O"] == "airmass" and fit.scales["H2O"] == 1.1


def test_hopeless_noise_keeps_the_fallback():
    rng = np.random.default_rng(5)
    flux = 1e-16 * np.exp(-_absorption(WAVE)) * (1 + rng.normal(0, 1.5, len(WAVE)))
    _, fit = tl.correct(_spec(WAVE, flux), _template(), airmass_ratio=1.1)
    assert fit.shift == 0.0
    assert fit.scales == {"O2": 1.1, "H2O": 1.1}
    assert set(fit.source.values()) == {"airmass"}
    # a night-level value replaces the airmass law when supplied
    _, fit = tl.correct(_spec(WAVE, flux), _template(), airmass_ratio=1.1,
                        fallback={"H2O": 1.6}, fallback_label="night")
    assert fit.scales == {"O2": 1.1, "H2O": 1.6} and fit.source["H2O"] == "night"


def test_emission_line_in_a_band_does_not_bias_the_scale():
    rng = np.random.default_rng(2)
    cont = 1e-16 * (1 + 0.2 * np.sin(WAVE / 500.0))
    line = 1.5 * np.exp(-0.5 * ((WAVE - 6900.0) / 12.0) ** 2)   # e.g. redshifted H-alpha
    flux = cont * (1 + line) * np.exp(-_absorption(WAVE, depth_o2=1.2))
    flux = flux * (1 + rng.normal(0, 0.02, len(WAVE)))
    _, fit = tl.correct(_spec(WAVE, flux), _template(), airmass_ratio=1.0)
    assert fit.source["O2"] == "fit"
    assert abs(fit.scales["O2"] - 1.2) < 0.08
    # the band holding the line is not fitted and keeps the airmass law
    assert fit.band_source["O2_6860"] == "airmass" and fit.band_scales["O2_6860"] == 1.0
    assert abs(fit.band_scales["O2_7590"] - 1.2) < 0.08


def test_bands_of_one_species_are_scaled_separately():
    # a standard whose A band is filled in (second-order light) but not its B band
    truth = np.exp(-_absorption(WAVE))
    diluted = truth.copy()
    a_band = (WAVE > 7500) & (WAVE < 7760)
    diluted[a_band] = (truth[a_band] + 0.25) / 1.25
    trans = tl.transmission(_spec(WAVE, STD_CONT * diluted), telescope="XLT")
    rng = np.random.default_rng(4)
    flux = 1e-16 * truth * (1 + rng.normal(0, 0.01, len(WAVE)))
    corrected, fit = tl.correct(_spec(WAVE, flux), trans, airmass_ratio=1.0)
    assert fit.band_scales["O2_7590"] > 1.2
    assert abs(fit.band_scales["O2_6860"] - 1.0) < 0.05
    for lo, hi in ((7590, 7700), (6860, 6960)):
        sel = (WAVE > lo) & (WAVE < hi)
        assert np.max(np.abs(corrected.flux[sel] / 1e-16 - 1)) < 0.08


def test_night_prior_pulls_noisy_targets_only():
    trans = _template()
    prior = {"O2_7590": (1.3, 0.03)}
    rng = np.random.default_rng(8)
    truth = np.exp(-_absorption(WAVE, depth_o2=1.0))
    noisy = 1e-16 * truth * (1 + rng.normal(0, 0.4, len(WAVE)))
    _, fit = tl.correct(_spec(WAVE, noisy), trans, airmass_ratio=1.0, prior=prior)
    assert abs(fit.band_scales["O2_7590"] - 1.3) < 0.05
    assert fit.band_source["O2_7590"] in ("night", "night+fit")
    bright = 1e-16 * truth * (1 + rng.normal(0, 0.005, len(WAVE)))
    _, fit = tl.correct(_spec(WAVE, bright), trans, airmass_ratio=1.0, prior=prior)
    assert abs(fit.band_scales["O2_7590"] - 1.0) < 0.03
    assert fit.band_source["O2_7590"] == "fit"


def test_template_has_local_continuum_and_drops_weak_bands():
    rng = np.random.default_rng(1)
    # curved continuum that a single power law over 6250-8000 A cannot follow
    cont = STD_CONT * (1 + 0.10 * np.sin(WAVE / 400.0))
    truth = np.exp(-_absorption(WAVE, depth_h2o=0.0))     # no water band at all
    std = _spec(WAVE, cont * truth * (1 + rng.normal(0, 0.004, len(WAVE))))
    tmpl = tl.build_template(std, telescope="XLT")
    status = {i["key"]: i["status"] for i in tmpl.info}
    assert status["O2_7590"] == "local" and status["O2_6860"] == "local"
    assert status["H2O_7150"] == "not significant"
    assert status["O2_6270"] == "blue of limit"
    for lo, hi in ((7590, 7700), (6860, 6960)):
        sel = (WAVE > lo) & (WAVE < hi)
        assert np.max(np.abs(tmpl.trans[sel] - truth[sel])) < 0.02
    assert np.all(tmpl.trans[(WAVE > 7150) & (WAVE < 7350)] == 1.0)
    # the template carries its own error inside the bands only
    err = tmpl.spec.bands[1]
    assert np.all(err[~tl._band_mask(WAVE)] == 0) and 0.002 < err[WAVE == 7630.0][0] < 0.01
