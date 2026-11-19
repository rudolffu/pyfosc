import numpy as np
import pytest
from astropy.io import fits

from pyfosc import fluxcal as fc
from pyfosc.specio import Spec1D


@pytest.mark.parametrize("obj,expected", [
    ("BD+28 4211", "bd284211"),
    ("BD+28d4211", "bd284211"),
    ("Feige 34", "feige34"),
    ("feige34_std", "feige34"),
    ("HZ 44", "hz44"),
    ("G191-B2B", "g191b2b"),
    ("Hiltner 600", "hilt600"),
])
def test_find_standard(obj, expected):
    entry = fc.find_standard(obj)
    assert entry is not None and entry[0] == expected


@pytest.mark.parametrize("obj", ["PG 2130+099", "J0141+7307", "", "FeAr"])
def test_find_standard_rejects_targets(obj):
    assert fc.find_standard(obj) is None


def test_standard_flux_scale():
    # Feige 34 has V ~ 11.2: f_lambda(5500 A) ~ 1.2e-13 erg/cm2/s/A
    tab = fc.load_standard("feige34", caldir="spec50cal")
    flam = np.interp(5500.0, tab.wave, tab.flam)
    assert 0.8e-13 < flam < 1.6e-13


def test_extinction_curves():
    for tel in ("XLT", "LJT", "HCT"):
        ext = fc.load_extinction(tel)
        assert ext(4000.0) > ext(7000.0) > 0


def test_airmass_from_header_and_computed():
    assert fc.get_airmass(fits.Header({"AIRMASS": 1.37})) == pytest.approx(1.37)
    hdr = fits.Header({"RA": "10:39:36.7", "DEC": "+43:06:09", "EXPTIME": 60.0,
                       "DATE-OBS": "2023-01-18T18:00:00"})
    value = fc.get_airmass(hdr, telescope="XLT")
    assert 1.0 <= value < 5.0
    with pytest.warns(UserWarning):
        assert fc.get_airmass(fits.Header({"OBJECT": "x"})) == 1.0


def _fake_observation(table, ext, airmass, exptime, sens_true, wave):
    """Counts per pixel for a star observed through a known sensitivity."""
    dw = wave[1] - wave[0]
    flam = np.interp(wave, table.wave, table.flam)
    rate = flam * 10 ** (0.4 * (sens_true(wave) - airmass * ext(wave)))
    return rate * exptime * dw


@pytest.mark.parametrize("caldir", ["spec50cal", "oke1990"])
def test_sensfunc_and_calibrate_recover_standard(caldir):
    table = fc.load_standard("feige34", caldir=caldir)
    ext = fc.load_extinction("XLT")
    wave = 3900.0 + 2.8 * np.arange(1700)

    def sens_true(w):
        return 38.0 - 2.0 * ((w - 6000.0) / 2500.0) ** 2

    counts = _fake_observation(table, ext, 1.3, 120.0, sens_true, wave)
    std = Spec1D(bands=[counts, counts, counts * 0, np.sqrt(counts)],
                 crval=wave[0], cdelt=2.8, calibrated=True)
    sf = fc.sensfunc(std, table, ext, airmass=1.3, exptime=120.0)
    inner = (wave > 4200) & (wave < 8300)
    assert np.max(np.abs(sf(wave[inner]) - sens_true(wave[inner]))) < 0.05

    # a second "target" at another airmass and exposure time
    counts2 = _fake_observation(table, ext, 1.8, 600.0, sens_true, wave)
    obj = Spec1D(bands=[counts2, counts2, counts2 * 0, np.sqrt(counts2)],
                 crval=wave[0], cdelt=2.8, calibrated=True)
    cal = fc.calibrate(obj, sf.wave, sf.sens, ext, airmass=1.8, exptime=600.0)
    truth = np.interp(wave, table.wave, table.flam)
    ratio = cal.flux[inner] / truth[inner]
    assert abs(np.median(ratio) - 1.0) < 0.02
    assert np.percentile(np.abs(ratio - 1.0), 95) < 0.05
    assert cal.header["BUNIT"] == "erg/cm2/s/A"


def test_sensfunc_ignores_telluric_absorption():
    table = fc.load_standard("feige34", caldir="spec50cal")
    ext = fc.load_extinction("LJT")
    wave = 5000.0 + 2.0 * np.arange(2200)
    counts = _fake_observation(table, ext, 1.1, 100.0, lambda w: 37.0 + 0 * w, wave)
    dip = 1 - 0.6 * np.exp(-0.5 * ((wave - 7630.0) / 20.0) ** 2)
    std = Spec1D(bands=[counts * dip], crval=wave[0], cdelt=2.0, calibrated=True)
    sf = fc.sensfunc(std, table, ext, airmass=1.1, exptime=100.0)
    assert abs(sf(7630.0) - 37.0) < 0.05


def test_zero_point_matches_iraf_standard():
    # IRAF "standard" converts magnitudes with fnu(0) = 3.68e-20 erg/cm2/s/Hz;
    # an IRAF std file lists 7.3732e-13 for BD+28 4211 (subaru) at 4048.6 A.
    table = fc.load_standard("bd284211", "subaru")
    i = int(np.argmin(np.abs(np.asarray(table.wave) - 4048.6)))
    assert abs(np.asarray(table.flam)[i] / 7.3732e-13 - 1) < 1e-3


def test_odd_standards_flags_the_wrong_shape():
    wave = np.linspace(4200, 8500, 60)
    real = 36.0 - 1.5 * ((wave - 6300.0) / 2500.0) ** 2
    grey = real - 0.4                                  # thin cloud: same shape
    red = real + 2.5 * (wave - 4200.0) / 4300.0        # not the star it is labelled as
    assert fc.odd_standards([(wave, real)]) == [False]
    assert fc.odd_standards([(wave, real), (wave, grey)]) == [False, False]
    assert fc.odd_standards([(wave, real), (wave, red)]) == [False, True]
    assert fc.odd_standards([(wave, red), (wave, real), (wave, grey)]) == [True, False, False]
