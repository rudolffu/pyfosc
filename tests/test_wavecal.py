import numpy as np
import pytest

from pyfosc import wavecal as wc
from pyfosc.specio import Spec1D


def _templates():
    out = []
    for tel in ("XLT", "LJT", "HCT"):
        out.extend((tel, t) for t in wc.list_templates(tel))
    return out


TEMPLATES = _templates()


def test_templates_are_packaged():
    assert len(TEMPLATES) >= 10
    assert {t.kind for _, t in TEMPLATES} == {"fear", "hene"}


def test_recenter_flat_topped_line():
    x = np.arange(200.0)
    flux = 100.0 + 5000.0 * ((x > 95.3) & (x < 104.3))  # slit image, 9 px wide
    flux += 3000.0 * np.exp(-0.5 * ((x - 40.6) / 1.5) ** 2)
    centres = wc.recenter_lines(flux, [98.0, 42.0, 150.0], fwidth=10)
    assert abs(centres[0] - 100.0) < 0.6
    assert abs(centres[1] - 40.6) < 0.1
    assert np.isnan(centres[2])


@pytest.mark.parametrize("tel,tmpl", TEMPLATES, ids=[t.name for _, t in TEMPLATES])
def test_shifted_template_is_recovered(tel, tmpl):
    rng = np.random.default_rng(42)
    n = len(tmpl.flux)
    x = np.arange(n, dtype=float)
    shift = 17.3
    lamp = 0.6 * np.interp(x - shift, x, tmpl.flux)
    lamp += rng.normal(0, 3e-4 * np.nanmax(tmpl.flux), n)

    match = wc.match_template(lamp, tmpl.flux)
    assert not match["flipped"]
    assert abs(match["shift"] - shift) < 3

    sol = wc.solve_from_template(lamp, tmpl)
    assert sol is not None and sol.is_good(), sol.summary() if sol else "no solution"
    # lines must sit where the template had them, moved by the shift
    ref = wc.recenter_lines(tmpl.flux, tmpl.pixels)
    offsets = []
    for pix, wave in zip(sol.pixels[sol.used], sol.wavelengths[sol.used]):
        j = int(np.argmin(np.abs(tmpl.wavelengths - wave)))
        if np.isfinite(ref[j]):
            offsets.append(pix - (ref[j] + shift))
    offsets = np.abs(offsets)
    assert len(offsets) >= 8
    assert np.median(offsets) < 0.2
    assert np.percentile(offsets, 90) < 0.5


def test_flipped_lamp_is_detected():
    tmpl = TEMPLATES[0][1]
    sol = wc.solve_from_template(tmpl.flux[::-1].copy(), tmpl)
    direct = wc.solve_from_template(tmpl.flux, tmpl)
    assert sol is not None and direct is not None
    assert np.sign(sol.dispersion) == -np.sign(direct.dispersion)
    assert abs(abs(sol.dispersion) - abs(direct.dispersion)) < 0.02


def test_auto_identify_picks_a_good_solution():
    tel, tmpl = next((tel, t) for tel, t in TEMPLATES if tel == "LJT" and "G3_202212" in t.name)
    sol = wc.auto_identify(tmpl.flux, wc.list_templates(tel, kind="hene", grism="G3"))
    assert sol is not None and sol.is_good()
    assert "G3" in sol.template


def test_solution_json_roundtrip(tmp_path):
    tmpl = TEMPLATES[0][1]
    sol = wc.solve_from_template(tmpl.flux, tmpl)
    sol.lamp = "aflamp.fits"
    loaded = wc.WaveSolution.load(sol.save(tmp_path / "db" / "wavesol.json"))
    grid = np.arange(sol.npix)
    np.testing.assert_allclose(loaded(grid), sol(grid))
    assert loaded.n_used == sol.n_used and loaded.lamp == "aflamp.fits"


def test_linelists():
    hene, _ = wc.load_linelist("hene")
    fear, _ = wc.load_linelist("fear")
    assert len(hene) > 100 and len(fear) > 30
    assert np.all(np.diff(hene) > 0) and np.all(np.diff(fear) > 0)


def test_add_lines_from_list():
    tel, tmpl = next((tel, t) for tel, t in TEMPLATES if "slit1.8_G8" in t.name)
    sol = wc.solve_from_template(tmpl.flux, tmpl)
    waves, names = wc.load_linelist("hene")
    pix, wav, _ = wc.add_lines_from_list(tmpl.flux, sol, waves, names)
    assert len(pix) >= sol.n_used
    refit = wc.fit_solution(pix, wav, npix=sol.npix)
    assert refit.rms_pix < 0.7


@pytest.mark.parametrize("sign", [1, -1])
def test_dispcor_conserves_flux(sign):
    n = 500
    x = np.arange(n, dtype=float)
    wave = 6000 + sign * (2.0 * (x - n / 2) + 2e-4 * (x - n / 2) ** 2)
    sol = wc.fit_solution(x[::25], wave[::25], npix=n, degree=3)
    flux = 100 + 50 * np.sin(x / 30.0)
    line_pix = 180.0
    flux += 2000 * np.exp(-0.5 * ((x - line_pix) / 2.0) ** 2)
    spec = Spec1D(bands=[flux, flux, np.ones(n), np.full(n, 3.0)])
    out = wc.dispcor(spec, sol)
    assert out.calibrated and out.cdelt > 0
    assert out.npix == n
    np.testing.assert_allclose(out.flux.sum(), flux.sum(), rtol=2e-3)
    # the emission line lands at its wavelength
    peak_wave = out.wave[np.argmax(out.flux)]
    assert abs(peak_wave - sol(line_pix)) < abs(out.cdelt)
    # uniform noise stays close to its input level
    assert abs(np.median(out.sigma) - 3.0) < 0.3
