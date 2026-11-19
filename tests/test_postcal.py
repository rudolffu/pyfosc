"""End-to-end run of the post-extraction steps on synthetic spectra."""
import json
import sys

import numpy as np
import pytest
from astropy.io import fits

from pyfosc import fluxcal as fc
from pyfosc import wavecal as wc
from pyfosc.postcal import PostExtraction
from pyfosc.specio import read_spec
from pyfosc.steps import reidentlamp2m, telluric_base2m, wavecal2m

SHIFT = 6.0  # pixels between the packaged template and the "observed" lamp


def _tau(wave):
    return (0.8 * np.exp(-0.5 * ((wave - 7630.0) / 18.0) ** 2)
            + 0.3 * np.exp(-0.5 * ((wave - 6890.0) / 14.0) ** 2))


def _sens_true(wave):
    return 37.5 - 1.5 * ((wave - 6300.0) / 2500.0) ** 2


@pytest.fixture(scope="module")
def night(tmp_path_factory):
    """A fake XLT/G4 night: one lamp, one standard, one target."""
    root = tmp_path_factory.mktemp("night")
    data = root / "data"
    data.mkdir()
    settings = {"mysettings": {"telescope": "XLT", "slit": "slit1.8", "Grism": "G4", "side": ""}}
    (root / "myfosc.json").write_text(json.dumps(settings))

    tmpl = next(t for t in wc.list_templates("XLT") if "newG4" in t.name)
    n = len(tmpl.flux)
    x = np.arange(n, dtype=float)
    lamp = np.interp(x - SHIFT, x, tmpl.flux)
    hdr = fits.Header({"OBJECT": "FeAr", "DATE-OBS": "2023-10-14T12:00:00", "EXPTIME": 5.0})
    fits.PrimaryHDU(lamp.astype(np.float32), header=hdr).writeto(data / "aflamp0001_FeAr.fits")

    # the true pixel -> wavelength relation of the shifted lamp
    ref = wc.recenter_lines(tmpl.flux, tmpl.pixels)
    ok = np.isfinite(ref)
    truth = wc.fit_solution(ref[ok] + SHIFT, tmpl.wavelengths[ok], npix=n)
    wave = truth(x)
    dw = np.abs(np.gradient(wave))

    table = fc.load_standard("feige34", caldir="spec50cal")
    ext = fc.load_extinction("XLT")
    # a smooth stand-in for the star: the published table is not telluric
    # corrected and has few-percent features inside the atmospheric bands
    smooth = np.polynomial.Chebyshev.fit(np.log(table.wave), np.log(table.flam), 3)
    std_flam = np.exp(smooth(np.log(wave)))
    obj_flam = 2e-15 * (1 + 0.2 * np.sin(wave / 500.0))

    def write(name, obj, flam, airmass, exptime, clock):
        rate = flam * 10 ** (0.4 * (_sens_true(wave) - airmass * ext(wave)))
        counts = rate * exptime * dw * np.exp(-airmass * _tau(wave))
        cube = np.array([counts, counts, np.full(n, 50.0), np.sqrt(counts)])[:, None, :]
        h = fits.Header({"OBJECT": obj, "AIRMASS": airmass, "EXPTIME": exptime,
                         "DATE-OBS": f"2023-10-14T{clock}"})
        fits.PrimaryHDU(cube.astype(np.float32), header=h).writeto(data / name)

    write("acrfspec0002_Feige34_bc.fits", "Feige 34", std_flam, 1.10, 60.0, "12:10:00")
    write("acrfspec0003_J1234+5678_bc.fits", "J1234+5678", obj_flam, 1.45, 900.0, "12:30:00")
    return {"root": root, "data": data, "wave": wave, "obj_flam": obj_flam,
            "std_flam": std_flam, "truth": truth}


def test_steps_run_in_order(night):
    root = str(night["root"])
    assert reidentlamp2m.main(["--cwd", root, "--interactive", "never"]) == 0
    assert wavecal2m.main(["--cwd", root]) == 0
    assert telluric_base2m.main(["--cwd", root]) == 0
    data = night["data"]
    assert (data / "database" / "wavesol_aflamp0001_FeAr.json").exists()
    assert (data / "sens.fits").exists()
    assert (night["root"] / "QAplots" / "sensfunc.pdf").exists()
    assert len(list(data.glob("*_tel.fits"))) == 2


def test_wavelength_solution_matches_truth(night):
    sol = wc.WaveSolution.load(night["data"] / "database" / "wavesol_aflamp0001_FeAr.json")
    # XLT Fe-Ar solutions are refitted on the BFOSC atlas lines, so they are
    # judged on those: many lines, small scatter, no wide range without lines.
    assert sol.is_good() and sol.n_used >= 25
    assert sol.rms_pix < 0.5
    assert sol.max_gap < 0.2 * sol.npix
    # The fixture's wavelength scale comes from the template's own record,
    # whose blue identifications differ from the atlas by up to ~5 A (checked
    # against sky and Balmer lines on a real night); redward they agree.
    truth = night["truth"]
    red = truth.used & (truth.wavelengths > 4700)
    at_lines = sol(truth.pixels[red]) - truth.wavelengths[red]
    assert np.max(np.abs(at_lines)) < 2.8
    assert sol.dispersion < 0  # XLT G4 runs red to blue on the detector

    w = read_spec(night["data"] / "wacrfspec0003_J1234+5678_bc.fits")
    assert w.calibrated and w.cdelt > 0 and len(w.bands) == 4
    assert w.header["REFSPEC1"] == "aflamp0001_FeAr"
    # the O2 A band lands at the right wavelength
    sel = (w.wave > 7500) & (w.wave < 7760)
    smooth = np.polyval(np.polyfit(w.wave[sel], w.flux[sel], 2), w.wave[sel])
    assert abs(w.wave[sel][np.argmin(w.flux[sel] / smooth)] - 7630.0) < 6.0


def test_flux_calibration_recovers_target(night):
    c = read_spec(night["data"] / "cwacrfspec0003_J1234+5678_bc.fits")
    truth = np.interp(c.wave, night["wave"][::-1], night["obj_flam"][::-1])
    sel = (c.wave > 4500) & (c.wave < 6700)   # away from edges and telluric bands
    ratio = c.flux[sel] / truth[sel]
    assert abs(np.median(ratio) - 1.0) < 0.04
    assert np.percentile(np.abs(ratio - 1), 95) < 0.08
    assert c.header["BUNIT"] == "erg/cm2/s/A"


def test_telluric_correction_fills_the_bands(night):
    tel = next(p for p in night["data"].glob("J1234+5678_XL*_feige34_tel.fits"))
    before = read_spec(night["data"] / "cwacrfspec0003_J1234+5678_bc.fits")
    after = read_spec(tel)
    truth = np.interp(after.wave, night["wave"][::-1], night["obj_flam"][::-1])
    for lo, hi in ((7590, 7700), (6860, 6960)):
        sel = (after.wave > lo) & (after.wave < hi)
        depth_before = np.max(np.abs(before.flux[sel] / truth[sel] - 1))
        depth_after = np.max(np.abs(after.flux[sel] / truth[sel] - 1))
        assert depth_before > 0.25
        assert depth_after < 0.08
    assert after.header["TELSTD"] == "feige34"


def test_rerun_is_idempotent(night):
    data = night["data"]
    target = data / "cwacrfspec0003_J1234+5678_bc.fits"
    mtime = target.stat().st_mtime_ns
    assert wavecal2m.main(["--cwd", str(night["root"])]) == 0
    assert target.stat().st_mtime_ns == mtime
    # the data directory itself is accepted as --cwd too
    pe = PostExtraction(data)
    assert pe.data_dir == data and pe.telescope == "XLT"


def test_missing_standard_is_not_fatal(tmp_path, night, monkeypatch):
    import shutil

    work = tmp_path / "proj"
    shutil.copytree(night["root"], work)
    for path in (work / "data").glob("*Feige34*"):
        path.unlink()
    for path in list((work / "data").glob("*_tel.fits")) + [work / "data" / "sens.fits"]:
        path.unlink()
    assert wavecal2m.main(["--cwd", str(work), "--overwrite"]) == 0
    assert not (work / "data" / "sens.fits").exists()
    assert (work / "data" / "wacrfspec0003_J1234+5678_bc.fits").exists()


def test_steps_import_without_pyraf():
    assert "pyraf" not in sys.modules
    import pyfosc.cli  # noqa: F401
    import pyfosc.steps.identlamp2m  # noqa: F401


def test_night_with_diluted_standard_and_mislabelled_frame(tmp_path, night):
    """Several targets, a standard whose A band is filled in, and a wrong 'standard'."""
    import csv
    import shutil

    work = tmp_path / "proj"
    shutil.copytree(night["root"], work)
    data = work / "data"
    for path in list(data.glob("w*")) + list(data.glob("c*")) + list(data.glob("*_tel.fits")) \
            + list(data.glob("telluric_*")) + [data / "sens.fits"]:
        path.unlink()
    for path in data.glob("acrf*"):
        path.unlink()
    wave = night["wave"]
    n = len(wave)
    dw = np.abs(np.gradient(wave))
    ext = fc.load_extinction("XLT")
    rng = np.random.default_rng(0)

    def write(name, obj, flam, airmass, exptime, clock, dilute=0.0, noise=0.0):
        trans = np.exp(-airmass * _tau(wave))
        a_band = (wave > 7500) & (wave < 7760)
        trans[a_band] = (trans[a_band] + dilute) / (1 + dilute)   # second-order light
        rate = flam * 10 ** (0.4 * (_sens_true(wave) - airmass * ext(wave)))
        counts = rate * exptime * dw * trans
        if noise:
            counts = counts * (1 + rng.normal(0, noise, n))
        cube = np.array([counts, counts, np.full(n, 50.0), np.sqrt(np.abs(counts))])[:, None, :]
        h = fits.Header({"OBJECT": obj, "AIRMASS": airmass, "EXPTIME": exptime,
                         "DATE-OBS": f"2023-10-14T{clock}"})
        fits.PrimaryHDU(cube.astype(np.float32), header=h).writeto(data / name)

    red = night["std_flam"] * (wave / 5000.0) ** 6          # not a hot standard at all
    write("acrfspec0002_Feige34_bc.fits", "Feige 34", night["std_flam"], 1.0, 60.0,
          "12:10:00", dilute=0.25)
    write("acrfspec0003_G191B2B_bc.fits", "G191-B2B", red, 1.1, 60.0, "12:20:00")
    flat = np.full(n, 2e-15)
    for i, (noise, airmass) in enumerate(((0.01, 1.0), (0.02, 1.2), (0.3, 1.1)), start=4):
        write(f"acrfspec000{i}_T{i}_bc.fits", f"T{i}", flat, airmass, 900.0, f"12:{i}0:00",
              noise=noise)

    pe = PostExtraction(work, overwrite=True)
    pe.apply_wavelength()
    assert [name for _, name, _ in pe.checked_standards()] == ["feige34"]
    pe.build_sensfunc()
    assert read_spec(data / "sens.fits").header["STDNAME"] == "feige34"
    pe.flux_calibrate()
    outputs = pe.telluric_correct()
    assert len(outputs) == 5
    with open(work / "QAplots" / "telluric_summary.csv") as fh:
        rows = {(r["object"], r["band"]): r for r in csv.DictReader(fh)}
    # the bright targets measure the A band themselves; the faint one gets the
    # night's value, which makes up for the diluted standard
    for obj, airmass in (("T4", 1.0), ("T5", 1.2), ("T6", 1.1)):
        a_scale = float(rows[(obj, "O2_7590")]["scale"]) / airmass
        assert 1.25 < a_scale < 1.75, (obj, a_scale)
        assert abs(float(rows[(obj, "O2_6860")]["scale"]) / airmass - 1.0) < 0.12
    assert rows[("T6", "O2_7590")]["origin"] in ("night", "night+fit")
    assert rows[("T4", "O2_7590")]["origin"] == "fit"
    t4 = read_spec(next(data.glob("T4_XL*_feige34_tel.fits")))
    sel = (t4.wave > 7590) & (t4.wave < 7700)
    assert abs(np.median(t4.flux[sel]) / 2e-15 - 1) < 0.05
    assert (work / "QAplots" / "telluric_template_feige34.pdf").exists()
