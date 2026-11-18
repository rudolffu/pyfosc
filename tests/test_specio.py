import numpy as np
from astropy.io import fits

from pyfosc.specio import read_spec, write_spec


def test_roundtrip_1d_pixel(tmp_path):
    flux = np.arange(50, dtype=float)
    path = write_spec(tmp_path / "lamp.fits", [flux])
    spec = read_spec(path)
    assert len(spec.bands) == 1 and spec.npix == 50
    assert not spec.calibrated
    np.testing.assert_allclose(spec.flux, flux)
    np.testing.assert_allclose(spec.wave, np.arange(50) + 1.0)


def test_roundtrip_multiband_linear(tmp_path):
    rng = np.random.default_rng(0)
    bands = [rng.normal(size=40) for _ in range(4)]
    hdr = fits.Header({"OBJECT": "Feige 34", "EXPTIME": 300.0, "CTYPE1": "PIXEL"})
    path = write_spec(tmp_path / "w.fits", bands, header=hdr, crval=4000.0, cdelt=2.5,
                      calibrated=True)
    with fits.open(path) as hdul:
        assert hdul[0].data.shape == (4, 1, 40)
        assert hdul[0].header["CTYPE1"] == "LINEAR"
        assert hdul[0].header["DC-FLAG"] == 0
    spec = read_spec(path)
    assert spec.calibrated and len(spec.bands) == 4
    assert spec.header["OBJECT"] == "Feige 34"
    np.testing.assert_allclose(spec.wave, 4000.0 + 2.5 * np.arange(40))
    np.testing.assert_allclose(spec.sigma, bands[3], rtol=1e-6)


def test_reads_extraction_layout(tmp_path):
    # layout written by Extract1dSpec: (4, 1, N) with a pixel WCS
    data = np.ones((4, 1, 30), dtype=np.float32)
    hdu = fits.PrimaryHDU(data)
    hdu.header["CRVAL1"] = 1
    hdu.header["CRPIX1"] = 1
    hdu.header["CD1_1"] = 1
    hdu.writeto(tmp_path / "a.fits")
    spec = read_spec(tmp_path / "a.fits")
    assert spec.npix == 30 and len(spec.bands) == 4 and not spec.calibrated
