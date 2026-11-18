"""Read and write the IRAF-style 1D spectrum files used across pyfosc.

Two layouts are handled:

* plain 1D images (lamp spectra, sensitivity functions), and
* "multispec"-like cubes of shape ``(nband, 1, npix)`` written by
  :class:`pyfosc.extraction.Extract1dSpec` (optimal, raw, background, sigma).

In both cases the first axis carries a linear WCS (``CRVAL1``, ``CRPIX1``,
``CD1_1``/``CDELT1``) that is either in pixels or, after dispersion
correction, in Angstroms.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Sequence

import numpy as np
from astropy.io import fits

BAND_NAMES = ("flux", "raw", "background", "sigma")

_BANDID = (
    "spectrum - background fit, weights variance, clean no",
    "raw - background fit, weights none, clean no",
    "background - background fit",
    "sigma - background fit, weights variance, clean no",
)

# Keywords that describe the old WCS and must not leak into rewritten files.
_WCS_KEYS = (
    "CDELT1", "CD1_1", "CRVAL1", "CRPIX1", "CTYPE1", "CUNIT1", "DC-FLAG",
    "WAT0_001", "WAT1_001", "WAT2_001", "WAT3_001", "WCSDIM",
    "CD2_2", "CD3_3", "LTM1_1", "LTM2_2", "LTM3_3", "CTYPE2", "CTYPE3",
    "DISPAXIS",
)


@dataclass
class Spec1D:
    """A 1D spectrum with one or more bands on a common linear axis."""

    bands: List[np.ndarray]
    header: fits.Header = field(default_factory=fits.Header)
    crval: float = 1.0
    cdelt: float = 1.0
    calibrated: bool = False  # True when the axis is wavelength (Angstrom)
    path: Optional[Path] = None

    @property
    def npix(self) -> int:
        return len(self.bands[0])

    @property
    def wave(self) -> np.ndarray:
        """Coordinate of each pixel centre (pixels or Angstroms)."""
        return self.crval + self.cdelt * np.arange(self.npix)

    @property
    def flux(self) -> np.ndarray:
        return self.bands[0]

    @property
    def sigma(self) -> Optional[np.ndarray]:
        return self.bands[3] if len(self.bands) >= 4 else None

    def copy(self) -> "Spec1D":
        return Spec1D(
            bands=[b.copy() for b in self.bands],
            header=self.header.copy(),
            crval=self.crval,
            cdelt=self.cdelt,
            calibrated=self.calibrated,
            path=self.path,
        )


def read_spec(path) -> Spec1D:
    """Read a 1D or ``(nband, 1, npix)`` spectrum file."""
    path = Path(path)
    with fits.open(path) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header.copy()
    if data.ndim == 1:
        bands = [data]
    elif data.ndim == 2:
        bands = [row for row in data]
    elif data.ndim == 3:
        bands = [plane[0] for plane in data]
    else:
        raise ValueError(f"{path}: unsupported data shape {data.shape}")
    cdelt = float(header.get("CD1_1", header.get("CDELT1", 1.0)))
    crpix = float(header.get("CRPIX1", 1.0))
    crval = float(header.get("CRVAL1", 1.0)) + (1.0 - crpix) * cdelt
    calibrated = "DC-FLAG" in header
    return Spec1D(bands=bands, header=header, crval=crval, cdelt=cdelt,
                  calibrated=calibrated, path=path)


def write_spec(path, bands: Sequence[np.ndarray], header=None, crval: float = 1.0,
               cdelt: float = 1.0, calibrated: bool = False,
               overwrite: bool = True) -> Path:
    """Write bands to an IRAF-readable file with a linear axis.

    One band gives a plain 1D image; several bands give the
    ``(nband, 1, npix)`` layout produced by the extraction step.
    """
    path = Path(path)
    bands = [np.asarray(b, dtype=np.float32) for b in bands]
    hdr = fits.Header() if header is None else header.copy()
    for key in _WCS_KEYS:
        hdr.remove(key, ignore_missing=True)
    for key in [k for k in hdr if k.startswith(("BANDID", "APNUM"))]:
        if len(bands) == 1 or key.startswith("BANDID"):
            hdr.remove(key, ignore_missing=True)
    if len(bands) == 1:
        data = bands[0]
    else:
        data = np.expand_dims(np.array(bands), 1)
    hdu = fits.PrimaryHDU(data=data, header=hdr)
    h = hdu.header
    h["CRVAL1"] = float(crval)
    h["CRPIX1"] = 1.0
    h["CD1_1"] = float(cdelt)
    h["CDELT1"] = float(cdelt)
    if calibrated:
        h["CTYPE1"] = "LINEAR"
        h["DC-FLAG"] = 0
        h["WAT1_001"] = "wtype=linear label=Wavelength units=angstroms"
    else:
        h["CTYPE1"] = "PIXEL"
        h["WAT1_001"] = "wtype=linear label=Pixel"
    if len(bands) > 1:
        h["WCSDIM"] = 3
        h["WAT0_001"] = "system=equispec"
        h["WAT2_001"] = "wtype=linear"
        h["WAT3_001"] = "wtype=linear"
        h["CTYPE2"] = "LINEAR"
        h["CTYPE3"] = "LINEAR"
        h["CD2_2"] = 1.0
        h["CD3_3"] = 1.0
        h["LTM1_1"] = 1.0
        h["LTM2_2"] = 1.0
        h["LTM3_3"] = 1.0
        for i, text in enumerate(_BANDID[: len(bands)], start=1):
            h[f"BANDID{i}"] = text
        h.setdefault("APNUM1", "1 1")
    hdu.writeto(path, overwrite=overwrite)
    return path
