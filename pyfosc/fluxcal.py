"""Flux calibration without IRAF (``standard`` + ``sensfunc`` + ``calibrate``).

Standard-star tables and extinction curves are the IRAF-format files
shipped under ``iraf_data/onedstds`` and ``extinction``. The conventions
follow IRAF so that results are comparable:

* sensitivity ``S = 2.5 log10[ C / (T * B * F) ] + airmass * ext`` where
  ``C`` are the counts in a bandpass of width ``B`` (Angstrom), ``T`` the
  exposure time and ``F`` the standard's flux in erg/cm2/s/A;
* calibrated flux ``= counts / (T * dlambda) * 10**(-0.4 * (S - airmass * ext))``.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.interpolate import LSQUnivariateSpline

from .paths import find_data_root
from .specio import Spec1D

C_ANGSTROM = 2.99792458e18  # speed of light in Angstrom/s
AB_ZEROPOINT = 48.60

# Preferred order when a star is tabulated in several directories.
CALDIR_PREFERENCE = (
    "spec50cal", "oke1990", "ctionewcal", "subaru", "irscal", "redcal",
    "iidscal", "spechayescal", "spec16cal", "ctiocal", "bstdscal",
)

EXTINCTION_FILES = {
    "XLT": "extinction/baoextinct.dat",
    "LJT": "extinction/LJextinct.dat",
    "HCT": "extinction/LJextinct.dat",
    "P200": "iraf_data/onedstds/kpnoextinct.dat",
}

# Geodetic longitude (deg E), latitude (deg), height (m)
SITES = {
    "XLT": (117.575, 40.3942, 950.0),
    "LJT": (100.031, 26.7092, 3193.0),
    "HCT": (78.9641, 32.7794, 4500.0),
    "P200": (-116.865, 33.3563, 1712.0),
}

_NAME_ALIASES = (("hiltner", "hilt"), ("wolf", "wolf"), ("feige", "feige"))


def _data_root(data_root=None) -> Path:
    return Path(data_root) if data_root else find_data_root(module_path=Path(__file__))


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------

def _read_columns(path, ncol: int) -> np.ndarray:
    rows = []
    for raw in Path(path).read_text().splitlines():
        parts = raw.split()
        if len(parts) < ncol or parts[0].startswith("#"):
            continue
        try:
            rows.append([float(p) for p in parts[:ncol]])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"No numeric rows with {ncol} columns in {path}")
    return np.asarray(rows)


def load_extinction(telescope_or_path, data_root=None) -> Callable[[np.ndarray], np.ndarray]:
    """Return ``ext(wavelength)`` in mag/airmass from an IRAF extinction file."""
    rel = EXTINCTION_FILES.get(str(telescope_or_path))
    path = _data_root(data_root) / rel if rel else Path(telescope_or_path)
    tab = _read_columns(path, 2)
    order = np.argsort(tab[:, 0])
    wave, ext = tab[order, 0], tab[order, 1]

    def extinction(w):
        return np.interp(w, wave, ext)

    return extinction


@dataclass
class StandardTable:
    """Tabulated fluxes of a spectrophotometric standard."""

    name: str
    caldir: str
    wave: np.ndarray   # bandpass centres (Angstrom)
    mag: np.ndarray    # AB magnitudes
    width: np.ndarray  # bandpass widths (Angstrom)
    path: Optional[Path] = None

    @property
    def flam(self) -> np.ndarray:
        """Flux density in erg/cm2/s/A."""
        fnu = 10.0 ** (-0.4 * (self.mag + AB_ZEROPOINT))
        return fnu * C_ANGSTROM / self.wave ** 2


def normalize_star_name(name: str) -> str:
    """Reduce a star name to lowercase alphanumerics (``BD+28 4211`` -> ``bd284211``)."""
    s = str(name).lower().strip()
    for old, new in _NAME_ALIASES:
        s = s.replace(old, new)
    s = re.sub(r"(?<=\d)d(?=\d)", "", s)  # BD+28d4211
    return re.sub(r"[^a-z0-9]", "", s)


def list_standards(data_root=None) -> List[Tuple[str, str, Path]]:
    """All ``(normalized name, caldir, path)`` entries, preferred dirs first."""
    base = _data_root(data_root) / "iraf_data" / "onedstds"
    out = []
    for caldir in CALDIR_PREFERENCE:
        for path in sorted((base / caldir).glob("*.dat")):
            out.append((normalize_star_name(path.stem), caldir, path))
    return out


def find_standard(object_name: str, caldir: Optional[str] = None,
                  data_root=None) -> Optional[Tuple[str, str, Path]]:
    """Match an ``OBJECT`` value to a tabulated standard.

    An exact match of the normalized name wins; otherwise the longest
    table name that the object name starts with (e.g. ``Feige34_std``).
    Returns ``(star, caldir, path)`` or None.
    """
    target = normalize_star_name(object_name)
    if not target:
        return None
    entries = [e for e in list_standards(data_root) if caldir in (None, e[1])]
    for entry in entries:
        if entry[0] == target:
            return entry
    prefixed = [e for e in entries
                if len(e[0]) >= 4 and target.startswith(e[0])
                and not target[len(e[0]):len(e[0]) + 1].isdigit()]
    if prefixed:
        longest = max(len(e[0]) for e in prefixed)
        return next(e for e in prefixed if len(e[0]) == longest)
    return None


def load_standard(star: str, caldir: Optional[str] = None, data_root=None) -> StandardTable:
    """Load a standard-star table by name (and optionally directory)."""
    path = Path(star)
    if path.suffix == ".dat" and path.exists():
        entry = (normalize_star_name(path.stem), path.parent.name, path)
    else:
        entry = find_standard(star, caldir=caldir, data_root=data_root)
        if entry is None:
            where = f" in {caldir}" if caldir else ""
            raise FileNotFoundError(f"Standard star '{star}' not found{where}.")
    tab = _read_columns(entry[2], 3)
    order = np.argsort(tab[:, 0])
    tab = tab[order]
    return StandardTable(name=entry[0], caldir=entry[1], wave=tab[:, 0], mag=tab[:, 1],
                         width=tab[:, 2], path=entry[2])


# ---------------------------------------------------------------------------
# Header helpers
# ---------------------------------------------------------------------------

def get_exptime(header) -> float:
    for key in ("EXPTIME", "EXPOSURE", "ITIME", "EXP_TIME"):
        if key in header:
            try:
                value = float(header[key])
            except (TypeError, ValueError):
                continue
            if value > 0:
                return value
    raise KeyError("No positive exposure time (EXPTIME/EXPOSURE/ITIME) in header.")


def get_obstime(header):
    """Observation start as ``astropy.time.Time`` or None."""
    from astropy.time import Time

    date = str(header.get("DATE-OBS", "")).strip()
    try:
        if date and "T" in date:
            return Time(date, scale="utc")
        if date:
            for key in ("TIME-OBS", "UT", "UTC", "UTSTART"):
                clock = str(header.get(key, "")).strip()
                if ":" in clock:
                    return Time(f"{date}T{clock}", scale="utc")
        for key, fmt in (("MJD-OBS", "mjd"), ("MJD", "mjd"), ("JD", "jd")):
            if key in header:
                return Time(float(header[key]), format=fmt, scale="utc")
        if date:
            return Time(date, scale="utc")
    except (ValueError, TypeError):
        return None
    return None


def compute_airmass(header, telescope: str) -> Optional[float]:
    """Airmass at mid-exposure from coordinates, time and site (or None)."""
    from astropy import units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.utils import iers

    site = SITES.get(telescope)
    time = get_obstime(header)
    if site is None or time is None:
        return None
    coord = None
    for ra_key, dec_key in (("RA", "DEC"), ("CAT-RA", "CAT-DEC"), ("OBJRA", "OBJDEC")):
        if ra_key in header and dec_key in header:
            ra, dec = header[ra_key], header[dec_key]
            try:
                if isinstance(ra, str) and (":" in ra or " " in ra.strip()):
                    coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
                else:
                    coord = SkyCoord(float(ra), float(dec), unit=(u.deg, u.deg))
                break
            except (ValueError, TypeError):
                continue
    if coord is None:
        return None
    try:
        time = time + 0.5 * get_exptime(header) * u.s
    except KeyError:
        pass
    location = EarthLocation.from_geodetic(site[0] * u.deg, site[1] * u.deg, site[2] * u.m)
    with iers.conf.set_temp("auto_download", False), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        altaz = coord.transform_to(AltAz(obstime=time, location=location))
    alt = altaz.alt.deg
    if alt <= 3:
        return None
    return float(1.0 / np.sin(np.radians(alt + 244.0 / (165.0 + 47.0 * alt ** 1.1))))


def get_airmass(header, telescope: Optional[str] = None, default: float = 1.0) -> float:
    """Airmass from the header, else computed, else ``default`` (with a warning)."""
    for key in ("AIRMASS", "SECZ", "AMSTART"):
        if key in header:
            try:
                value = float(header[key])
            except (TypeError, ValueError):
                continue
            if 0.99 <= value < 20:
                return value
    if telescope:
        try:
            value = compute_airmass(header, telescope)
        except Exception:  # coordinates in an unexpected format
            value = None
        if value is not None:
            return value
    warnings.warn(f"No usable airmass for {header.get('OBJECT', 'spectrum')}; "
                  f"assuming {default}.")
    return default


# ---------------------------------------------------------------------------
# Sensitivity function
# ---------------------------------------------------------------------------

def band_counts(spec: Spec1D, centres, widths) -> np.ndarray:
    """Counts of band 0 inside each bandpass (NaN when not fully covered)."""
    n = spec.npix
    edges = spec.crval + spec.cdelt * (np.arange(n + 1) - 0.5)
    cum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(spec.flux))])
    if edges[0] > edges[-1]:
        edges, cum = edges[::-1], (cum[-1] - cum)[::-1]
    lo = np.asarray(centres) - 0.5 * np.asarray(widths)
    hi = np.asarray(centres) + 0.5 * np.asarray(widths)
    counts = np.interp(hi, edges, cum) - np.interp(lo, edges, cum)
    inside = (lo >= edges[0]) & (hi <= edges[-1])
    return np.where(inside, counts, np.nan)


def _in_bands(wave, bands) -> np.ndarray:
    wave = np.asarray(wave)
    mask = np.zeros(wave.shape, dtype=bool)
    for lo, hi in bands:
        mask |= (wave >= lo) & (wave <= hi)
    return mask


@dataclass
class SensFunc:
    """Fitted sensitivity curve with the points it came from."""

    wave: np.ndarray      # evaluation grid (Angstrom)
    sens: np.ndarray      # 2.5 log10(counts/s/A per erg/cm2/s/A) at zero airmass
    points_wave: np.ndarray
    points_sens: np.ndarray
    points_used: np.ndarray
    rms: float
    star: str = ""
    caldir: str = ""
    airmass: float = 1.0

    def __call__(self, wave):
        return np.interp(wave, self.wave, self.sens)


def _fit_smooth(x, y, order: int, good: np.ndarray):
    """Cubic spline with ``order`` pieces, falling back to a polynomial."""
    xs, ys = x[good], y[good]
    pieces = int(max(1, min(order, len(xs) // 4)))
    if pieces >= 2:
        knots = np.linspace(xs[0], xs[-1], pieces + 1)[1:-1]
        try:
            spline = LSQUnivariateSpline(xs, ys, knots, k=3, ext=3)
            return lambda w: spline(np.clip(w, xs[0], xs[-1]))
        except ValueError:
            pass  # knots without data in between (large masked gaps)
    deg = int(max(1, min(order, 5, len(xs) - 2)))
    poly = np.polynomial.Chebyshev.fit(xs, ys, deg)
    return lambda w: poly(np.clip(w, xs[0], xs[-1]))


def sensfunc(std_spec: Spec1D, table: StandardTable, extinction, airmass: float,
             exptime: float, order: int = 6, sigma: float = 3.0, niter: int = 5,
             mask_bands: Optional[Sequence[Tuple[float, float]]] = None) -> SensFunc:
    """Derive the sensitivity function from a wavelength-calibrated standard.

    ``mask_bands`` (default: the telluric bands) are excluded from the fit so
    that the curve interpolates smoothly across atmospheric absorption.
    """
    if not std_spec.calibrated:
        raise ValueError("The standard-star spectrum must be wavelength calibrated.")
    if mask_bands is None:
        from .telluric import TELLURIC_BANDS
        mask_bands = [(lo, hi) for _, lo, hi in TELLURIC_BANDS]
    counts = band_counts(std_spec, table.wave, table.width)
    rate = counts / table.width / exptime
    flam = table.flam
    valid = np.isfinite(rate) & (rate > 0) & (flam > 0)
    if valid.sum() < 4:
        raise ValueError(
            f"Only {int(valid.sum())} standard-star bandpasses overlap the spectrum.")
    wave = table.wave[valid]
    sens = 2.5 * np.log10(rate[valid] / flam[valid]) + airmass * extinction(wave)
    good = ~_in_bands(wave, mask_bands)
    if good.sum() < 4:
        good = np.ones(len(wave), dtype=bool)
    usable = good.copy()
    model = _fit_smooth(wave, sens, order, good)
    for _ in range(niter):
        res = sens - model(wave)
        scatter = 1.4826 * np.median(np.abs(res[good] - np.median(res[good])))
        new_good = usable & (np.abs(res) <= sigma * max(scatter, 1e-4))
        if new_good.sum() < 4 or np.array_equal(new_good, good):
            break
        good = new_good
        model = _fit_smooth(wave, sens, order, good)
    rms = float(np.sqrt(np.mean((sens[good] - model(wave[good])) ** 2)))
    grid = std_spec.wave
    return SensFunc(wave=grid, sens=model(grid), points_wave=wave, points_sens=sens,
                    points_used=good, rms=rms, star=table.name, caldir=table.caldir,
                    airmass=airmass)


def calibrate(spec: Spec1D, sens_wave, sens_values, extinction, airmass: float,
              exptime: float) -> Spec1D:
    """Convert counts to erg/cm2/s/A, correcting for atmospheric extinction."""
    if not spec.calibrated:
        raise ValueError("The spectrum must be wavelength calibrated first.")
    wave = spec.wave
    sens = np.interp(wave, sens_wave, sens_values)
    factor = 10.0 ** (-0.4 * (sens - airmass * extinction(wave))) / (exptime * abs(spec.cdelt))
    out = spec.copy()
    out.bands = [band * factor for band in spec.bands]
    hdr = out.header
    hdr["BUNIT"] = "erg/cm2/s/A"
    hdr["CA-FLAG"] = (0, "flux calibrated")
    hdr["EX-FLAG"] = (0, "extinction corrected")
    hdr["FCAIRM"] = (round(float(airmass), 4), "airmass used for extinction correction")
    return out


def plot_sensfunc(sf: SensFunc, path=None):
    """Sensitivity points and fitted curve."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(9, 5))
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    used = sf.points_used
    ax.plot(sf.points_wave[used], sf.points_sens[used], ".", ms=3, label="used")
    ax.plot(sf.points_wave[~used], sf.points_sens[~used], "x", ms=4, color="C3",
            label="masked / rejected")
    ax.plot(sf.wave, sf.sens, "k-", lw=1, label="fit")
    ax.set_xlabel("Wavelength (A)")
    ax.set_ylabel("Sensitivity (mag)")
    ax.set_title(f"{sf.star} ({sf.caldir}), airmass {sf.airmass:.2f}, rms {sf.rms:.3f} mag",
                 fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()
    if path is not None:
        fig.savefig(path)
    return fig
