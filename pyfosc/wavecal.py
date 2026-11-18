"""Wavelength calibration without IRAF.

The automatic path mirrors IRAF ``reidentify``: a new arc-lamp spectrum is
aligned with a packaged template (``database/lamp_*.fits``) by
cross-correlation, the features stored in the matching ``identify`` record
(``database/idlamp_*``) are re-centred on the new spectrum, and a
sigma-clipped Legendre polynomial is fitted to pixel -> wavelength.
:func:`dispcor` then rebins spectra onto a linear wavelength grid.

Pixel coordinates in this module are 0-based array indices; IRAF records
are 1-based and converted on reading.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from astropy.io import fits
from scipy.ndimage import median_filter
from scipy.signal import correlate, find_peaks

from .paths import find_data_root
from .specio import Spec1D

__all__ = [
    "Template", "WaveSolution", "WaveCalibrator", "read_identify_records",
    "read_identify_feature_list", "lamp_kind", "list_templates", "load_linelist",
    "gaussian", "subpix_peaks", "recenter_lines", "match_template", "fit_solution",
    "solve_from_template", "auto_identify", "add_lines_from_list", "refine_with_linelist",
    "rebin_flux",
    "dispcor", "plot_solution", "identify_interactive",
]

TELESCOPE_PREFIX = {"XLT": "XL", "LJT": "LJ", "HCT": "HCT", "P200": "P200"}


# ---------------------------------------------------------------------------
# IRAF identify records and packaged templates
# ---------------------------------------------------------------------------

def read_identify_records(filename) -> List[Dict]:
    """Parse every record of an IRAF ``identify`` database file.

    Returns a list of dicts with keys ``image``, ``pixels`` (0-based),
    ``wavelengths`` (reference values) and ``names``.
    """
    records: List[Dict] = []
    current: Optional[Dict] = None
    remaining = 0
    with open(filename) as fh:
        for raw in fh:
            parts = raw.split()
            if not parts:
                continue
            if parts[0] == "begin":
                current = {"image": "", "pixels": [], "wavelengths": [], "names": []}
                records.append(current)
                remaining = 0
            elif current is None:
                continue
            elif parts[0] == "id" and len(parts) > 1:
                current["image"] = parts[1]
            elif parts[0] == "features" and len(parts) > 1:
                remaining = int(parts[1])
            elif remaining > 0:
                remaining -= 1
                try:
                    pix, ref = float(parts[0]), float(parts[2])
                except (ValueError, IndexError):
                    continue
                if ref <= 0:  # INDEF or unset reference wavelength
                    continue
                current["pixels"].append(pix - 1.0)
                current["wavelengths"].append(ref)
                current["names"].append(parts[6] if len(parts) > 6 else "")
    for rec in records:
        for key in ("pixels", "wavelengths"):
            rec[key] = np.asarray(rec[key], dtype=float)
    return [rec for rec in records if len(rec["pixels"]) > 0]


def read_identify_feature_list(filename):
    """Return the last feature list of an IRAF identify file as a DataFrame.

    Kept for backward compatibility; positions are IRAF 1-based pixels.
    """
    import pandas as pd

    records = read_identify_records(filename)
    if not records:
        raise ValueError("Feature list not found in the file.")
    rec = records[-1]
    return pd.DataFrame({
        "Position": rec["pixels"] + 1.0,
        "ObservedWavelength": rec["wavelengths"],
        "ReferenceWavelength": rec["wavelengths"],
        "LineName": rec["names"],
    })


@dataclass
class Template:
    """A packaged lamp spectrum with its identified features."""

    name: str
    flux: np.ndarray
    pixels: np.ndarray
    wavelengths: np.ndarray
    names: List[str]
    kind: str


def lamp_kind(header, telescope: str) -> str:
    """Classify a lamp spectrum as ``'fear'`` or ``'hene'``."""
    text = " ".join(str(header.get(k, "")) for k in ("FILTER", "OBJECT", "LAMP")).lower()
    if "neon" in text or "helium" in text or "hene" in text or "he+ne" in text:
        return "hene"
    if "fe_argon" in text or "fear" in text or "fe-ar" in text:
        return "fear"
    return "hene" if telescope == "LJT" else "fear"


def _template_kind(name: str) -> str:
    return "hene" if "he+ne" in name.lower() else "fear"


def list_templates(telescope: str, kind: Optional[str] = None,
                   grism: Optional[str] = None, data_root=None) -> List[Template]:
    """Load packaged templates for a telescope, best name matches first.

    Every record of an ``identify`` file becomes one template, so files
    holding several feature lists contribute several candidates.
    """
    root = Path(data_root) if data_root else find_data_root(module_path=Path(__file__))
    prefix = TELESCOPE_PREFIX.get(telescope, telescope)
    templates: List[Tuple[int, Template]] = []
    for fpath in sorted((root / "database").glob(f"lamp_{prefix}_*.fits")):
        idfile = fpath.with_name("id" + fpath.stem)
        if not idfile.exists():
            continue
        tkind = _template_kind(fpath.stem)
        if kind and tkind != kind:
            continue
        flux = np.asarray(fits.getdata(fpath), dtype=float).squeeze()
        tokens = [t.lower() for t in re.split(r"[_+]", fpath.stem)]
        rank = 0
        if grism:
            g = grism.lower()
            rank = 0 if any(t in (g, "new" + g) for t in tokens) else 1
        for i, rec in enumerate(read_identify_records(idfile)):
            name = fpath.stem if i == 0 else f"{fpath.stem}#{i + 1}"
            templates.append((rank, Template(
                name=name, flux=flux, pixels=rec["pixels"],
                wavelengths=rec["wavelengths"], names=rec["names"], kind=tkind)))
    templates.sort(key=lambda item: item[0])
    return [t for _, t in templates]


def load_linelist(kind: str, data_root=None) -> Tuple[np.ndarray, List[str]]:
    """Return ``(wavelengths, names)`` of known lines for a lamp kind.

    He-Ne(-Ar) lines come from ``iraf_data/newlinelists/henear_scidoc2215.txt``;
    Fe-Ar lines are the union of the reference wavelengths in the packaged
    identify records (originally taken from IRAF ``fear.dat``).
    """
    root = Path(data_root) if data_root else find_data_root(module_path=Path(__file__))
    lines: List[Tuple[float, str]] = []
    if kind == "hene":
        path = root / "iraf_data" / "newlinelists" / "henear_scidoc2215.txt"
        if path.exists():
            for raw in path.read_text().splitlines():
                parts = raw.split()
                if not parts or parts[0].startswith("#"):
                    continue
                try:
                    lines.append((float(parts[0]), parts[1] if len(parts) > 1 else ""))
                except ValueError:
                    continue
    for idfile in sorted((root / "database").glob("idlamp_*")):
        if _template_kind(idfile.name) != kind:
            continue
        for rec in read_identify_records(idfile):
            lines.extend(zip(rec["wavelengths"].tolist(), rec["names"]))
    lines.sort(key=lambda item: item[0])
    waves: List[float] = []
    names: List[str] = []
    for wave, name in lines:
        if waves and abs(wave - waves[-1]) < 0.3:
            continue
        waves.append(wave)
        names.append(name)
    return np.asarray(waves), names


# ---------------------------------------------------------------------------
# Line centring and template matching
# ---------------------------------------------------------------------------

def gaussian(x, amp, mean, sigma):
    return amp * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))


def subpix_peaks(lamp_spectrum_flux, height=4000):
    """Gaussian-refined positions of all peaks above ``height``."""
    flux = np.asarray(lamp_spectrum_flux, dtype=float)
    peaks, _ = find_peaks(flux, height=height)
    centres = recenter_lines(flux, peaks.astype(float), fwidth=10)
    return [float(c) for c in centres if np.isfinite(c)]


def _noise_level(flux: np.ndarray) -> float:
    diff = np.diff(flux)
    return 1.4826 * np.median(np.abs(diff - np.median(diff))) / np.sqrt(2.0)


def recenter_lines(flux, guesses, fwidth: float = 10.0, min_snr: float = 5.0,
                   search: Optional[float] = None) -> np.ndarray:
    """Re-centre lines near ``guesses`` (0-based pixels).

    Lamp lines taken through a wide slit are flat-topped, so instead of a
    Gaussian fit the centre is the centroid of the flux above half maximum
    (similar in spirit to IRAF ``center1d``). The peak is searched within
    ``search`` pixels of each guess (default ``fwidth/2``). Entries are NaN
    where no significant, isolated line is found within that distance.
    """
    flux = np.nan_to_num(np.asarray(flux, dtype=float))
    n = len(flux)
    if search is None:
        search = fwidth / 2.0
    half = max(1, int(round(search)))
    reach = max(4, int(round(1.5 * fwidth)))
    noise = _noise_level(flux)
    out = np.full(len(guesses), np.nan)
    for i, guess in enumerate(np.asarray(guesses, dtype=float)):
        if not np.isfinite(guess) or guess < 1 or guess > n - 2:
            continue
        g = int(round(guess))
        lo, hi = max(0, g - half), min(n, g + half + 1)
        peak = lo + int(np.argmax(flux[lo:hi]))
        wlo, whi = max(0, peak - reach), min(n, peak + reach + 1)
        bg = float(np.min(flux[wlo:whi]))
        amp = float(flux[peak] - bg)
        if amp <= min_snr * noise or amp <= 0:
            continue
        level = bg + 0.5 * amp
        left = peak
        while left > wlo and flux[left - 1] > level:
            left -= 1
        right = peak
        while right < whi - 1 and flux[right + 1] > level:
            right += 1
        if left == wlo or right == whi - 1:
            continue  # no half-maximum crossing: blend, edge or broad feature
        x = np.arange(left, right + 1, dtype=float)
        w = flux[left:right + 1] - level
        if w.sum() <= 0:
            continue
        centre = float(np.sum(w * x) / np.sum(w))
        if abs(centre - guess) > search:
            continue
        out[i] = centre
    return out


def _prep_for_xcorr(flux: np.ndarray) -> np.ndarray:
    """Baseline-subtract and compress a lamp spectrum for correlation.

    The square root keeps a handful of very bright lines from dominating.
    """
    flux = np.nan_to_num(np.asarray(flux, dtype=float))
    size = min(len(flux) // 2 * 2 - 1, 151)
    line = np.clip(flux - median_filter(flux, size=max(size, 3), mode="nearest"), 0, None)
    top = np.percentile(line, 99.5)
    if top <= 0:
        top = line.max() or 1.0
    return np.sqrt(np.clip(line / top, 0, 1))


def match_template(flux, tmpl_flux, stretches: Optional[Sequence[float]] = None,
                   max_shift: Optional[int] = None) -> Dict:
    """Find the linear map from template pixels to lamp pixels.

    Returns a dict with ``shift`` and ``stretch`` such that
    ``lamp_pixel = shift + stretch * template_pixel``, the normalised
    correlation ``score`` and ``flipped`` (True when the lamp runs in the
    opposite direction; the map then applies to ``flux[::-1]``).
    """
    lamp = _prep_for_xcorr(flux)
    tmpl = _prep_for_xcorr(tmpl_flux)
    n = len(lamp)
    if stretches is None:
        stretches = np.linspace(0.96, 1.04, 17)
    if max_shift is None:
        max_shift = int(0.4 * max(n, len(tmpl)))
    best = {"score": -np.inf}
    grid = np.arange(n, dtype=float)
    for flipped in (False, True):
        target = lamp[::-1] if flipped else lamp
        norm_t = np.sqrt(np.sum(target ** 2))
        for b in stretches:
            warped = np.interp(grid / b, np.arange(len(tmpl)), tmpl, left=0.0, right=0.0)
            norm_w = np.sqrt(np.sum(warped ** 2))
            if norm_t == 0 or norm_w == 0:
                continue
            cc = correlate(target, warped, mode="full", method="fft")
            lags = np.arange(-(n - 1), n)
            keep = np.abs(lags) <= max_shift
            j = int(np.argmax(np.where(keep, cc, -np.inf)))
            score = cc[j] / (norm_t * norm_w)
            if score > best["score"]:
                best = {"score": float(score), "shift": float(lags[j]),
                        "stretch": float(b), "flipped": flipped}
    return best


# ---------------------------------------------------------------------------
# Wavelength solution
# ---------------------------------------------------------------------------

@dataclass
class WaveSolution:
    """Pixel -> wavelength polynomial with the lines it was fitted to."""

    coeffs: np.ndarray
    npix: int
    pixels: np.ndarray
    wavelengths: np.ndarray
    used: np.ndarray
    names: List[str] = field(default_factory=list)
    rms: float = np.nan
    template: str = ""
    lamp: str = ""
    kind: str = ""
    method: str = "auto"
    score: float = np.nan  # template cross-correlation (1 = identical lamp)

    def __post_init__(self):
        self.coeffs = np.asarray(self.coeffs, dtype=float)
        self.pixels = np.asarray(self.pixels, dtype=float)
        self.wavelengths = np.asarray(self.wavelengths, dtype=float)
        self.used = np.asarray(self.used, dtype=bool)
        self._poly = np.polynomial.Legendre(self.coeffs, domain=[0, self.npix - 1])

    def __call__(self, pixel):
        return self._poly(np.asarray(pixel, dtype=float))

    @property
    def degree(self) -> int:
        return len(self.coeffs) - 1

    @property
    def dispersion(self) -> float:
        """Mean dispersion in Angstrom per pixel (signed)."""
        return float((self(self.npix - 1) - self(0)) / (self.npix - 1))

    @property
    def n_used(self) -> int:
        return int(self.used.sum())

    @property
    def rms_pix(self) -> float:
        return float(self.rms / abs(self.dispersion)) if self.dispersion else np.inf

    @property
    def residuals(self) -> np.ndarray:
        return self.wavelengths - self(self.pixels)

    def is_monotonic(self) -> bool:
        d = np.diff(self(np.arange(self.npix)))
        return bool(np.all(d > 0) or np.all(d < 0))

    def is_good(self, max_rms_pix: float = 0.7, min_lines: int = 6) -> bool:
        return (self.n_used >= max(min_lines, self.degree + 2)
                and self.rms_pix < max_rms_pix and self.is_monotonic())

    def pixel_of(self, wavelength) -> np.ndarray:
        """Invert the solution (valid inside the detector range)."""
        grid = np.arange(self.npix, dtype=float)
        waves = self(grid)
        if waves[0] > waves[-1]:
            grid, waves = grid[::-1], waves[::-1]
        return np.interp(wavelength, waves, grid, left=np.nan, right=np.nan)

    def summary(self) -> str:
        return (f"{self.lamp or 'lamp'}: {self.n_used}/{len(self.pixels)} lines, "
                f"degree {self.degree}, rms {self.rms:.3f} A ({self.rms_pix:.2f} pix), "
                f"dispersion {self.dispersion:.3f} A/pix, "
                f"range {self(0):.0f}-{self(self.npix - 1):.0f} A"
                + (f", template {self.template}" if self.template else "")
                + (f" (match {self.score:.2f})" if np.isfinite(self.score) else ""))

    def save(self, path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "model": "Legendre", "domain": [0, self.npix - 1],
            "coeffs": self.coeffs.tolist(), "npix": self.npix,
            "pixel_origin": 0, "rms": self.rms, "template": self.template,
            "lamp": self.lamp, "kind": self.kind, "method": self.method,
            "score": None if not np.isfinite(self.score) else float(self.score),
            "lines": [
                {"pixel": float(p), "wavelength": float(w), "name": n, "used": bool(u)}
                for p, w, n, u in zip(self.pixels, self.wavelengths, self._names(), self.used)
            ],
        }
        path.write_text(json.dumps(payload, indent=1))
        return path

    def _names(self) -> List[str]:
        names = list(self.names)
        return names + [""] * (len(self.pixels) - len(names))

    @classmethod
    def load(cls, path) -> "WaveSolution":
        d = json.loads(Path(path).read_text())
        lines = d.get("lines", [])
        return cls(
            coeffs=d["coeffs"], npix=d["npix"],
            pixels=[ln["pixel"] for ln in lines],
            wavelengths=[ln["wavelength"] for ln in lines],
            used=[ln["used"] for ln in lines],
            names=[ln.get("name", "") for ln in lines],
            rms=d.get("rms", np.nan), template=d.get("template", ""),
            lamp=d.get("lamp", ""), kind=d.get("kind", ""),
            method=d.get("method", "auto"),
            score=np.nan if d.get("score") is None else d["score"])


def fit_solution(pixels, wavelengths, npix: int, degree: int = 4, sigma: float = 3.0,
                 niter: int = 5, names: Optional[Sequence[str]] = None,
                 **meta) -> WaveSolution:
    """Sigma-clipped Legendre fit of wavelength against pixel."""
    x = np.asarray(pixels, dtype=float)
    y = np.asarray(wavelengths, dtype=float)
    if len(x) < 2:
        raise ValueError("At least two lines are needed for a wavelength solution.")
    order = np.argsort(x)
    x, y = x[order], y[order]
    names = [list(names)[i] for i in order] if names is not None else [""] * len(x)
    used = np.ones(len(x), dtype=bool)
    poly = None
    for _ in range(max(1, niter)):
        deg = int(min(degree, max(1, (used.sum() - 1) // 2)))
        poly = np.polynomial.Legendre.fit(x[used], y[used], deg, domain=[0, npix - 1])
        res = y - poly(x)
        scatter = 1.4826 * np.median(np.abs(res[used] - np.median(res[used])))
        scatter = max(scatter, 1e-6)
        new_used = np.abs(res) <= sigma * max(scatter, np.std(res[used]) / 2.0)
        if new_used.sum() < deg + 2 or np.array_equal(new_used, used):
            break
        used = new_used
    deg = int(min(degree, max(1, (used.sum() - 1) // 2)))
    poly = np.polynomial.Legendre.fit(x[used], y[used], deg, domain=[0, npix - 1])
    rms = float(np.sqrt(np.mean((y[used] - poly(x[used])) ** 2)))
    return WaveSolution(coeffs=poly.coef, npix=npix, pixels=x, wavelengths=y,
                        used=used, names=names, rms=rms, **meta)


def solve_from_template(flux, template: Template, fwidth: float = 10.0,
                        degree: int = 4, min_score: float = 0.2) -> Optional[WaveSolution]:
    """Transfer a template's features to ``flux`` and fit a solution."""
    flux = np.nan_to_num(np.asarray(flux, dtype=float))
    n = len(flux)
    match = match_template(flux, template.flux)
    if not np.isfinite(match["score"]) or match["score"] < min_score:
        return None
    work = flux[::-1] if match["flipped"] else flux
    guesses = match["shift"] + match["stretch"] * template.pixels
    centres = recenter_lines(work, guesses, fwidth=fwidth)
    ok = np.isfinite(centres)
    if ok.sum() >= 4:
        # refine the template->lamp map with the lines found so far
        coef = np.polyfit(template.pixels[ok], centres[ok], 2 if ok.sum() >= 6 else 1)
        res = centres[ok] - np.polyval(coef, template.pixels[ok])
        keep = np.abs(res) <= max(3 * 1.4826 * np.median(np.abs(res)), 0.5)
        if keep.sum() >= 4:
            coef = np.polyfit(template.pixels[ok][keep], centres[ok][keep],
                              2 if keep.sum() >= 6 else 1)
        guesses = np.polyval(coef, template.pixels)
        centres = recenter_lines(work, guesses, fwidth=fwidth, search=max(2.0, fwidth / 4.0))
        ok = np.isfinite(centres)
    if ok.sum() < 3:
        return None
    pix = centres[ok]
    if match["flipped"]:
        pix = (n - 1) - pix
    names = [template.names[i] for i in np.flatnonzero(ok)]
    return fit_solution(pix, template.wavelengths[ok], npix=n, degree=degree,
                        names=names, template=template.name, kind=template.kind,
                        method="auto", score=match["score"])


def auto_identify(flux, templates: Sequence[Template], fwidth: float = 10.0,
                  degree: int = 4) -> Optional[WaveSolution]:
    """Try every template and return the best solution (or None).

    Good solutions from the template that resembles the lamp most (highest
    cross-correlation, in steps of 0.05) win; the number of lines and the
    scatter break ties. A template of the same grism and slit shares the
    lamp's line blends, so its identifications transfer most reliably.
    """
    best: Optional[WaveSolution] = None
    for tmpl in templates:
        sol = solve_from_template(flux, tmpl, fwidth=fwidth, degree=degree)
        if sol is None:
            continue
        if best is None:
            best = sol
            continue
        if _rank(sol) > _rank(best):
            best = sol
    return best


def _rank(sol: "WaveSolution"):
    score = sol.score if np.isfinite(sol.score) else 0.0
    return (sol.is_good(), round(score / 0.05), sol.n_used, -sol.rms_pix)


def add_lines_from_list(flux, solution: WaveSolution, list_waves, list_names=None,
                        fwidth: float = 10.0, tol_pix: float = 1.5):
    """Find extra catalogue lines predicted by ``solution``.

    Returns ``(pixels, wavelengths, names)`` of the current lines plus the
    newly matched ones. Catalogue lines closer together than ``fwidth/2``
    pixels are skipped as ambiguous.
    """
    list_waves = np.asarray(list_waves, dtype=float)
    list_names = list(list_names) if list_names is not None else [""] * len(list_waves)
    pred = solution.pixel_of(list_waves)
    pix = list(solution.pixels[solution.used])
    wav = list(solution.wavelengths[solution.used])
    nam = [n for n, u in zip(solution._names(), solution.used) if u]
    order = np.argsort(pred)
    pred_sorted = pred[order]
    for rank, j in enumerate(order):
        p = pred[j]
        if not np.isfinite(p):
            continue
        neighbours = [pred_sorted[k] for k in (rank - 1, rank + 1)
                      if 0 <= k < len(order) and np.isfinite(pred_sorted[k])]
        if any(abs(p - q) < fwidth / 2.0 for q in neighbours):
            continue
        if any(abs(p - q) < fwidth / 2.0 for q in pix):
            continue
        centre = recenter_lines(flux, [p], fwidth=fwidth, search=tol_pix)[0]
        if np.isfinite(centre):
            pix.append(float(centre))
            wav.append(float(list_waves[j]))
            nam.append(list_names[j])
    return np.asarray(pix), np.asarray(wav), nam


def refine_with_linelist(flux, solution: WaveSolution, list_waves, list_names=None,
                         fwidth: float = 10.0, tol_pix: float = 1.0,
                         degree: Optional[int] = None) -> WaveSolution:
    """Add catalogue lines to a solution and refit, keeping it only if still good.

    Templates can lack identified lines over wide wavelength ranges; filling
    those gaps from the line list constrains the polynomial there.
    """
    if list_waves is None or len(list_waves) == 0:
        return solution
    pix, wav, nam = add_lines_from_list(flux, solution, list_waves, list_names,
                                        fwidth=fwidth, tol_pix=tol_pix)
    if len(pix) <= solution.n_used:
        return solution
    new = fit_solution(pix, wav, npix=solution.npix,
                       degree=solution.degree if degree is None else degree, names=nam,
                       template=solution.template, lamp=solution.lamp, kind=solution.kind,
                       method=solution.method)
    limit = max(1.3 * solution.rms_pix, 0.3)
    if new.is_good() and new.n_used > solution.n_used and new.rms_pix <= limit:
        return new
    return solution


# ---------------------------------------------------------------------------
# Dispersion correction
# ---------------------------------------------------------------------------

def rebin_flux(edges_in, values, edges_out) -> np.ndarray:
    """Flux-conserving rebinning between two sets of increasing bin edges."""
    cum = np.concatenate([[0.0], np.cumsum(np.nan_to_num(values))])
    return np.diff(np.interp(edges_out, edges_in, cum))


def dispcor(spec: Spec1D, solution: WaveSolution) -> Spec1D:
    """Rebin a pixel-space spectrum onto a linear wavelength grid.

    Counts are conserved (as IRAF ``dispcor flux=yes``). The sigma band is
    rebinned in quadrature, which ignores the mild correlation between
    neighbouring output pixels.
    """
    n = spec.npix
    if n != solution.npix:
        raise ValueError(f"Spectrum has {n} pixels but the solution expects {solution.npix}.")
    edges = solution(np.arange(n + 1) - 0.5)
    centres = solution(np.arange(n))
    reverse = centres[0] > centres[-1]
    if reverse:
        edges = edges[::-1]
    if not np.all(np.diff(edges) > 0):
        raise ValueError("Wavelength solution is not monotonic; cannot resample.")
    w0, w1 = centres.min(), centres.max()
    dw = (w1 - w0) / (n - 1)
    out_edges = w0 + dw * (np.arange(n + 1) - 0.5)
    bands = []
    for i, band in enumerate(spec.bands):
        band = band[::-1] if reverse else band
        if i == 3:
            bands.append(np.sqrt(np.clip(rebin_flux(edges, band ** 2, out_edges), 0, None)))
        else:
            bands.append(rebin_flux(edges, band, out_edges))
    header = spec.header.copy()
    header["WAVERMS"] = (round(float(solution.rms), 4), "RMS of wavelength fit (A)")
    header["WAVENLIN"] = (solution.n_used, "Lines used in wavelength fit")
    return Spec1D(bands=bands, header=header, crval=float(w0), cdelt=float(dw),
                  calibrated=True)


# ---------------------------------------------------------------------------
# QA plot and interactive fallback
# ---------------------------------------------------------------------------

def plot_solution(flux, solution: WaveSolution, path=None, title: str = ""):
    """Lamp spectrum with labelled lines, plus fit residuals."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    flux = np.asarray(flux, dtype=float)
    fig = Figure(figsize=(11, 7))
    FigureCanvasAgg(fig)
    ax0, ax1 = fig.subplots(2, 1, gridspec_kw={"height_ratios": [2, 1]})
    waves = solution(np.arange(len(flux)))
    ax0.plot(waves, flux, lw=0.7, color="k")
    top = np.nanmax(flux) if len(flux) else 1.0
    for p, w, u in zip(solution.pixels, solution.wavelengths, solution.used):
        y = flux[int(np.clip(round(p), 0, len(flux) - 1))]
        ax0.annotate(f"{w:.1f}", (solution(p), y), xytext=(0, 12), textcoords="offset points",
                     rotation=90, ha="center", va="bottom", fontsize=6,
                     color="C0" if u else "C3",
                     arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "0.5"})
    ax0.set_ylim(-0.02 * top, 1.25 * top)
    ax0.set_ylabel("Counts")
    ax0.set_title(title or solution.summary(), fontsize=9)
    res = solution.residuals
    ax1.axhline(0, color="0.5", lw=0.8)
    ax1.plot(solution.wavelengths[solution.used], res[solution.used], "o", ms=4, label="used")
    if (~solution.used).any():
        ax1.plot(solution.wavelengths[~solution.used], res[~solution.used], "x",
                 color="C3", label="rejected")
    ax1.set_xlim(ax0.get_xlim())
    ax1.set_xlabel("Wavelength (A)")
    ax1.set_ylabel("Residual (A)")
    ax1.legend(fontsize=7)
    fig.tight_layout()
    if path is not None:
        fig.savefig(path)
    return fig


def identify_interactive(flux, list_waves=None, list_names=None,
                         initial: Optional[WaveSolution] = None, fwidth: float = 10.0,
                         degree: int = 4, title: str = "") -> Optional[WaveSolution]:
    """Mark lamp lines by hand in a matplotlib window.

    Keys (with the cursor over the plot):
      m  mark the line nearest the cursor; type its wavelength in the terminal
         (Enter accepts the suggested catalogue line once a fit exists)
      d  delete the marked line nearest the cursor
      f  fit and print the residuals
      l  add catalogue lines predicted by the current fit, then refit
      q  accept and close
    Returns the final solution, or None if fewer than two lines were marked.
    """
    import matplotlib.pyplot as plt

    flux = np.nan_to_num(np.asarray(flux, dtype=float))
    n = len(flux)
    state = {"pix": [], "wav": [], "nam": [], "sol": None}
    if initial is not None:
        state["pix"] = list(initial.pixels[initial.used])
        state["wav"] = list(initial.wavelengths[initial.used])
        state["nam"] = [nm for nm, u in zip(initial._names(), initial.used) if u]
        state["sol"] = initial
    have_list = list_waves is not None and len(list_waves) > 0

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(np.arange(n), flux, lw=0.8, color="k")
    ax.set_xlabel("Pixel")
    ax.set_ylabel("Counts")
    ax.set_title((title + "  " if title else "")
                 + "m: mark  d: delete  f: fit  l: add from list  q: accept", fontsize=9)
    artists: list = []

    def refit(report: bool = True):
        if len(state["pix"]) < 2:
            state["sol"] = None
            return
        state["sol"] = fit_solution(state["pix"], state["wav"], npix=n, degree=degree,
                                    names=state["nam"], method="interactive")
        if report:
            sol = state["sol"]
            print(sol.summary())
            for p, w, r, u in zip(sol.pixels, sol.wavelengths, sol.residuals, sol.used):
                print(f"  {p:9.2f}  {w:10.3f}  {r:+7.3f} {'' if u else ' (rejected)'}")

    def redraw():
        for art in artists:
            art.remove()
        artists.clear()
        for p, w in zip(state["pix"], state["wav"]):
            y = flux[int(np.clip(round(p), 0, n - 1))]
            artists.append(ax.axvline(p, color="C0", lw=0.6, alpha=0.6))
            artists.append(ax.text(p, y, f" {w:.1f}", rotation=90, fontsize=7,
                                   va="bottom", ha="center", color="C0"))
        fig.canvas.draw_idle()

    def on_key(event):
        if event.key == "q":
            plt.close(fig)
            return
        if event.key == "f":
            refit()
        elif event.key == "l":
            if state["sol"] is None or not have_list:
                print("Fit at least two lines (and provide a line list) first.")
            else:
                pix, wav, nam = add_lines_from_list(
                    flux, state["sol"], list_waves, list_names, fwidth=fwidth)
                state["pix"], state["wav"], state["nam"] = list(pix), list(wav), list(nam)
                refit()
        elif event.inaxes is ax and event.xdata is not None:
            if event.key == "m":
                centre = recenter_lines(flux, [event.xdata], fwidth=fwidth, min_snr=3.0)[0]
                if not np.isfinite(centre):
                    print(f"No line found near pixel {event.xdata:.1f}")
                    return
                suggestion, name = None, ""
                if state["sol"] is not None and have_list:
                    guess = float(state["sol"](centre))
                    j = int(np.argmin(np.abs(np.asarray(list_waves) - guess)))
                    suggestion = float(list_waves[j])
                    name = list_names[j] if list_names is not None else ""
                prompt = f"Wavelength of line at pixel {centre:.2f}"
                prompt += f" [{suggestion:.3f} {name}]: " if suggestion else ": "
                text = input(prompt).strip()
                try:
                    wave = float(text) if text else suggestion
                except ValueError:
                    wave = None
                if wave is None:
                    print("Skipped.")
                    return
                state["pix"].append(float(centre))
                state["wav"].append(float(wave))
                state["nam"].append(name if not text else "")
            elif event.key == "d" and state["pix"]:
                j = int(np.argmin(np.abs(np.asarray(state["pix"]) - event.xdata)))
                for key in ("pix", "wav", "nam"):
                    state[key].pop(j)
        redraw()

    fig.canvas.mpl_connect("key_press_event", on_key)
    redraw()
    plt.show()
    refit(report=False)
    return state["sol"]


# ---------------------------------------------------------------------------
# gwcs-based helper kept from the earlier API
# ---------------------------------------------------------------------------

class WaveCalibrator:
    """Fit pixel-wavelength pairs and expose the result as a gwcs object."""

    def __init__(self, line_pixels, line_wavelengths, wave_unit="Angstrom",
                 model="Legendre1D", degree=7, sigma=2, plot=False):
        from astropy import units as u

        self.line_pixels = line_pixels
        self.line_wavelengths = line_wavelengths
        self.wave_unit = getattr(u, wave_unit)
        self.model = model
        self.degree = degree
        self.sigma = sigma
        self.plot = plot
        self.calibrate()

    def calibrate(self):
        from astropy import units as u
        from gwcs import coordinate_frames as cf
        from gwcs import wcs

        from .utils import linear_fit_poly1D

        self.model_fit = linear_fit_poly1D(
            self.line_pixels, self.line_wavelengths, model=self.model,
            degree=self.degree, sigma=self.sigma, plot=self.plot)
        pixel_frame = cf.CoordinateFrame(1, "SPECTRAL", [0, ], axes_names=["x", ], unit=[u.pix, ])
        spectral_frame = cf.SpectralFrame(axes_names=["wavelength", ], unit=[self.wave_unit, ])
        self.wcs = wcs.WCS([(pixel_frame, self.model_fit), (spectral_frame, None)])

    def apply_to_spectrum(self, spectrum, resample=False):
        from specutils import Spectrum1D

        updated = Spectrum1D(spectrum.flux, wcs=self.wcs, mask=spectrum.mask,
                             uncertainty=spectrum.uncertainty)
        if updated.spectral_axis[0] > updated.spectral_axis[-1]:
            new_wave = updated.spectral_axis[::-1]
            new_flux = updated.flux[::-1]
            if updated.uncertainty is not None:
                updated.uncertainty.array = updated.uncertainty.array[::-1]
        else:
            new_wave = updated.spectral_axis
            new_flux = updated.flux
        return Spectrum1D(flux=new_flux, spectral_axis=new_wave,
                          uncertainty=updated.uncertainty)
