"""Telluric correction without IRAF (``continuum`` + ``telluric``).

The transmission spectrum is the standard star divided by a smooth
continuum inside the atmospheric bands (and exactly 1 elsewhere). For
each target it is shifted and scaled as ``T(lambda - shift) ** scale``
(Beer's law); shift and scale are chosen to leave the corrected spectrum
as smooth as possible across the bands, starting from the airmass ratio.
O2 and H2O get separate scale factors because water vapour varies
independently of airmass.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize

from .specio import Spec1D

# (species, lower, upper) in Angstrom
TELLURIC_BANDS: Tuple[Tuple[str, float, float], ...] = (
    ("O2", 6270.0, 6330.0),
    ("O2", 6860.0, 6960.0),
    ("H2O", 7150.0, 7350.0),
    ("O2", 7590.0, 7700.0),
    ("H2O", 8100.0, 8400.0),
    ("H2O", 8900.0, 9800.0),
)

# Stellar lines of hot standards that must not shape the continuum fit
# (H-alpha and the Paschen series): (centre, half width).
STELLAR_LINES: Tuple[Tuple[float, float], ...] = (
    (6562.8, 40.0), (8545.4, 20.0), (8598.4, 20.0), (8665.0, 20.0), (8750.5, 20.0),
    (8862.8, 20.0), (9014.9, 25.0), (9229.0, 25.0), (9546.0, 30.0), (10049.4, 30.0),
)

RED_LIMIT = {"LJT": 10000.0, "P200": 10000.0}
BLUE_LIMIT = 6750.0


def _band_mask(wave, bands=TELLURIC_BANDS) -> np.ndarray:
    mask = np.zeros(np.shape(wave), dtype=bool)
    for _, lo, hi in bands:
        mask |= (wave >= lo) & (wave <= hi)
    return mask


def fit_continuum(wave, flux, degree: int = 4, exclude: Optional[np.ndarray] = None,
                  low_reject: float = 2.0, high_reject: float = 3.0,
                  niter: int = 8) -> np.ndarray:
    """Smooth continuum from a sigma-clipped polynomial in log-log space.

    ``exclude`` marks pixels (telluric bands, stellar lines) that never
    enter the fit; the returned continuum is evaluated everywhere.
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    good = np.isfinite(flux) & (flux > 0)
    if exclude is not None:
        good &= ~exclude
    if good.sum() < degree + 2:
        raise ValueError("Too few continuum pixels to fit.")
    x = np.log(wave)
    y = np.full(flux.shape, np.nan)
    y[good] = np.log(flux[good])
    use = good.copy()
    poly = None
    for _ in range(niter):
        deg = int(min(degree, max(1, use.sum() // 20)))
        poly = np.polynomial.Chebyshev.fit(x[use], y[use], deg)
        res = y - poly(x)
        scatter = 1.4826 * np.nanmedian(np.abs(res[use] - np.nanmedian(res[use])))
        scatter = max(scatter, 1e-5)
        new_use = good & (res > -low_reject * scatter) & (res < high_reject * scatter)
        if new_use.sum() < deg + 2 or np.array_equal(new_use, use):
            break
        use = new_use
    return np.exp(poly(x))


def transmission(std_spec: Spec1D, telescope: Optional[str] = None,
                 blue_limit: float = BLUE_LIMIT, red_limit: Optional[float] = None,
                 degree: int = 4) -> Spec1D:
    """Atmospheric transmission from a standard-star spectrum.

    The continuum is fitted redward of ``blue_limit - 500`` A; the result is
    the normalised standard inside the telluric bands between ``blue_limit``
    and ``red_limit`` (8000 A, or 10000 A for LJT/P200) and 1 elsewhere.
    """
    if not std_spec.calibrated:
        raise ValueError("The standard-star spectrum must be wavelength calibrated.")
    if red_limit is None:
        red_limit = RED_LIMIT.get(telescope or "", 8000.0)
    wave = std_spec.wave
    flux = std_spec.flux
    region = wave >= blue_limit - 500.0
    trans = np.ones_like(flux)
    if region.sum() < 30:
        return Spec1D(bands=[trans], header=std_spec.header.copy(), crval=std_spec.crval,
                      cdelt=std_spec.cdelt, calibrated=True)
    exclude = _band_mask(wave)
    for centre, half in STELLAR_LINES:
        exclude |= np.abs(wave - centre) <= half
    cont = fit_continuum(wave[region], flux[region], degree=degree, exclude=exclude[region])
    ratio = np.ones_like(flux)
    ratio[region] = flux[region] / cont
    inside = _band_mask(wave) & (wave >= blue_limit) & (wave <= red_limit)
    trans[inside] = np.clip(ratio[inside], 0.02, 1.1)
    header = std_spec.header.copy()
    header["TELBLUE"] = (blue_limit, "blue limit of telluric template (A)")
    header["TELRED"] = (red_limit, "red limit of telluric template (A)")
    return Spec1D(bands=[trans], header=header, crval=std_spec.crval, cdelt=std_spec.cdelt,
                  calibrated=True)


def _roughness(wave, flux, regions: Sequence[np.ndarray]) -> float:
    """Summed squared residual about a local quadratic in each region."""
    total, count = 0.0, 0
    for sel in regions:
        y = flux[sel]
        x = wave[sel] - wave[sel].mean()
        ok = np.isfinite(y)
        if ok.sum() < 8:
            continue
        scale = np.nanmedian(np.abs(y[ok])) or 1.0
        coef = np.polyfit(x[ok], y[ok] / scale, 2)
        res = y[ok] / scale - np.polyval(coef, x[ok])
        total += float(np.sum(res ** 2))
        count += int(ok.sum())
    return total / count if count else np.inf


@dataclass
class TelluricFit:
    shift: float
    scales: Dict[str, float]
    model: np.ndarray  # transmission applied to the spectrum


def fit_telluric(wave, flux, trans_wave, trans, airmass_ratio: float = 1.0,
                 max_shift_pix: float = 3.0, scale_bounds: Tuple[float, float] = (0.2, 3.0),
                 fit: bool = True) -> TelluricFit:
    """Find the shift and per-species scale of the transmission for a target."""
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    step = abs(float(np.median(np.diff(wave))))
    x0 = float(np.clip(airmass_ratio, *scale_bounds))
    max_shift = max_shift_pix * step

    def species_model(species: str, shift: float, scale: float) -> np.ndarray:
        t = np.interp(wave - shift, trans_wave, trans, left=1.0, right=1.0)
        sel = _band_mask(wave, [b for b in TELLURIC_BANDS if b[0] == species])
        out = np.ones_like(t)
        out[sel] = np.clip(t[sel], 0.02, None) ** scale
        return out

    def regions_for(species: str) -> List[np.ndarray]:
        regions = []
        t0 = np.interp(wave, trans_wave, trans, left=1.0, right=1.0)
        for name, lo, hi in TELLURIC_BANDS:
            if name != species:
                continue
            inband = (wave >= lo) & (wave <= hi)
            if inband.sum() < 4 or np.min(t0[inband]) > 0.97:
                continue
            pad = min(0.5 * (hi - lo), 80.0)
            sel = (wave >= lo - pad) & (wave <= hi + pad)
            if sel.sum() >= 12:
                regions.append(sel)
        return regions

    shift = 0.0
    scales = {"O2": x0, "H2O": x0}
    if fit:
        o2_regions = regions_for("O2")
        if o2_regions:
            def cost_o2(p):
                return _roughness(wave, flux / species_model("O2", p[0], p[1]), o2_regions)
            res = minimize(cost_o2, [0.0, x0], method="Nelder-Mead",
                           bounds=[(-max_shift, max_shift), scale_bounds],
                           options={"xatol": 1e-3, "fatol": 1e-10, "maxiter": 400})
            if np.isfinite(res.fun) and res.fun <= cost_o2([0.0, x0]):
                shift, scales["O2"] = float(res.x[0]), float(res.x[1])
        h2o_regions = regions_for("H2O")
        if h2o_regions:
            def cost_h2o(p):
                return _roughness(wave, flux / species_model("H2O", shift, p[0]), h2o_regions)
            res = minimize(cost_h2o, [x0], method="Nelder-Mead", bounds=[scale_bounds],
                           options={"xatol": 1e-3, "fatol": 1e-10, "maxiter": 200})
            if np.isfinite(res.fun) and res.fun <= cost_h2o([x0]):
                scales["H2O"] = float(res.x[0])
    model = species_model("O2", shift, scales["O2"]) * species_model("H2O", shift, scales["H2O"])
    return TelluricFit(shift=shift, scales=scales, model=model)


def correct(spec: Spec1D, trans: Spec1D, airmass_ratio: float = 1.0,
            fit: bool = True) -> Tuple[Spec1D, TelluricFit]:
    """Divide a spectrum by the fitted transmission.

    The object, raw and sigma bands are corrected; the sky band is left alone.
    """
    result = fit_telluric(spec.wave, spec.flux, trans.wave, trans.flux,
                          airmass_ratio=airmass_ratio, fit=fit)
    out = spec.copy()
    for i in range(len(out.bands)):
        if i != 2:
            out.bands[i] = spec.bands[i] / result.model
    hdr = out.header
    hdr["TELSHIFT"] = (round(result.shift, 3), "telluric template shift (A)")
    hdr["TELSCO2"] = (round(result.scales["O2"], 3), "telluric O2 scale")
    hdr["TELSCH2O"] = (round(result.scales["H2O"], 3), "telluric H2O scale")
    return out, result


def plot_correction(spec: Spec1D, corrected: Spec1D, model: np.ndarray, path=None,
                    title: str = ""):
    """Before/after comparison over the red part of the spectrum."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(10, 6))
    FigureCanvasAgg(fig)
    ax0, ax1 = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    wave = spec.wave
    sel = wave >= 6200
    if sel.sum() < 10:
        sel = np.ones_like(wave, dtype=bool)
    ax0.plot(wave[sel], spec.flux[sel], lw=0.7, color="0.6", label="observed")
    ax0.plot(wave[sel], corrected.flux[sel], lw=0.8, color="k", label="corrected")
    ax0.legend(fontsize=7)
    ax0.set_ylabel("Flux")
    ax0.set_title(title, fontsize=9)
    ax1.plot(wave[sel], model[sel], lw=0.8)
    ax1.set_ylabel("Transmission")
    ax1.set_xlabel("Wavelength (A)")
    fig.tight_layout()
    if path is not None:
        fig.savefig(path)
    return fig
