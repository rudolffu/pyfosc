"""Telluric correction without IRAF (``continuum`` + ``telluric``).

The transmission spectrum is the standard star divided by a local
continuum in each atmospheric band (a power law through clean windows on
both sides of the band) and exactly 1 elsewhere. For each target it is
shifted and scaled as ``T(lambda - shift) ** scale`` (Beer's law). Scale
and shift are measured by a chi-square fit of ``T ** scale`` times a
smooth local continuum to the target, separately for O2 and H2O because
water vapour varies independently of airmass. A measured scale is used
only when it is significant; otherwise the airmass law (or a night-level
value supplied by the caller) is applied.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from .specio import Spec1D


@dataclass(frozen=True)
class Band:
    """A telluric band with the clean windows used as its local continuum."""

    species: str
    lo: float
    hi: float
    blue: Tuple[float, float]
    red: Tuple[float, float]

    @property
    def key(self) -> str:
        return f"{self.species}_{self.lo:.0f}"


BANDS: Tuple[Band, ...] = (
    Band("O2", 6270.0, 6330.0, (6180.0, 6265.0), (6335.0, 6450.0)),
    Band("O2", 6860.0, 6960.0, (6700.0, 6850.0), (6970.0, 7130.0)),
    Band("H2O", 7150.0, 7350.0, (6990.0, 7140.0), (7360.0, 7560.0)),
    Band("O2", 7590.0, 7700.0, (7400.0, 7580.0), (7710.0, 7880.0)),
    Band("H2O", 8100.0, 8400.0, (7900.0, 8090.0), (8410.0, 8640.0)),
    Band("H2O", 8900.0, 9800.0, (8700.0, 8890.0), (9810.0, 10020.0)),
)

# (species, lower, upper) in Angstrom
TELLURIC_BANDS: Tuple[Tuple[str, float, float], ...] = tuple(
    (b.species, b.lo, b.hi) for b in BANDS)

# Stellar lines of hot standards that must not shape the continuum fit
# (H-alpha and the Paschen series): (centre, half width).
STELLAR_LINES: Tuple[Tuple[float, float], ...] = (
    (6562.8, 40.0), (8545.4, 20.0), (8598.4, 20.0), (8665.0, 20.0), (8750.5, 20.0),
    (8862.8, 20.0), (9014.9, 25.0), (9229.0, 25.0), (9546.0, 30.0), (10049.4, 30.0),
)

# Reddest wavelength at which a band without clean windows on both sides
# may still be taken from the global continuum fit.
RED_LIMIT = {"LJT": 10000.0, "P200": 10000.0}
BLUE_LIMIT = 6750.0


def _band_mask(wave, bands=TELLURIC_BANDS) -> np.ndarray:
    mask = np.zeros(np.shape(wave), dtype=bool)
    for _, lo, hi in bands:
        mask |= (wave >= lo) & (wave <= hi)
    return mask


def _robust_std(values) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return np.nan
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


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


def band_continuum(wave, flux, band: Band, exclude: Optional[np.ndarray] = None,
                   min_pix: int = 8):
    """Local continuum of one band from the clean windows on both sides.

    A parabola in log-log (a straight line when a window has fewer than 15
    pixels) is fitted, with sigma clipping, to the pixels of the two flanking
    windows. Returns ``(region, continuum,
    scatter)`` where ``region`` selects the pixels from the blue edge of the
    blue window to the red edge of the red one, ``continuum`` is evaluated on
    them and ``scatter`` is the relative scatter of the flanks about the
    fit; ``None`` when either window has fewer than ``min_pix`` usable pixels.
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    good = np.isfinite(flux) & (flux > 0)
    if exclude is not None:
        good &= ~exclude
    blue = good & (wave >= band.blue[0]) & (wave <= band.blue[1])
    red = good & (wave >= band.red[0]) & (wave <= band.red[1])
    if blue.sum() < min_pix or red.sum() < min_pix:
        return None
    region = (wave >= band.blue[0]) & (wave <= band.red[1])
    flank = (blue | red)[region]
    degree = 2 if min(blue.sum(), red.sum()) >= 15 else 1
    cont = fit_continuum(wave[region], flux[region], degree=degree, exclude=~flank,
                         low_reject=2.5, high_reject=3.0)
    scatter = _robust_std(flux[region][flank] / cont[flank] - 1.0)
    return region, cont, scatter


@dataclass
class TelluricTemplate:
    """Transmission spectrum with what is needed to judge it."""

    spec: Spec1D                 # bands: transmission, its 1-sigma error
    continuum: np.ndarray        # continuum used inside each band region (NaN elsewhere)
    info: List[dict] = field(default_factory=list)   # one entry per band

    @property
    def wave(self) -> np.ndarray:
        return self.spec.wave

    @property
    def trans(self) -> np.ndarray:
        return self.spec.bands[0]

    def used(self, key: str) -> bool:
        return any(i["key"] == key and i["status"] in ("local", "global") for i in self.info)


def build_template(std_spec: Spec1D, telescope: Optional[str] = None,
                   blue_limit: float = BLUE_LIMIT, red_limit: Optional[float] = None,
                   degree: int = 4, min_sigma: float = 5.0,
                   stellar_lines: Sequence[Tuple[float, float]] = STELLAR_LINES
                   ) -> TelluricTemplate:
    """Atmospheric transmission from a standard-star spectrum.

    Each band redward of ``blue_limit`` gets a local continuum from its
    flanking windows (:func:`band_continuum`), so the transmission reaches 1
    at the band edges by construction. A band is kept only if its
    equivalent width exceeds ``min_sigma`` times the noise expected from the
    scatter of its flanks. Bands without usable windows on both sides fall
    back to a global continuum fit (polynomial of ``degree`` in log-log),
    and only blueward of 8000 A (10000 A for LJT/P200). ``red_limit`` cuts
    the template explicitly.
    """
    if not std_spec.calibrated:
        raise ValueError("The standard-star spectrum must be wavelength calibrated.")
    wave = std_spec.wave
    flux = np.asarray(std_spec.flux, dtype=float)
    step = abs(float(std_spec.cdelt))
    trans = np.ones_like(flux)
    err = np.zeros_like(flux)
    continuum = np.full(flux.shape, np.nan)
    stellar = np.zeros(flux.shape, dtype=bool)
    for centre, half in stellar_lines:
        stellar |= np.abs(wave - centre) <= half
    exclude = _band_mask(wave) | stellar
    global_limit = RED_LIMIT.get(telescope or "", 8000.0)
    global_cont = None
    info: List[dict] = []
    for band in BANDS:
        inside = (wave >= band.lo) & (wave <= band.hi)
        entry = {"key": band.key, "species": band.species, "lo": band.lo, "hi": band.hi,
                 "status": "uncovered", "ew": np.nan, "ew_err": np.nan,
                 "scatter": np.nan, "min": np.nan}
        info.append(entry)
        if inside.sum() < 4 or band.hi < blue_limit:
            if inside.sum() >= 4:
                entry["status"] = "blue of limit"
            continue
        if red_limit is not None and band.lo >= red_limit:
            entry["status"] = "red of limit"
            continue
        local = band_continuum(wave, flux, band, exclude=exclude)
        if local is not None:
            region, cont, scatter = local
            ratio = flux[region] / cont
            sel = inside[region]
            continuum[region] = cont
            status = "local"
        else:
            if band.lo >= global_limit:
                entry["status"] = "no continuum"
                continue
            if global_cont is None:
                fitreg = wave >= blue_limit - 500.0
                if fitreg.sum() < 30:
                    entry["status"] = "no continuum"
                    continue
                global_cont = np.full(flux.shape, np.nan)
                global_cont[fitreg] = fit_continuum(wave[fitreg], flux[fitreg], degree=degree,
                                                    exclude=exclude[fitreg])
            region = inside & np.isfinite(global_cont)
            if region.sum() < 4:
                entry["status"] = "no continuum"
                continue
            ratio = flux[region] / global_cont[region]
            sel = np.ones(region.sum(), dtype=bool)
            continuum[region] = global_cont[region]
            clean = np.isfinite(global_cont) & ~exclude
            scatter = _robust_std(flux[clean] / global_cont[clean] - 1.0)
            status = "global"
        values = ratio[sel]
        if red_limit is not None:
            values = np.where(wave[region][sel] <= red_limit, values, 1.0)
        scatter = scatter if np.isfinite(scatter) and scatter > 0 else 0.0
        ew = float(np.nansum(1.0 - values) * step)
        ew_err = float(scatter * np.sqrt(len(values)) * step)
        entry.update(ew=ew, ew_err=ew_err, scatter=float(scatter),
                     min=float(np.nanmin(values)))
        if ew <= min_sigma * ew_err:
            entry["status"] = "not significant"
            continue
        entry["status"] = status
        idx = np.flatnonzero(region)[sel]
        trans[idx] = np.clip(np.nan_to_num(values, nan=1.0), 0.02, 1.1)
        err[idx] = scatter
    header = std_spec.header.copy()
    header["TELBLUE"] = (blue_limit, "blue limit of telluric template (A)")
    kept = [i for i in info if i["status"] in ("local", "global")]
    header["TELRED"] = (max([i["hi"] for i in kept], default=blue_limit),
                        "red limit of telluric template (A)")
    header["TELBANDS"] = (" ".join(i["key"] for i in kept)[:68], "bands in the template")
    spec = Spec1D(bands=[trans, err], header=header, crval=std_spec.crval,
                  cdelt=std_spec.cdelt, calibrated=True)
    return TelluricTemplate(spec=spec, continuum=continuum, info=info)


def transmission(std_spec: Spec1D, telescope: Optional[str] = None,
                 blue_limit: float = BLUE_LIMIT, red_limit: Optional[float] = None,
                 degree: int = 4, min_sigma: float = 5.0) -> Spec1D:
    """Transmission spectrum of :func:`build_template` (band 0), with its error (band 1)."""
    return build_template(std_spec, telescope=telescope, blue_limit=blue_limit,
                          red_limit=red_limit, degree=degree, min_sigma=min_sigma).spec


@dataclass
class TelluricFit:
    """Shift and scales applied to a target, and where they came from.

    ``scales``, ``errors``, ``snr`` and ``source`` are keyed by species and
    describe its best-measured band; the ``band_*`` dictionaries hold every
    band (keys such as ``"O2_7590"``). ``source`` is ``"fit"`` or the label
    of the fallback (``"airmass"``, ``"night"``).
    """

    shift: float
    scales: Dict[str, float]
    model: np.ndarray  # transmission applied to the spectrum
    errors: Dict[str, float] = field(default_factory=dict)
    snr: Dict[str, float] = field(default_factory=dict)
    source: Dict[str, str] = field(default_factory=dict)
    shift_fitted: bool = False
    band_scales: Dict[str, float] = field(default_factory=dict)
    band_errors: Dict[str, float] = field(default_factory=dict)
    band_source: Dict[str, str] = field(default_factory=dict)
    band_law: Dict[str, float] = field(default_factory=dict)   # airmass-law value per band
    # raw (scale, error) of every band the model could describe, significant or not
    band_measured: Dict[str, Tuple[float, float]] = field(default_factory=dict)


class _BandFit:
    """Chi-square of ``T**scale`` times a smooth local continuum.

    Each selected band contributes the pixels of the band and of its
    two flanking windows. For a trial scale and shift the continuum (a
    quadratic) is solved linearly, so only scale and shift are searched.
    Noise is estimated per band from the target's own flanks.
    """

    def __init__(self, wave, flux, trans_wave, trans, species: Optional[str] = None,
                 other: Optional[np.ndarray] = None, degree: int = 2,
                 min_depth: float = 0.03, keys: Optional[Sequence[str]] = None):
        if other is None:
            other = np.ones_like(np.asarray(flux, dtype=float))
        self.wave, self.flux = wave, flux
        self.trans_wave, self.trans = trans_wave, trans
        self.regions: List[dict] = []
        t0 = np.interp(wave, trans_wave, trans, left=1.0, right=1.0)
        for band in BANDS:
            if species is not None and band.species != species:
                continue
            if keys is not None and band.key not in keys:
                continue
            inband = (wave >= band.lo) & (wave <= band.hi)
            if inband.sum() < 4 or np.min(t0[inband]) > 1.0 - min_depth:
                continue
            sel = (wave >= band.blue[0]) & (wave <= band.red[1]) & np.isfinite(flux)
            flank = sel & ~inband & ~_band_mask(wave)
            if flank.sum() < 2 * (degree + 1) or (sel & inband).sum() < 4:
                continue
            idx = np.flatnonzero(sel)
            x = wave[idx]
            xn = 2.0 * (x - x.min()) / max(np.ptp(x), 1e-9) - 1.0
            noise = _robust_std(np.diff(flux[np.flatnonzero(flank)])) / np.sqrt(2.0)
            if not np.isfinite(noise) or noise <= 0:
                continue
            self.regions.append({
                "key": band.key, "idx": idx, "inband": inband[idx], "noise": float(noise),
                "vander": np.polynomial.legendre.legvander(xn, degree),
                "other": other[idx], "good": np.ones(len(idx), dtype=bool)})

    def __bool__(self) -> bool:
        return bool(self.regions)

    def _power(self, region, scale: float, shift: float) -> np.ndarray:
        t = np.interp(self.wave[region["idx"]] - shift, self.trans_wave, self.trans,
                      left=1.0, right=1.0)
        out = np.ones_like(t)
        out[region["inband"]] = np.clip(t[region["inband"]], 0.02, None) ** scale
        return out * region["other"]

    def _residuals(self, region, scale: float, shift: float) -> np.ndarray:
        design = region["vander"] * self._power(region, scale, shift)[:, None]
        good = region["good"]
        coef, *_ = np.linalg.lstsq(design[good], self.flux[region["idx"]][good], rcond=None)
        return (self.flux[region["idx"]] - design @ coef) / region["noise"]

    def chi2(self, scale: float, shift: float = 0.0, keys: Optional[Sequence[str]] = None):
        total, dof = 0.0, 0
        for region in self.regions:
            if keys is not None and region["key"] not in keys:
                continue
            res = self._residuals(region, scale, shift)[region["good"]]
            total += float(np.sum(res ** 2))
            dof += len(res) - region["vander"].shape[1]
        return total, max(dof - 1, 1)

    def clip(self, scale: float, shift: float = 0.0, nsigma: float = 4.0, grow: int = 2,
             niter: int = 3):
        """Mask pixels far from the model (emission lines, cosmic rays).

        Call it at the best-fitting scale and shift. Residuals are judged
        against their own spread, separately inside the band and in the
        flanks, so that a template mismatch across the whole band of a bright
        target is not mistaken for outliers.
        """
        for region in self.regions:
            region["good"][:] = True
            inband = region["inband"]
            for _ in range(niter):
                res = self._residuals(region, scale, shift)
                bad = np.zeros(len(res), dtype=bool)
                for part in (inband, ~inband):
                    use = part & region["good"]
                    spread = max(1.0, _robust_std(res[use])) if use.sum() > 3 else 1.0
                    bad |= part & (np.abs(res) > nsigma * spread)
                bad = np.convolve(bad, np.ones(2 * grow + 1), mode="same") > 0
                if bad.sum() > 0.4 * len(bad) or np.array_equal(~bad, region["good"]):
                    break
                region["good"] = ~bad

    def kept_fraction(self) -> float:
        """Smallest fraction of in-band pixels that survived the clipping."""
        return min((float(r["good"][r["inband"]].mean()) for r in self.regions), default=0.0)

    def best_scale(self, bounds, shift: float = 0.0, keys=None):
        """``(scale, error, chi2, reduced chi2)`` minimising chi-square within ``bounds``."""
        res = minimize_scalar(lambda v: self.chi2(v, shift, keys)[0], bounds=bounds,
                              method="bounded", options={"xatol": 1e-4})
        scale = float(res.x)
        chi, dof = self.chi2(scale, shift, keys)
        h = 0.02 * max(scale, 0.1)
        curv = (self.chi2(scale + h, shift, keys)[0] - 2 * chi
                + self.chi2(scale - h, shift, keys)[0]) / h ** 2
        err = np.sqrt(2.0 / curv) if curv > 0 else np.inf
        red = chi / dof
        return scale, float(err * np.sqrt(max(red, 1.0))), chi, red

    def keys(self) -> List[str]:
        return [r["key"] for r in self.regions]


def _band_model(wave, trans_wave, trans, shift: float, band_scales: Dict[str, float]):
    t = np.interp(wave - shift, trans_wave, trans, left=1.0, right=1.0)
    out = np.ones_like(t)
    for band in BANDS:
        if band.key not in band_scales:
            continue
        sel = (wave >= band.lo) & (wave <= band.hi)
        out[sel] = np.clip(t[sel], 0.02, None) ** band_scales[band.key]
    return out


def fit_telluric(wave, flux, trans_wave, trans, airmass_ratio: float = 1.0,
                 max_shift_pix: float = 3.0, o2_range: Tuple[float, float] = (0.6, 1.7),
                 h2o_range: Tuple[float, float] = (0.4, 2.5), min_snr: float = 5.0,
                 o2_exponent: float = 1.0, fallback: Optional[Dict[str, float]] = None,
                 fallback_label: str = "airmass",
                 prior: Optional[Dict[str, Tuple[float, float]]] = None,
                 max_chi2: float = 25.0, min_shift_gain: float = 10.0,
                 fit: bool = True) -> TelluricFit:
    """Find the shift and the scale of each band of the transmission for a target.

    Every band is scaled as ``T ** scale``. The fallback is the airmass law
    with no shift: ``airmass_ratio ** o2_exponent`` for O2 bands and
    ``airmass_ratio`` for H2O bands, or the values in ``fallback`` (keyed by
    band, e.g. ``"O2_7590"``, or by species). For each band ``T ** scale``
    times a smooth local continuum is fitted to the target (chi-square, noise
    from the target's own flanking windows, outliers such as emission lines
    clipped at the best fit). The fitted scale replaces the fallback only if
    it is at least ``min_snr`` times its error and lies inside
    ``o2_range``/``h2o_range`` (relative to the airmass ratio), not on the
    edge; a band whose reduced chi-square exceeds ``max_chi2`` or that lost
    more than a quarter of its pixels to the clipping (the model cannot
    describe it, e.g. an emission line of the target sits in it) is never
    fitted. ``prior`` gives ``(value, sigma)`` per band from the
    rest of the night: the value is the fallback, and a usable measurement
    is combined with it by inverse-variance weighting instead of replacing
    it, so a noisy target stays close to the night's value and a bright one
    follows its own fit. Bands are scaled separately because the template can be wrong by
    a different factor in each (second-order light of a blue standard fills
    in the A band but not the B band). One shift is fitted to the usable O2
    bands together (their band heads are sharp); it is accepted when it
    lowers their chi-square, in units of the scatter about the fit, by more
    than ``min_shift_gain`` (half of that for shifts below one pixel) and
    is limited to ``max_shift_pix``.
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    step = abs(float(np.median(np.diff(wave))))
    x0 = float(np.clip(airmass_ratio, 0.2, 5.0))
    law = {"O2": x0 ** o2_exponent, "H2O": x0}
    ranges = {"O2": (o2_range[0] * x0, o2_range[1] * x0),
              "H2O": (h2o_range[0] * x0, h2o_range[1] * x0)}
    fallback = dict(fallback or {})
    prior = prior or {}
    if prior:
        fallback.update({k: v[0] for k, v in prior.items()})
        fallback_label = "night"
    covered = [b for b in BANDS
               if ((wave >= b.lo) & (wave <= b.hi)).sum() >= 4]
    band_law = {b.key: law[b.species] for b in covered}
    band_scales = {b.key: float(fallback.get(b.key, fallback.get(b.species, law[b.species])))
                   for b in covered}
    band_source = {b.key: (fallback_label if (b.key in fallback or b.species in fallback)
                           else "airmass") for b in covered}
    band_errors = {b.key: np.nan for b in covered}
    band_snr = {b.key: 0.0 for b in covered}
    band_measured: Dict[str, Tuple[float, float]] = {}
    shift, shift_fitted = 0.0, False
    if fit:
        max_shift = max_shift_pix * step
        fits = {}
        for band in covered:
            bf = _BandFit(wave, flux, trans_wave, trans, keys=[band.key])
            if bf:
                fits[band.key] = (band, bf)
        # first pass without shift: clip outliers and see which bands the model
        # can describe at all
        sound = {}
        for key, (band, bf) in fits.items():
            bounds = ranges[band.species]
            scale = bf.best_scale(bounds, 0.0)[0]
            bf.clip(scale, 0.0)
            red = bf.best_scale(bounds, 0.0)[3]
            if red <= max_chi2 and bf.kept_fraction() >= 0.75:
                sound[key] = (band, bf)
        # one shift for all bands, from the summed chi-square of the sound O2
        # bands: their band heads are sharp, while the shallow water bands
        # follow structure in the target's own continuum
        sharp = {k: v for k, v in sound.items() if v[0].species == "O2"}
        if max_shift > 0 and sharp:
            grid = np.arange(-max_shift, max_shift + 1e-9, 0.2 * step)
            chis = np.zeros(len(grid))
            at_zero, quality = 0.0, 1.0
            for key, (band, bf) in sharp.items():
                bounds = ranges[band.species]
                chis += np.array([bf.best_scale(bounds, g)[2] for g in grid])
                zero = bf.best_scale(bounds, 0.0)
                at_zero += zero[2]
                quality = max(quality, zero[3])
            i = int(np.argmin(chis))
            gain = (at_zero - chis[i]) / quality   # in units of the actual scatter
            small = abs(grid[i]) <= step
            if 0 < i < len(grid) - 1 and gain > (0.5 * min_shift_gain if small
                                                 else min_shift_gain):
                denom = chis[i - 1] - 2 * chis[i] + chis[i + 1]
                offset = (0.5 * (chis[i - 1] - chis[i + 1]) / denom * (grid[1] - grid[0])
                          if denom > 0 else 0.0)
                shift, shift_fitted = float(grid[i] + offset), True
        fitted_any = False
        for key, (band, bf) in sound.items():
            bounds = ranges[band.species]
            if shift_fitted:
                bf.clip(bf.best_scale(bounds, shift)[0], shift)
            scale, err, _, red = bf.best_scale(bounds, shift)
            band_errors[key] = err
            band_snr[key] = scale / err if err > 0 and np.isfinite(err) else 0.0
            lo, hi = bounds
            inside = min(scale - lo, hi - scale) > 0.02 * (hi - lo)
            usable = (inside and red <= max_chi2 and np.isfinite(err) and err > 0
                      and bf.kept_fraction() >= 0.75)
            if usable:
                band_measured[key] = (scale, err)
            if usable and key in prior and prior[key][1] > 0:
                value, sigma = prior[key]
                w_fit, w_prior = 1.0 / err ** 2, 1.0 / sigma ** 2
                band_scales[key] = float((scale * w_fit + value * w_prior) / (w_fit + w_prior))
                band_errors[key] = float(1.0 / np.sqrt(w_fit + w_prior))
                band_source[key] = "fit" if w_fit > w_prior else "night+fit"
                fitted_any = True
            elif usable and band_snr[key] >= min_snr:
                band_scales[key], band_source[key] = scale, "fit"
                fitted_any = True
        if not fitted_any:
            shift, shift_fitted = 0.0, False
    model = _band_model(wave, trans_wave, trans, shift, band_scales)
    scales, errors, snr, source = {}, {}, {}, {}
    for species in ("O2", "H2O"):
        keys = [b.key for b in covered if b.species == species]
        fitted = [k for k in keys if "fit" in band_source[k]]
        pick = (min(fitted, key=lambda k: band_errors[k]) if fitted
                else (max(keys, key=lambda k: band_snr[k]) if keys else None))
        scales[species] = band_scales[pick] if pick else law[species]
        errors[species] = band_errors[pick] if pick else np.nan
        snr[species] = band_snr[pick] if pick else 0.0
        source[species] = band_source[pick] if pick else "airmass"
    return TelluricFit(shift=shift, scales=scales, model=model, errors=errors, snr=snr,
                       source=source, shift_fitted=shift_fitted, band_scales=band_scales,
                       band_errors=band_errors, band_source=band_source, band_law=band_law,
                       band_measured=band_measured)


# Scatter of the scale factor from target to target that is real, not noise:
# small for O2 (constant column, only resolution effects), larger for water.
INTRINSIC_SCATTER = {"O2": 0.03, "H2O": 0.10}


def night_prior(fits: Sequence[TelluricFit], min_targets: int = 2,
                max_error: float = 0.15) -> Dict[str, Tuple[float, float]]:
    """Scale of each band relative to the airmass law, from all targets of a night.

    Returns ``{band: (factor, sigma)}``: the inverse-variance weighted mean of
    ``fitted scale / airmass law`` over the targets (3-sigma clipped) and the
    spread expected for a single target (measured excess scatter, at least
    :data:`INTRINSIC_SCATTER`). A band is included only if at least
    ``min_targets`` targets measured it and the mean is known to better than
    ``max_error``. Multiply by a target's airmass law to get its prior for
    :func:`fit_telluric`.
    """
    out: Dict[str, Tuple[float, float]] = {}
    for band in BANDS:
        values, errors = [], []
        for fit in fits:
            if band.key in fit.band_measured and fit.band_law.get(band.key, 0) > 0:
                scale, err = fit.band_measured[band.key]
                values.append(scale / fit.band_law[band.key])
                errors.append(err / fit.band_law[band.key])
        values, errors = np.asarray(values), np.asarray(errors)
        floor = INTRINSIC_SCATTER[band.species]
        keep = np.ones(len(values), dtype=bool)
        mean = err_mean = np.nan
        for _ in range(5):
            if keep.sum() < min_targets:
                break
            w = 1.0 / (errors[keep] ** 2 + floor ** 2)
            mean = float(np.sum(w * values[keep]) / np.sum(w))
            err_mean = float(1.0 / np.sqrt(np.sum(w)))
            new_keep = np.abs(values - mean) <= 3.0 * np.hypot(errors, floor)
            if np.array_equal(new_keep, keep):
                break
            keep = new_keep
        if keep.sum() < min_targets or not np.isfinite(mean) or err_mean > max_error:
            continue
        excess = np.mean((values[keep] - mean) ** 2) - np.mean(errors[keep] ** 2)
        scatter = float(np.sqrt(max(excess, 0.0)))
        out[band.key] = (mean, float(np.hypot(max(scatter, floor), err_mean)))
    return out


def blend_water(base: Spec1D, other: Spec1D, weight: float, airmass_base: float,
                airmass_other: float) -> Spec1D:
    """Template ``base`` with its H2O bands interpolated towards ``other``.

    ``weight`` is the share of ``other`` (0 to 1, e.g. from the observing
    times). The optical depths are brought to the airmass of ``base`` first,
    so the result is used with airmass ratios relative to ``base``. Bands
    that either template lacks are left as in ``base``.
    """
    weight = float(np.clip(weight, 0.0, 1.0))
    out = base.copy()
    wave = base.wave
    t_other = np.interp(wave, other.wave, other.bands[0], left=1.0, right=1.0)
    for band in BANDS:
        if band.species != "H2O":
            continue
        sel = (wave >= band.lo) & (wave <= band.hi)
        if sel.sum() < 4 or np.all(base.bands[0][sel] == 1.0) or np.all(t_other[sel] == 1.0):
            continue
        tau_base = -np.log(np.clip(base.bands[0][sel], 0.02, None))
        tau_other = -np.log(np.clip(t_other[sel], 0.02, None)) * airmass_base / airmass_other
        out.bands[0][sel] = np.exp(-((1.0 - weight) * tau_base + weight * tau_other))
    return out


def correct(spec: Spec1D, trans: Spec1D, airmass_ratio: float = 1.0,
            fit: bool = True, **fit_kwargs) -> Tuple[Spec1D, TelluricFit]:
    """Divide a spectrum by the fitted transmission.

    The object, raw and sigma bands are corrected; the sky band is left
    alone. When the template carries its own error (band 1), that error is
    added in quadrature to the sigma band. Extra keywords go to
    :func:`fit_telluric`.
    """
    result = fit_telluric(spec.wave, spec.flux, trans.wave, trans.flux,
                          airmass_ratio=airmass_ratio, fit=fit, **fit_kwargs)
    out = spec.copy()
    for i in range(len(out.bands)):
        if i != 2:
            out.bands[i] = spec.bands[i] / result.model
    if len(out.bands) > 3 and len(trans.bands) > 1:
        t = np.interp(spec.wave - result.shift, trans.wave, trans.bands[0], left=1.0, right=1.0)
        e = np.interp(spec.wave - result.shift, trans.wave, trans.bands[1], left=0.0, right=0.0)
        power = np.ones_like(t)
        for band in BANDS:
            if band.key in result.band_scales:
                power[(spec.wave >= band.lo) & (spec.wave <= band.hi)] = \
                    result.band_scales[band.key]
        rel = power * e / np.clip(t, 0.02, None)
        out.bands[3] = np.hypot(out.bands[3], np.abs(out.bands[0]) * rel)
    hdr = out.header
    hdr["TELSHIFT"] = (round(result.shift, 3), "telluric template shift (A)")
    hdr["TELSCO2"] = (round(result.scales["O2"], 3), "telluric O2 scale")
    hdr["TELSCH2O"] = (round(result.scales["H2O"], 3), "telluric H2O scale")
    for species, key in (("O2", "TELEO2"), ("H2O", "TELEH2O")):
        err = result.errors.get(species, np.nan)
        if np.isfinite(err):
            hdr[key] = (round(float(err), 3), f"error of fitted {species} scale")
    hdr["TELSRC"] = (" ".join(f"{k}:{v}" for k, v in result.source.items()),
                     "origin of the telluric scales")
    for i, (key, scale) in enumerate(result.band_scales.items(), start=1):
        err = result.band_errors.get(key, np.nan)
        text = f"{key} {scale:.3f} " + (f"{err:.3f}" if np.isfinite(err) else "-")
        hdr[f"TELBND{i}"] = (f"{text} {result.band_source[key]}", "band, scale, error, origin")
    return out, result


def plot_template(std_spec: Spec1D, template: TelluricTemplate, path=None, title: str = ""):
    """Standard star with the local continua, and the resulting transmission."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(11, 7))
    FigureCanvasAgg(fig)
    ax0, ax1 = fig.subplots(2, 1, sharex=True, gridspec_kw={"height_ratios": [3, 2]})
    wave = std_spec.wave
    sel = wave >= 6100
    if sel.sum() < 10:
        sel = np.ones_like(wave, dtype=bool)
    ax0.plot(wave[sel], std_spec.flux[sel], lw=0.7, color="k", label="standard")
    ax0.plot(wave[sel], template.continuum[sel], lw=1.0, color="C3", label="continuum used")
    colours = {"local": "C0", "global": "C1"}
    for band, entry in zip(BANDS, template.info):
        if not (wave[sel].min() < band.hi and wave[sel].max() > band.lo):
            continue
        used = entry["status"] in colours
        for ax in (ax0, ax1):
            ax.axvspan(band.lo, band.hi, color=colours.get(entry["status"], "0.6"),
                       alpha=0.15 if used else 0.08)
        if entry["status"] == "local":
            for lo, hi in (band.blue, band.red):
                ax0.axvspan(lo, hi, ymin=0.0, ymax=0.04, color="C2")
        label = f"{band.key}\n{entry['status']}"
        if np.isfinite(entry["ew"]):
            label += f"\nEW {entry['ew']:.1f}+/-{entry['ew_err']:.1f} A"
        ax1.text(0.5 * (band.lo + band.hi), 0.03, label, ha="center", va="bottom", fontsize=6,
                 transform=ax1.get_xaxis_transform())
    ax0.legend(fontsize=7)
    ax0.set_ylabel("Flux")
    ax0.set_title(title + "   (green: continuum windows)", fontsize=9)
    ax1.plot(wave[sel], template.trans[sel], lw=0.8, color="k")
    ax1.set_ylabel("Transmission")
    ax1.set_xlabel("Wavelength (A)")
    fig.tight_layout()
    if path is not None:
        fig.savefig(path)
    return fig


def plot_correction(spec: Spec1D, corrected: Spec1D, model: np.ndarray, path=None,
                    title: str = "", fit: Optional[TelluricFit] = None):
    """Before/after comparison over the red part of the spectrum, with band zooms."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(figsize=(11, 8))
    FigureCanvasAgg(fig)
    grid = fig.add_gridspec(3, 2, height_ratios=[2, 1, 1.6])
    ax0 = fig.add_subplot(grid[0, :])
    ax1 = fig.add_subplot(grid[1, :], sharex=ax0)
    wave = spec.wave
    sel = wave >= 6200
    if sel.sum() < 10:
        sel = np.ones_like(wave, dtype=bool)
    ax0.plot(wave[sel], spec.flux[sel], lw=0.7, color="0.6", label="observed")
    ax0.plot(wave[sel], corrected.flux[sel], lw=0.8, color="k", label="corrected")
    ax0.legend(fontsize=7)
    ax0.set_ylabel("Flux")
    if fit is not None:
        parts = [f"shift {fit.shift:+.2f} A"]
        for key, scale in fit.band_scales.items():
            if np.all(model[(wave >= float(key.split("_")[1]))][:3] == 1.0) and \
                    fit.band_source[key] == "airmass" and key not in fit.band_measured:
                continue
            err = fit.band_errors.get(key, np.nan)
            parts.append(f"{key} x{scale:.2f}" + (f"+/-{err:.2f}" if np.isfinite(err) else "")
                         + f" ({fit.band_source[key]})")
        title = title + "\n" + ", ".join(parts)
    ax0.set_title(title, fontsize=8)
    ax1.plot(wave[sel], model[sel], lw=0.8)
    ax1.set_ylabel("Transmission")
    ax1.set_xlabel("Wavelength (A)")
    for column, (lo, hi, name) in enumerate(((6780.0, 7040.0, "O2 B band"),
                                             (7480.0, 7800.0, "O2 A band"))):
        ax = fig.add_subplot(grid[2, column])
        zoom = (wave >= lo) & (wave <= hi)
        if zoom.sum() > 5:
            ax.plot(wave[zoom], spec.flux[zoom], lw=0.7, color="0.6")
            ax.plot(wave[zoom], corrected.flux[zoom], lw=0.9, color="k")
        ax.set_title(name, fontsize=8)
        ax.set_xlabel("Wavelength (A)")
    fig.tight_layout()
    if path is not None:
        fig.savefig(path)
    return fig
