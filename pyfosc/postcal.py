"""Post-extraction calibration: wavelength, flux and telluric (no IRAF).

Example (notebook, from the directory holding ``myfosc.json``)::

    from pyfosc.postcal import PostExtraction
    pe = PostExtraction('.')
    pe.identify_lamps()          # af*.fits  -> database/wavesol_*.json
    pe.apply_wavelength()        # acrf*.fits -> wacrf*.fits
    pe.build_sensfunc()          # standard star -> sens.fits
    pe.flux_calibrate()          # wacrf*.fits -> cwacrf*.fits
    pe.telluric_correct()        # cwacrf*.fits -> <OBJECT>_<tel><n>_<star>_tel.fits

Every step skips outputs that already exist unless ``overwrite=True``.
"""
from __future__ import annotations

import json
import logging
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from . import fluxcal, telluric
from . import newlog as log
from . import wavecal as wc
from .specio import Spec1D, read_spec, write_spec

logger = log.getLogger("pyfosc.postcal", level=logging.INFO)
logger.propagate = False

MIDNAME = {"XLT": "XL", "LJT": "LJ", "HCT": "HCT"}


class StandardNotFound(RuntimeError):
    """No spectrum could be matched to a tabulated standard star."""


def _can_interact() -> bool:
    """True when a terminal and an interactive matplotlib backend are available."""
    if not sys.stdin or not sys.stdin.isatty():
        return False
    import matplotlib

    return matplotlib.get_backend().lower() not in {"agg", "pdf", "svg", "ps", "cairo", "template"}


class PostExtraction:
    """Wavelength, flux and telluric calibration of extracted 1D spectra."""

    def __init__(self, work_dir: str | Path = ".", overwrite: bool = False):
        self.work_dir = Path(work_dir)
        data = self.work_dir / "data"
        # Accept both the project root (with ./data) and the data directory itself.
        self.data_dir = data if data.is_dir() and not list(self.work_dir.glob("a*.fit*")) \
            else self.work_dir
        self.overwrite = overwrite
        self.settings = self._load_settings()
        my = self.settings.get("mysettings", {})
        self.telescope: str = my.get("telescope", "")
        self.grism: str = my.get("Grism", "")
        if self.telescope == "P200":
            self.midname = "P200" + my.get("side", "")
        else:
            self.midname = MIDNAME.get(self.telescope, self.telescope)
        self.db_dir = self.data_dir / "database"
        self.qa_dir = (self.data_dir.parent if self.data_dir.name == "data"
                       else self.data_dir) / "QAplots"
        self.solutions: Dict[str, wc.WaveSolution] = {}
        self._session_templates: List[wc.Template] = []

    # -- setup ------------------------------------------------------------
    def _load_settings(self) -> Dict:
        for base in (self.data_dir, self.work_dir, self.data_dir.parent):
            path = base / "myfosc.json"
            if path.exists():
                return json.loads(path.read_text())
        raise FileNotFoundError(
            f"myfosc.json not found in {self.data_dir} or {self.work_dir}; run `pyfosc init`.")

    def _glob(self, pattern: str) -> List[Path]:
        return sorted(self.data_dir.glob(pattern))

    @property
    def lamp_files(self) -> List[Path]:
        return self._glob("af*.fit*")

    @property
    def science_files(self) -> List[Path]:
        return self._glob("acrf*.fit*")

    def _solution_path(self, lamp: Path) -> Path:
        return self.db_dir / f"wavesol_{lamp.stem}.json"

    # -- wavelength -------------------------------------------------------
    def identify_lamps(self, interactive: str = "auto", fwidth: float = 10.0,
                       degree: int = 4,
                       refine: Optional[bool] = None) -> Dict[str, wc.WaveSolution]:
        """Find a wavelength solution for every lamp spectrum.

        interactive : ``'auto'`` (packaged templates; open the matplotlib
            identifier only when that fails and a terminal is attached),
            ``'never'``, or ``'always'`` (mark the first lamp of each kind by
            hand, then transfer that solution to the other lamps).
        refine : refit on the line list after template matching. ``None``
            (default) does so only where a curated list exists (XLT Fe-Ar,
            from the BFOSC atlas); ``True`` also adds lines from the denser
            He-Ne-Ar list, which can introduce misidentifications, so check
            the QA plot; ``False`` keeps the template's own lines.
        """
        if interactive not in ("auto", "never", "always"):
            raise ValueError("interactive must be 'auto', 'never' or 'always'")
        lamps = self.lamp_files
        if not lamps:
            raise FileNotFoundError(f"No lamp spectra (af*.fit*) in {self.data_dir}.")
        self.qa_dir.mkdir(parents=True, exist_ok=True)
        done_by_hand = set()
        for lamp in lamps:
            sol_path = self._solution_path(lamp)
            spec = read_spec(lamp)
            kind = wc.lamp_kind(spec.header, self.telescope)
            if sol_path.exists() and not self.overwrite:
                sol = wc.WaveSolution.load(sol_path)
                if sol.npix == spec.npix:
                    self._register(lamp, spec, sol)
                    logger.info(f"Using saved solution {sol_path.name}")
                    continue
            templates = [t for t in self._session_templates if t.kind == kind]
            packaged = wc.list_templates(self.telescope, kind=kind, grism=self.grism)
            by_hand = interactive == "always" and kind not in done_by_hand
            sol = wc.auto_identify(spec.flux, templates + packaged, fwidth=fwidth, degree=degree)
            if sol is not None:
                sol.kind = kind
            good = sol is not None and sol.is_good()
            atlas = self.telescope == "XLT" and wc.has_atlas_linelist(kind)
            if good and refine is not False and atlas:
                sol = wc.refit_with_linelist(spec.flux, sol, *wc.load_linelist(kind),
                                             fwidth=fwidth, degree=degree)
            elif good and refine:
                sol = wc.refine_with_linelist(spec.flux, sol, *wc.load_linelist(kind),
                                              fwidth=fwidth, degree=degree)
            if by_hand or (not good and interactive != "never" and _can_interact()):
                if not good:
                    logger.warning(f"{lamp.name}: no reliable automatic solution; "
                                   "opening the interactive identifier.")
                waves, names = wc.load_linelist(kind)
                manual = wc.identify_interactive(
                    spec.flux, waves, names, initial=sol if good else None,
                    fwidth=fwidth, degree=degree, title=f"{lamp.name} ({kind})")
                if manual is not None:
                    sol, good = manual, manual.n_used >= 3 and manual.is_monotonic()
                    sol.kind = kind
                    done_by_hand.add(kind)
            if sol is None or not good:
                reason = ("no packaged template for this setup" if not packaged
                          else "template matching failed")
                detail = f" ({sol.summary()})" if sol is not None else ""
                logger.warning(f"{lamp.name}: {kind} lamp not calibrated, {reason}{detail}.")
                continue
            sol.lamp = lamp.name
            sol.save(sol_path)
            wc.plot_solution(spec.flux, sol, self.qa_dir / f"wavesol_{lamp.stem}.pdf")
            self._register(lamp, spec, sol)
            logger.info(sol.summary())
            if sol.max_gap > 0.2 * sol.npix:
                lo, hi = sol.gap_range()
                logger.warning(f"{lamp.name}: no lines between {lo:.0f} and {hi:.0f} A; "
                               "the wavelength scale is poorly constrained there.")
        if not self.solutions:
            raise RuntimeError(
                "No lamp spectrum could be wavelength calibrated. Run "
                "`python -m pyfosc.steps.identlamp2m` to identify lines by hand.")
        return self.solutions

    def _register(self, lamp: Path, spec: Spec1D, sol: wc.WaveSolution) -> None:
        self.solutions[lamp.name] = sol
        used = sol.used
        self._session_templates.append(wc.Template(
            name=lamp.stem, flux=spec.flux, pixels=sol.pixels[used],
            wavelengths=sol.wavelengths[used],
            names=[n for n, u in zip(sol._names(), used) if u],
            kind=sol.kind or wc.lamp_kind(spec.header, self.telescope)))

    def _load_solutions(self) -> None:
        if self.solutions:
            return
        for lamp in self.lamp_files:
            path = self._solution_path(lamp)
            if path.exists():
                self.solutions[lamp.name] = wc.WaveSolution.load(path)
        if not self.solutions:
            self.identify_lamps(interactive="auto")

    def _pick_solution(self, spec: Spec1D, name: str):
        """Solution of the lamp nearest in time (preferring the main lamp kind)."""
        if name in self.solutions and self.solutions[name].npix == spec.npix:
            return name, self.solutions[name]
        preferred = "hene" if self.telescope == "LJT" else "fear"
        t0 = fluxcal.get_obstime(spec.header)
        ranked = []
        for lamp_name, sol in self.solutions.items():
            if sol.npix != spec.npix:
                continue
            dt = np.inf
            lamp_path = self.data_dir / lamp_name
            if t0 is not None and lamp_path.exists():
                t1 = fluxcal.get_obstime(read_spec(lamp_path).header)
                if t1 is not None:
                    dt = abs((t1 - t0).sec)
            ranked.append((sol.kind != preferred, dt, lamp_name))
        if not ranked:
            raise RuntimeError(f"No wavelength solution with {spec.npix} pixels for {name}.")
        ranked.sort()
        best = ranked[0][2]
        return best, self.solutions[best]

    def apply_wavelength(self) -> List[Path]:
        """Dispersion-correct science and lamp spectra (``a*`` -> ``wa*``)."""
        self._load_solutions()
        outputs = []
        for path in self.science_files + self.lamp_files:
            out = path.with_name("w" + path.name)
            if out.exists() and not self.overwrite:
                outputs.append(out)
                continue
            spec = read_spec(path)
            try:
                lamp_name, sol = self._pick_solution(spec, path.name)
            except RuntimeError as exc:
                logger.warning(str(exc))
                continue
            cal = wc.dispcor(spec, sol)
            cal.header["REFSPEC1"] = Path(lamp_name).stem
            write_spec(out, cal.bands, header=cal.header, crval=cal.crval, cdelt=cal.cdelt,
                       calibrated=True)
            outputs.append(out)
            logger.info(f"{path.name} -> {out.name}  [{lamp_name}, "
                        f"{cal.crval:.1f}-{cal.wave[-1]:.1f} A, {cal.cdelt:.3f} A/pix]")
        return outputs

    # -- flux -------------------------------------------------------------
    @property
    def wavecal_science(self) -> List[Path]:
        return self._glob("wacrf*.fit*")

    def find_standards(self, caldir: Optional[str] = None) -> List[tuple]:
        """``(path, star, caldir)`` for every wavelength-calibrated standard."""
        found = []
        for path in self.wavecal_science:
            obj = str(read_spec(path).header.get("OBJECT", ""))
            entry = fluxcal.find_standard(obj, caldir=caldir)
            if entry is not None:
                found.append((path, entry[0], entry[1]))
        return found

    def checked_standards(self, caldir: Optional[str] = None) -> List[tuple]:
        """Standards of :meth:`find_standards` whose sensitivity curve is plausible.

        A frame whose curve has a different shape from the other standards of
        the night (wrong object in the slit, mislabelled frame) is dropped
        with a warning. With a single standard there is nothing to compare.
        """
        found = self.find_standards(caldir=caldir)
        if len(found) < 2:
            return found
        extinction = fluxcal.load_extinction(self.telescope)
        curves, usable = [], []
        for path, name, found_dir in found:
            spec = read_spec(path)
            try:
                table = fluxcal.load_standard(name, caldir=caldir or found_dir)
                sf = fluxcal.sensfunc(spec, table, extinction,
                                      airmass=fluxcal.get_airmass(spec.header, self.telescope),
                                      exptime=fluxcal.get_exptime(spec.header))
            except (KeyError, ValueError, FileNotFoundError) as exc:
                logger.warning(f"{path.name}: cannot be used as a standard ({exc}).")
                continue
            curves.append((sf.wave, sf.sens))
            usable.append((path, name, found_dir))
        flags = fluxcal.odd_standards(curves)
        for (path, name, _), odd in zip(usable, flags):
            if odd:
                logger.warning(
                    f"{path.name}: labelled {name}, but its sensitivity curve does not match "
                    "the other standard(s) of the night; not used as a standard.")
        return [item for item, odd in zip(usable, flags) if not odd]

    def build_sensfunc(self, std_file: Optional[str] = None, star: Optional[str] = None,
                       caldir: Optional[str] = None, order: int = 6) -> Path:
        """Fit the sensitivity function from a standard star and write ``sens.fits``."""
        sens_path = self.data_dir / "sens.fits"
        if sens_path.exists() and not self.overwrite and std_file is None and star is None:
            logger.info("Using existing sens.fits")
            return sens_path
        if std_file is not None:
            std_path = self.data_dir / Path(std_file).name
            if not std_path.exists():
                raise FileNotFoundError(f"Standard-star spectrum {std_path} not found.")
            name = star or str(read_spec(std_path).header.get("OBJECT", ""))
            table = fluxcal.load_standard(name, caldir=caldir)
        else:
            candidates = self.checked_standards(caldir=caldir)
            if star is not None:
                wanted = fluxcal.normalize_star_name(star)
                candidates = [c for c in candidates if c[1] == wanted]
            if not candidates:
                raise StandardNotFound(
                    "No wavelength-calibrated spectrum matches a tabulated standard star; "
                    "pass --std-file and --star.")
            if len(candidates) > 1:
                logger.info("Standards found: " + ", ".join(
                    f"{p.name} ({s})" for p, s, _ in candidates) + "; using the first.")
            std_path, name, found_dir = candidates[0]
            table = fluxcal.load_standard(name, caldir=caldir or found_dir)
        spec = read_spec(std_path)
        extinction = fluxcal.load_extinction(self.telescope)
        airmass = fluxcal.get_airmass(spec.header, self.telescope)
        exptime = fluxcal.get_exptime(spec.header)
        sf = fluxcal.sensfunc(spec, table, extinction, airmass=airmass, exptime=exptime,
                              order=order)
        header = spec.header.copy()
        header["STDNAME"] = (table.name, "standard star")
        header["STDFILE"] = (std_path.name, "standard-star spectrum")
        header["CALDIR"] = (table.caldir, "onedstds directory")
        header["STDAIRM"] = (round(airmass, 4), "airmass of the standard")
        header["SENSRMS"] = (round(sf.rms, 4), "rms of sensitivity fit (mag)")
        header["BUNIT"] = "2.5 log10(counts/s/A / (erg/cm2/s/A))"
        write_spec(sens_path, [sf.sens], header=header, crval=spec.crval, cdelt=spec.cdelt,
                   calibrated=True)
        self.qa_dir.mkdir(parents=True, exist_ok=True)
        fluxcal.plot_sensfunc(sf, self.qa_dir / "sensfunc.pdf")
        logger.info(f"Sensitivity function from {std_path.name} ({table.name}, {table.caldir}), "
                    f"{int(sf.points_used.sum())} points, rms {sf.rms:.3f} mag -> sens.fits")
        return sens_path

    def flux_calibrate(self) -> List[Path]:
        """Flux-calibrate all wavelength-calibrated science spectra (``w*`` -> ``cw*``)."""
        sens_path = self.data_dir / "sens.fits"
        if not sens_path.exists():
            self.build_sensfunc()
        sens = read_spec(sens_path)
        extinction = fluxcal.load_extinction(self.telescope)
        outputs = []
        for path in self.wavecal_science:
            out = path.with_name("c" + path.name)
            if out.exists() and not self.overwrite:
                outputs.append(out)
                continue
            spec = read_spec(path)
            try:
                exptime = fluxcal.get_exptime(spec.header)
            except KeyError as exc:
                logger.warning(f"{path.name}: {exc}; skipped.")
                continue
            airmass = fluxcal.get_airmass(spec.header, self.telescope)
            cal = fluxcal.calibrate(spec, sens.wave, sens.flux, extinction, airmass, exptime)
            cal.header["STDNAME"] = sens.header.get("STDNAME", "")
            write_spec(out, cal.bands, header=cal.header, crval=cal.crval, cdelt=cal.cdelt,
                       calibrated=True)
            outputs.append(out)
            logger.info(f"{path.name} -> {out.name}  [airmass {airmass:.2f}, {exptime:.0f} s]")
        return outputs

    # -- telluric ---------------------------------------------------------
    def _standard_for_telluric(self, std_file: Optional[str], star: Optional[str]):
        """Return ``(wavelength-calibrated standard path, star name)``."""
        if std_file is not None:
            path = self.data_dir / Path(std_file).name
            name = star or str(read_spec(path).header.get("OBJECT", "std"))
            return path, fluxcal.normalize_star_name(name)
        sens_path = self.data_dir / "sens.fits"
        if sens_path.exists():
            hdr = read_spec(sens_path).header
            path = self.data_dir / str(hdr.get("STDFILE", ""))
            if path.is_file():
                return path, str(hdr.get("STDNAME", star or "std"))
        candidates = self.find_standards()
        if not candidates:
            raise StandardNotFound(
                "No standard-star spectrum found for the telluric correction; "
                "pass --std-file and --star.")
        return candidates[0][0], candidates[0][1]

    def _telluric_standards(self, std_file: Optional[str], star: Optional[str],
                            prefix: str) -> List[dict]:
        """Standard exposures to build telluric templates from, primary first."""
        primary, name = self._standard_for_telluric(std_file, star)
        chosen = [(primary, name)]
        if std_file is None:
            for path, other, _ in self.checked_standards():
                if path != primary and (star is None
                                        or other == fluxcal.normalize_star_name(star)):
                    chosen.append((path, other))
        standards = []
        for path, starname in chosen:
            cal = path.with_name(prefix + path.name)
            use = cal if (prefix and cal.exists()) else path
            spec = read_spec(use)
            standards.append({
                "wfile": path.name, "path": use, "star": starname, "spec": spec,
                "calibrated": use == cal and bool(prefix),
                "airmass": fluxcal.get_airmass(spec.header, self.telescope),
                "time": fluxcal.get_obstime(spec.header)})
        return standards

    @staticmethod
    def _template_for(spec_time, standards: List[dict]):
        """``(index of the nearest standard, its template blended in water)``."""
        timed = [i for i, s in enumerate(standards) if s["time"] is not None]
        if spec_time is None or len(timed) < 2:
            return 0, standards[0]["template"].spec
        gaps = {i: (spec_time - standards[i]["time"]).to_value("s") for i in timed}
        near = min(timed, key=lambda i: abs(gaps[i]))
        other_side = [i for i in timed if i != near and gaps[i] * gaps[near] < 0]
        base = standards[near]
        if not other_side:
            return near, base["template"].spec
        far = min(other_side, key=lambda i: abs(gaps[i]))
        weight = abs(gaps[near]) / (abs(gaps[near]) + abs(gaps[far]))
        return near, telluric.blend_water(base["template"].spec, standards[far]["template"].spec,
                                          weight, base["airmass"], standards[far]["airmass"])

    def telluric_correct(self, std_file: Optional[str] = None, star: Optional[str] = None,
                         fit: bool = True, min_snr: float = 5.0, o2_exponent: float = 1.0,
                         red_limit: Optional[float] = None) -> List[Path]:
        """Remove telluric absorption using the standard star(s).

        Works on flux-calibrated spectra (``cwacrf*``) when present, otherwise
        on the wavelength-calibrated ones. Outputs are named
        ``<OBJECT>_<tel><n>_<star>_tel.fits``.

        A transmission template is built from every usable standard exposure;
        each target uses the one nearest in time (water interpolated in time
        when two standards bracket it). The scale of every band is measured
        on all targets first; their weighted mean per band is the night's
        value, which faint targets receive and bright targets are combined
        with (see :func:`pyfosc.telluric.fit_telluric`). ``fit=False`` applies
        the airmass law only. ``o2_exponent`` sets that law for O2
        (``airmass_ratio ** o2_exponent``); ``red_limit`` cuts the templates.
        """
        calibrated = bool(self._glob("cwacrf*.fit*"))
        prefix = "c" if calibrated else ""
        standards = self._telluric_standards(std_file, star, prefix)
        self.qa_dir.mkdir(parents=True, exist_ok=True)
        for i, std in enumerate(standards):
            # a flux-calibrated hot star has a near power-law continuum; raw counts
            # still carry the instrument response and need a more flexible global fit
            std["template"] = telluric.build_template(
                std["spec"], telescope=self.telescope, red_limit=red_limit,
                degree=4 if std["calibrated"] else 6)
            tag = std["star"] if i == 0 else f"{std['star']}_{i + 1}"
            tmpl = std["template"]
            write_spec(self.data_dir / f"telluric_{tag}.fits", tmpl.spec.bands,
                       header=tmpl.spec.header, crval=tmpl.spec.crval, cdelt=tmpl.spec.cdelt,
                       calibrated=True)
            telluric.plot_template(std["spec"], tmpl,
                                   self.qa_dir / f"telluric_template_{tag}.pdf",
                                   title=f"{std['path'].name} ({std['star']})")
            bands = ", ".join(f"{b['key']} ({b['ew']:.1f} A)" for b in tmpl.info
                              if b["status"] in ("local", "global"))
            dropped = ", ".join(f"{b['key']}: {b['status']}" for b in tmpl.info
                                if b["status"] in ("not significant", "no continuum"))
            logger.info(f"Telluric template from {std['path'].name}: {bands or 'no bands'}"
                        + (f"; left out: {dropped}" if dropped else ""))
        std_files = {s["path"].name for s in standards}
        targets = self._glob("cwacrf*.fit*") if calibrated else self.wavecal_science

        # first pass: measure every band on every target
        jobs = []
        for idx, path in enumerate(targets, start=1):
            spec = read_spec(path)
            which, trans = self._template_for(fluxcal.get_obstime(spec.header), standards)
            ratio = fluxcal.get_airmass(spec.header, self.telescope) / standards[which]["airmass"]
            first = telluric.fit_telluric(spec.wave, spec.flux, trans.wave, trans.flux,
                                          airmass_ratio=ratio, fit=fit, min_snr=min_snr,
                                          o2_exponent=o2_exponent)
            jobs.append({"idx": idx, "path": path, "spec": spec, "which": which,
                         "trans": trans, "ratio": ratio, "first": first})
        # the night's scale of each band, per standard, from the science targets
        night = {}
        if fit:
            for which in range(len(standards)):
                fits = [j["first"] for j in jobs
                        if j["which"] == which and j["path"].name not in std_files]
                night[which] = telluric.night_prior(fits)
                if night[which]:
                    logger.info(
                        f"Night-level scales relative to {standards[which]['star']}: "
                        + ", ".join(f"{k} x{v[0]:.2f} (+/-{v[1]:.2f} per target)"
                                    for k, v in night[which].items()))
                for key, (factor, _) in night[which].items():
                    if abs(factor - 1.0) > 0.15 and key.startswith("O2"):
                        logger.warning(
                            f"{key}: the targets need {factor:.2f}x the depth seen in "
                            f"{standards[which]['star']}. O2 does not vary, so the band of "
                            "the standard is probably filled in (second-order light of a "
                            "blue star); the night-level value is applied.")

        # second pass: correct, with the night's values as prior
        outputs, rows = [], []
        for job in jobs:
            spec, path, first = job["spec"], job["path"], job["first"]
            std = standards[job["which"]]
            obj = re.sub(r"[^\w.+-]+", "_", str(spec.header.get("OBJECT", "object")).strip())
            out = self.data_dir / f"{obj}_{self.midname}{job['idx']}_{std['star']}_tel.fits"
            prior = {key: (factor * first.band_law[key], sigma * first.band_law[key])
                     for key, (factor, sigma) in night.get(job["which"], {}).items()
                     if key in first.band_law}
            if path.name in std_files:
                prior = {}
            corrected, result = telluric.correct(
                spec, job["trans"], airmass_ratio=job["ratio"], fit=fit, min_snr=min_snr,
                o2_exponent=o2_exponent, prior=prior or None)
            for key, scale in result.band_scales.items():
                err = result.band_errors.get(key, np.nan)
                rows.append([out.name, obj, f"{job['ratio']:.4f}", f"{result.shift:.3f}", key,
                             f"{scale:.4f}", f"{err:.4f}" if np.isfinite(err) else "",
                             result.band_source[key], std["star"]])
            if out.exists() and not self.overwrite:
                outputs.append(out)
                continue
            corrected.header["TELSTD"] = (std["star"], "telluric standard")
            write_spec(out, corrected.bands, header=corrected.header, crval=corrected.crval,
                       cdelt=corrected.cdelt, calibrated=True)
            telluric.plot_correction(spec, corrected, result.model,
                                     self.qa_dir / f"telluric_{out.stem}.pdf", title=out.name,
                                     fit=result)
            outputs.append(out)
            detail = ", ".join(f"{k} x{v:.2f} ({result.band_source[k]})"
                               for k, v in result.band_scales.items()
                               if std["template"].used(k))
            logger.info(f"{path.name} -> {out.name}  [shift {result.shift:+.2f} A, {detail}]")
        if rows:
            import csv
            with open(self.qa_dir / "telluric_summary.csv", "w", newline="") as fh:
                writer = csv.writer(fh)
                writer.writerow(["file", "object", "airmass_ratio", "shift_A", "band", "scale",
                                 "error", "origin", "standard"])
                writer.writerows(rows)
        return outputs

    def run(self, interactive: str = "auto", std_file: Optional[str] = None,
            star: Optional[str] = None, caldir: Optional[str] = None) -> None:
        """Run all post-extraction steps in order."""
        self.identify_lamps(interactive=interactive)
        self.apply_wavelength()
        self.build_sensfunc(std_file=std_file, star=star, caldir=caldir)
        self.flux_calibrate()
        self.telluric_correct(std_file=std_file, star=star)
