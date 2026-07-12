from __future__ import annotations

"""
High-level pipeline API for use from notebooks and CLI.

Example (Notebook):
    from pyfosc.pipeline import PreWaveCal
    pw = PreWaveCal('.')              # working directory that has myfosc.json and FITS
    pw.discover()                     # classify files, set parameters
    pw.build_master_bias()
    pw.calibrate_bias_and_trim()      # applies to flats, lamps, science
    pw.build_master_flat()
    pw.normalize_flat(degree=25, gauss_sigma=20)
    pw.apply_flat_correction()
    pw.cosmic_ray_clean()
    pw.extract_1d()
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Dict
import json
import re
import warnings

import numpy as np
import ccdproc as ccdp
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.io import fits
from astropy import units as u

from .fosc import FOSCFileCollection
from .masters import MasterBias, MasterFlat, FlatNormalizer
from .extraction import Extract1dSpec
from .project import (
    MANIFEST_NAME,
    ProjectPaths,
    ProvenanceStore,
    atomic_ccd_write,
    prepare_iraf_handoff,
    validate_manifest,
)


@dataclass
class PipelineParams:
    trimsec: Optional[str] = None
    biassec: Optional[str] = None
    disp_axis: int = 0
    overscan: bool = False


class PreWaveCal:
    def __init__(self, work_dir: str | Path = '.', run_group: str | None = None):
        self.work_dir = Path(work_dir)
        self.settings = self._load_settings(self.work_dir / 'myfosc.json')
        self.paths = ProjectPaths.from_settings(self.work_dir, self.settings, run_group)
        self.paths.ensure()
        self.data_dir = self.paths.input
        self.staged = self.paths.mode == 'staged'
        if not self.staged:
            warnings.warn(
                "Legacy mixed data/ layout is deprecated; run "
                "'pyfosc migrate-layout' to preview migration.",
                FutureWarning,
                stacklevel=2,
            )
        self.manifest = None
        self.group_entries = []
        self.provenance = ProvenanceStore(self.paths)
        self.params = PipelineParams()
        self.ic_all: Optional[FOSCFileCollection] = None
        self.master_bias: Optional[CCDData] = None
        self.master_flat: Optional[CCDData] = None
        self.normalized_flat: Optional[CCDData] = None
        # Suffix for bias-corrected, trimmed frames (keep original untrimmed)
        self.bias_suffix = "_bc"

    @staticmethod
    def _load_settings(path: Path) -> Dict:
        with open(path) as f:
            return json.load(f)

    def discover(self) -> None:
        filenames = None
        if self.staged:
            manifest_path = self.work_dir / MANIFEST_NAME
            if not manifest_path.exists():
                raise FileNotFoundError(f"Staged manifest is missing: {manifest_path}")
            self.manifest = validate_manifest(
                self.work_dir, json.loads(manifest_path.read_text())
            )
            self.group_entries = [
                item for item in self.manifest['entries']
                if item['group'] == self.paths.group and item['role'] != 'other'
            ]
            if not self.group_entries:
                raise ValueError(f"Unknown or empty run group: {self.paths.group}")
            filenames = [item['canonical_name'] for item in self.group_entries]
        ic = FOSCFileCollection(location=str(self.data_dir), filenames=filenames)
        if ic.summary is None:
            raise FileNotFoundError(
                f"No FITS files found under {self.data_dir}. Ensure initialization "
                "populated the configured input directory."
            )
        ic.check_groups()
        if self.staged:
            by_role = {
                role: [item['canonical_name'] for item in self.group_entries if item['role'] == role]
                for role in ('bias', 'flat', 'lamp', 'science')
            }
            ic.list_bias = by_role['bias']
            ic.list_flat = by_role['flat']
            ic.list_cal = by_role['lamp']
            ic.list_sci = by_role['science']
            ic.list_allbutbias = ic.list_flat + ic.list_cal + ic.list_sci
            ic.list_sci_cal = ic.list_cal + ic.list_sci
        if not ic.list_bias:
            raise ValueError('No bias frames match the selected science-frame geometry')
        if self.staged:
            missing_roles = [
                role for role, values in (
                    ('flat', ic.list_flat), ('lamp', ic.list_cal), ('science', ic.list_sci)
                ) if not values
            ]
            if missing_roles:
                raise ValueError(
                    f"Run group {self.paths.group} lacks required roles: {missing_roles}"
                )
        ic.set_parameters()
        overrides = self.settings.get('prewavecal', {})
        if not isinstance(overrides, dict):
            raise ValueError("myfosc.json 'prewavecal' must be an object")
        allowed = {'trimsec', 'biassec', 'disp_axis', 'overscan'}
        unknown = set(overrides) - allowed
        if unknown:
            raise ValueError(f"Unknown prewavecal settings: {sorted(unknown)}")
        resolved = dict(ic.parameters)
        resolved.update(overrides)
        ic.parameters = resolved
        self.ic_all = ic
        p = PipelineParams(
            trimsec=resolved.get('trimsec'),
            biassec=resolved.get('biassec'),
            disp_axis=resolved.get('disp_axis', 0),
            overscan=resolved.get('overscan', False),
        )
        self.params = p
        self._validate_input_geometry()
        self._report_section_disagreements()
        print(f"Resolved prewavecal parameters: {self.params}")

    # Master Bias
    def build_master_bias(self, method: str = 'average', save: bool = True) -> CCDData:
        assert self.ic_all is not None, 'Call discover() first'
        master_path = self.paths.masters / 'master_bias.fits'
        if self.staged:
            calibrated = [
                self.paths.stage('calibrated') / name
                for name in self._filter_base(self.ic_all.list_allbutbias)
            ]
            if self._stage_is_current('calibrated', calibrated + [master_path]):
                self.master_bias = CCDData.read(master_path, unit='adu')
                return self.master_bias
        mb = MasterBias(
            work_dir=str(self.work_dir), list_bias=self.ic_all.list_bias,
            data_dir=self.data_dir, output_dir=self.paths.masters,
        )
        # Build without saving, then store under data/MasterFrames
        self.master_bias = mb.build(save=False, method=method)
        if save:
            outdir = self.paths.masters
            outdir.mkdir(parents=True, exist_ok=True)
            if self.staged:
                atomic_ccd_write(self.master_bias, master_path)
            else:
                self.master_bias.write(outdir / 'master_bias.fits', overwrite=True)
        return self.master_bias

    def calibrate_bias_and_trim(self) -> None:
        assert self.master_bias is not None and self.ic_all is not None
        params = self.params
        # Prepare a reference bias with same processing as the science frames
        bias_ref = self.master_bias
        if params.overscan and params.biassec:
            bias_ref = ccdp.subtract_overscan(bias_ref, fits_section=params.biassec, overscan_axis=1)
        if params.trimsec:
            bias_ref = ccdp.trim_image(bias_ref, fits_section=params.trimsec)
        # Only process base files; skip derived outputs from prior runs
        base_files = self._filter_base(self.ic_all.list_allbutbias)
        ic = ImageFileCollection(location=str(self.data_dir), filenames=base_files)
        output_dir = self.paths.stage('calibrated')
        stage_outputs = [output_dir / f for f in base_files]
        provenance_outputs = stage_outputs + [self.paths.masters / 'master_bias.fits']
        if self.staged and self._stage_is_current('calibrated', provenance_outputs):
            return
        for ccd, fname in ic.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            outname = fname if self.staged else self._add_suffix(fname, self.bias_suffix)
            outpath = output_dir / outname
            # Raw acquisition headers commonly contain TRIMSEC. Only a valid
            # pyfosc-derived output is evidence that this step is complete.
            if (
                not self.staged
                and outpath.exists()
                and self._processed_output_is_current(outpath)
            ):
                continue
            if params.overscan and params.biassec:
                ccd = ccdp.subtract_overscan(ccd, fits_section=params.biassec, overscan_axis=1)
            if params.trimsec:
                ccd = ccdp.trim_image(ccd, fits_section=params.trimsec)
            ccd = ccdp.subtract_bias(ccd, bias_ref)
            ccd.header['SUBTRACT_BIAS'] = 'subbias'
            ccd.header['PYFOSC'] = ('PREWAVECAL', 'pyfosc processing stage')
            if params.trimsec:
                ccd.header['TRIMSEC'] = params.trimsec
            if params.biassec and params.overscan:
                ccd.header['BIASSEC'] = params.biassec
            if self.staged:
                atomic_ccd_write(ccd, outpath)
            else:
                ccd.write(outpath, overwrite=True)
        if self.staged:
            self._write_stage_provenance(
                'calibrated', provenance_outputs
            )

    # Master Flat and normalization
    def build_master_flat(self, method: str = 'average', save: bool = True) -> CCDData:
        assert self.ic_all is not None
        # Prefer bias-corrected, trimmed flats if present
        flat_bc = [f if self.staged else self._add_suffix(f, self.bias_suffix) for f in self.ic_all.list_flat]
        calibrated_dir = self.paths.stage('calibrated')
        flat_list = [f for f in flat_bc if (calibrated_dir / f).exists()]
        missing = sorted(set(flat_bc) - set(flat_list))
        if missing:
            raise RuntimeError(
                "Processed flat inputs are missing; run calibrate_bias_and_trim() first: "
                + ", ".join(missing)
            )
        mf = MasterFlat(
            work_dir=str(self.work_dir), list_flat=flat_list,
            data_dir=calibrated_dir, output_dir=self.paths.masters,
        )
        # Build without saving, then store under data/MasterFrames
        self.master_flat = mf.build(save=False, method=method)
        if save:
            outdir = self.paths.masters
            outdir.mkdir(parents=True, exist_ok=True)
            if self.staged:
                atomic_ccd_write(self.master_flat, outdir / 'master_flat.fits')
            else:
                self.master_flat.write(outdir / 'master_flat.fits', overwrite=True)
        return self.master_flat

    def normalize_flat(self, degree: int = 25, gauss_sigma: int = 20, save: bool = True) -> CCDData:
        assert self.master_flat is not None
        self._flat_normalization_params = {
            'degree': degree,
            'gauss_sigma': gauss_sigma,
        }
        # Save normalized outputs under MasterFrames at repo root
        fnorm = FlatNormalizer(
            master_flat=self.master_flat, disp_axis=self.params.disp_axis,
            work_dir=str(self.work_dir), output_dir=self.paths.masters,
        )
        fnorm.response_correct(degree=degree)
        fnorm.illumination_correct(gauss_sigma=gauss_sigma, save=False)
        self.normalized_flat = fnorm.normalized_flat
        if save:
            outpath = self.paths.masters / 'normalized_flat.fits'
            if self.staged:
                atomic_ccd_write(self.normalized_flat, outpath)
            else:
                self.normalized_flat.write(outpath, overwrite=True)
        return self.normalized_flat

    # Apply flat correction
    def apply_flat_correction(self) -> None:
        assert self.ic_all is not None and self.normalized_flat is not None
        outputs = self._flat_stage_outputs() if self.staged else [
            self.paths.stage('flat_corrected') / f
            for f in self.ic_all.list_sci + self.ic_all.list_cal
        ]
        if self.staged and self._stage_is_current(
            'flat_corrected', outputs,
            getattr(self, '_flat_normalization_params', {}), 'calibrated'
        ):
            return
        for fname in self.ic_all.list_sci + self.ic_all.list_cal:
            # Apply flat to bias-corrected, trimmed frames (if present)
            inname = fname if self.staged else self._add_suffix(fname, self.bias_suffix)
            inpath = self.paths.stage('calibrated') / inname
            if not inpath.exists():
                raise FileNotFoundError(
                    f"Processed input {inpath} is missing; run calibrate_bias_and_trim() first"
                )
            ccd = CCDData.read(inpath, unit='adu')
            if ccd.shape != self.normalized_flat.shape:
                raise ValueError(
                    f"Shape mismatch for {inpath.name}: image {ccd.shape}, "
                    f"normalized flat {self.normalized_flat.shape}"
                )
            outpath = self.paths.stage('flat_corrected') / (
                fname if self.staged else 'f' + inpath.name
            )
            if outpath.exists() and not self.staged:
                continue
            if 'subbias' not in ccd.header.get('SUBTRACT_BIAS', '').lower():
                raise RuntimeError(f"{inpath.name} is not bias-subtracted")
            if 'FLAT_CORRECT' not in ccd.header:
                flat_corrected = ccdp.flat_correct(ccd, self.normalized_flat)
                flat_corrected.header['FLAT_CORRECT'] = 'flatcor'
                if self.staged:
                    atomic_ccd_write(flat_corrected, outpath)
                else:
                    flat_corrected.write(outpath, overwrite=True)
        if self.staged:
            self._write_stage_provenance(
                'flat_corrected', outputs,
                getattr(self, '_flat_normalization_params', {}), 'calibrated'
            )

    # Cosmic ray cleaning
    def cosmic_ray_clean(self, sigclip: float = 7.0) -> None:
        assert self.ic_all is not None
        outputs = [
            self.paths.stage('cosmic_cleaned') / f for f in self.ic_all.list_sci
        ]
        if self.staged and self._stage_is_current(
            'cosmic_cleaned', outputs, {'sigclip': sigclip}, 'flat_corrected'
        ):
            return
        # telescope-based defaults if header is missing
        tel = self.ic_all.telescope
        defaults = {'XLT': (2.2, 7.8), 'LJT': (1.1, 6.6), 'HCT': (1.22, 4.8)}
        d_gain, d_rd = defaults.get(tel, (1.0, 5.0))
        for fname in self._filter_base(self.ic_all.list_sci):
            # work from flat-corrected of bias-corrected trimmed inputs
            inname = fname if self.staged else 'f' + self._add_suffix(fname, self.bias_suffix)
            fpath = self.paths.stage('flat_corrected') / inname
            if not fpath.exists():
                fpath = Path(self.data_dir) / ('f' + fname)
            if not fpath.exists():
                raise FileNotFoundError(
                    f"Flat-corrected input for {fname} is missing; "
                    "run apply_flat_correction() first"
                )
            outpath = self.paths.stage('cosmic_cleaned') / (
                fname if self.staged else 'cr' + fpath.name
            )
            if outpath.exists() and not self.staged:
                continue
            ccd = CCDData.read(fpath, unit='adu')
            hdr = ccd.header
            # Prefer values from headers (case-insensitive, multiple aliases)
            gain = self._first_present(hdr, ['GAIN', 'gain']) or d_gain
            rd = self._first_present(hdr, ['RDNOISE', 'rdnoise', 'RON', 'ron', 'READNOISE', 'readnoise']) or d_rd
            new_ccd = ccdp.cosmicray_lacosmic(
                ccd,
                gain=gain,
                readnoise=rd,
                sigclip=sigclip,
                verbose=True,
            )
            # Ensure a safe uncertainty is present for downstream extraction
            # Do not attach a unit to the StdDevUncertainty per request
            new_ccd.mask = None
            # Drop existing uncertainty to avoid unit-mismatch propagation
            try:
                new_ccd.uncertainty = None
            except Exception:
                pass
            std = np.sqrt(np.clip(np.abs(new_ccd.data), 0, None))
            new_ccd.uncertainty = StdDevUncertainty(std)
            # Write alongside inputs under data_dir with 'cr' prefix
            if self.staged:
                atomic_ccd_write(new_ccd, outpath)
            else:
                new_ccd.write(outpath, overwrite=True)
        if self.staged:
            self._write_stage_provenance('cosmic_cleaned', outputs, {'sigclip': sigclip}, 'flat_corrected')

    # 1D extraction
    def extract_1d(self, guess: float | int | None = None) -> None:
        assert self.ic_all is not None
        if self.staged:
            science = list(self.ic_all.list_sci)
            lamps = list(self.ic_all.list_cal)
            outputs = [self.paths.stage('extracted') / f for f in science + lamps]
            if self._stage_is_current(
                'extracted', outputs, {'guess': guess}, 'cosmic_cleaned'
            ):
                prepare_iraf_handoff(self.paths, science, lamps)
                return
            missing = [f for f in science if not (self.paths.stage('cosmic_cleaned') / f).exists()]
            missing += [f for f in lamps if not (self.paths.stage('flat_corrected') / f).exists()]
            if missing:
                raise FileNotFoundError(f"Extraction inputs are missing: {missing}")
            science_ic = FOSCFileCollection(
                location=str(self.paths.stage('cosmic_cleaned')), filenames=science
            )
            extractor = Extract1dSpec(
                ic_2dspec=science_ic, disp_axis=self.params.disp_axis,
                qa_dir=self.paths.qa, output_dir=self.paths.stage('extracted'),
            )
            extractor.trace_and_extract(guess=guess)
            lamp_ic = FOSCFileCollection(
                location=str(self.paths.stage('flat_corrected')), filenames=lamps
            )
            extractor.extract_lamp(ic_lamp=lamp_ic)
            self._write_stage_provenance('extracted', outputs, {'guess': guess}, 'cosmic_cleaned')
            prepare_iraf_handoff(self.paths, science, lamps)
            return
        # Objects (resolve actual files on disk; avoid re-adding prefixes)
        crf_list = []
        for x in self._filter_base(self.ic_all.list_sci):
            cand1 = 'crf' + self._add_suffix(x, self.bias_suffix)
            cand2 = 'crf' + x
            if (self.data_dir / cand1).exists():
                crf_list.append(cand1)
            elif (self.data_dir / cand2).exists():
                crf_list.append(cand2)
            else:
                continue
        ic = FOSCFileCollection(location=str(self.data_dir), filenames=crf_list)
        extractor = Extract1dSpec(ic_2dspec=ic, disp_axis=self.params.disp_axis)
        extractor.trace_and_extract(guess=guess)
        # Lamps (resolve actual files on disk)
        flamp_list = []
        for x in self._filter_base(self.ic_all.list_cal):
            cand1 = 'f' + self._add_suffix(x, self.bias_suffix)
            cand2 = 'f' + x
            if (self.data_dir / cand1).exists():
                flamp_list.append(cand1)
            elif (self.data_dir / cand2).exists():
                flamp_list.append(cand2)
            else:
                continue
        ic_lamp = FOSCFileCollection(location=str(self.data_dir), filenames=flamp_list)
        extractor.extract_lamp(ic_lamp=ic_lamp)
        # Rename lamp outputs from a<lamp> -> af<lamp> inside data directory
        for lf in flamp_list:
            aout = self.data_dir / ('a' + lf)
            if aout.exists():
                aout.rename(self.data_dir / ('af' + lf[1:]))

    @staticmethod
    def _add_suffix(name: str, suffix: str) -> str:
        p = Path(name)
        return f"{p.stem}{suffix}{p.suffix}"

    def _is_derived_name(self, name: str) -> bool:
        base = Path(name).name
        if Path(base).stem.endswith(self.bias_suffix):
            return True
        # Instrument-renamed base frames can begin with letters that are also
        # used as downstream product prefixes (notably ``cal`` versus ``c``).
        if base.startswith(('bias', 'lampflat', 'cal', 'sci')):
            return False
        if base.startswith(('f', 'cr', 'a', 'af', 'w', 'c')):
            return True
        return False

    def _filter_base(self, names: List[str]) -> List[str]:
        if self.staged:
            return list(names)
        return [n for n in names if not self._is_derived_name(n)]

    def _write_stage_provenance(self, stage, outputs, extra=None, upstream=None):
        if not self.staged:
            return
        parameters = dict(self.params.__dict__)
        parameters.update(extra or {})
        upstream_record = self.provenance.read(upstream) if upstream else None
        upstream_fingerprint = upstream_record.get('fingerprint') if upstream_record else None
        fingerprint = self.provenance.fingerprint(
            stage, self.group_entries, parameters, upstream_fingerprint
        )
        self.provenance.write(stage, fingerprint, outputs, parameters)

    def _stage_is_current(self, stage, outputs, extra=None, upstream=None):
        parameters = dict(self.params.__dict__)
        parameters.update(extra or {})
        upstream_record = self.provenance.read(upstream) if upstream else None
        upstream_fingerprint = upstream_record.get('fingerprint') if upstream_record else None
        fingerprint = self.provenance.fingerprint(
            stage, self.group_entries, parameters, upstream_fingerprint
        )
        return self.provenance.current(stage, fingerprint, outputs)

    def _flat_stage_outputs(self):
        assert self.ic_all is not None
        outputs = [
            self.paths.stage('flat_corrected') / name
            for name in self.ic_all.list_sci + self.ic_all.list_cal
        ]
        outputs.extend([
            self.paths.masters / 'master_flat.fits',
            self.paths.masters / 'normalized_flat.fits',
        ])
        return outputs

    @staticmethod
    def _section_bounds(section: str):
        match = re.fullmatch(
            r'\[\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*\]',
            section,
        )
        if not match:
            raise ValueError(f"Invalid FITS section: {section!r}")
        x1, x2, y1, y2 = map(int, match.groups())
        if x1 > x2 or y1 > y2:
            raise ValueError(f"Reversed FITS section is not supported: {section!r}")
        return x1, x2, y1, y2

    @classmethod
    def _validate_section(cls, section: str, shape, label: str, filename: str):
        x1, x2, y1, y2 = cls._section_bounds(section)
        ny, nx = shape
        if x1 < 1 or y1 < 1 or x2 > nx or y2 > ny:
            raise ValueError(
                f"{label} {section} is outside {filename} with shape {shape}"
            )

    def _validate_input_geometry(self) -> None:
        assert self.ic_all is not None
        groups = {}
        for fname in self._filter_base(self.ic_all.list_allbutbias + self.ic_all.list_bias):
            path = self.data_dir / fname
            shape = (fits.getval(path, 'NAXIS2'), fits.getval(path, 'NAXIS1'))
            groups.setdefault(shape, []).append(fname)
            if self.params.trimsec:
                self._validate_section(self.params.trimsec, shape, 'TRIMSEC', fname)
            if self.params.overscan and self.params.biassec:
                self._validate_section(self.params.biassec, shape, 'BIASSEC', fname)
        if len(groups) > 1:
            detail = '; '.join(
                f"{shape}: {len(files)} file(s)" for shape, files in groups.items()
            )
            raise ValueError(f"Selected calibration inputs have mixed image geometries: {detail}")

    def _report_section_disagreements(self) -> None:
        assert self.ic_all is not None
        names = self._filter_base(self.ic_all.list_allbutbias + self.ic_all.list_bias)
        if not names:
            return
        header = fits.getheader(self.data_dir / names[0])
        for key, configured in (
            ('TRIMSEC', self.params.trimsec),
            ('BIASSEC', self.params.biassec if self.params.overscan else None),
        ):
            raw_value = header.get(key)
            if configured and raw_value and configured != raw_value:
                print(
                    f"Warning: configured {key} {configured} differs from raw FITS "
                    f"header {raw_value}; using configured pyfosc value."
                )

    def _processed_output_is_current(self, path: Path) -> bool:
        header = fits.getheader(path)
        if header.get('PYFOSC') != 'PREWAVECAL':
            return False
        if 'subbias' not in str(header.get('SUBTRACT_BIAS', '')).lower():
            return False
        if self.params.trimsec and header.get('TRIMSEC') != self.params.trimsec:
            return False
        if self.params.trimsec:
            x1, x2, y1, y2 = self._section_bounds(self.params.trimsec)
            if (header.get('NAXIS2'), header.get('NAXIS1')) != (
                y2 - y1 + 1,
                x2 - x1 + 1,
            ):
                return False
        if self.params.overscan and self.params.biassec:
            if header.get('BIASSEC') != self.params.biassec:
                return False
        return True

    @staticmethod
    def _first_present(header, keys: List[str]):
        for k in keys:
            if k in header:
                try:
                    return float(header[k])
                except Exception:
                    try:
                        return float(header.get(k))
                    except Exception:
                        continue
        return None
