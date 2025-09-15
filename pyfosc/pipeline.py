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

import numpy as np
import ccdproc as ccdp
from ccdproc import ImageFileCollection
from astropy.nddata import CCDData, StdDevUncertainty
from astropy import units as u

from .fosc import FOSCFileCollection
from .masters import MasterBias, MasterFlat, FlatNormalizer
from .extraction import Extract1dSpec


@dataclass
class PipelineParams:
    trimsec: Optional[str] = None
    biassec: Optional[str] = None
    disp_axis: int = 0
    overscan: bool = False


class PreWaveCal:
    def __init__(self, work_dir: str | Path = '.'):
        self.work_dir = Path(work_dir)
        # Work with renamed data under ./data by convention
        self.data_dir = self.work_dir / 'data'
        self.settings = self._load_settings(self.work_dir / 'myfosc.json')
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
        ic = FOSCFileCollection(location=str(self.data_dir))
        if ic.summary is None:
            raise FileNotFoundError(f"No FITS files found under {self.data_dir}. Ensure backup step populated ./data.")
        ic.check_groups()
        ic.set_parameters()
        self.ic_all = ic
        p = PipelineParams(
            trimsec=ic.parameters.get('trimsec'),
            biassec=ic.parameters.get('biassec'),
            disp_axis=ic.parameters.get('disp_axis', 0),
            overscan=ic.parameters.get('overscan', False),
        )
        self.params = p

    # Master Bias
    def build_master_bias(self, method: str = 'average', save: bool = True) -> CCDData:
        assert self.ic_all is not None, 'Call discover() first'
        mb = MasterBias(work_dir=str(self.work_dir), list_bias=self.ic_all.list_bias)
        # Build without saving, then store under data/MasterFrames
        self.master_bias = mb.build(save=False, method=method)
        if save:
            outdir = self.work_dir / 'MasterFrames'
            outdir.mkdir(parents=True, exist_ok=True)
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
        for ccd, fname in ic.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            # Idempotency: if already processed (headers present), skip
            if ('SUBTRACT_BIAS' in ccd.header) or ('TRIMSEC' in ccd.header):
                continue
            if params.overscan and params.biassec:
                ccd = ccdp.subtract_overscan(ccd, fits_section=params.biassec, overscan_axis=1)
            if params.trimsec:
                ccd = ccdp.trim_image(ccd, fits_section=params.trimsec)
            ccd = ccdp.subtract_bias(ccd, bias_ref)
            ccd.header['SUBTRACT_BIAS'] = 'subbias'
            if params.trimsec:
                ccd.header['TRIMSEC'] = params.trimsec
            # Save with suffix to keep original untrimmed version
            outname = self._add_suffix(fname, self.bias_suffix)
            outpath = Path(self.data_dir) / outname
            if outpath.exists():
                continue
            ccd.write(outpath, overwrite=True)

    # Master Flat and normalization
    def build_master_flat(self, method: str = 'average', save: bool = True) -> CCDData:
        assert self.ic_all is not None
        # Prefer bias-corrected, trimmed flats if present
        flat_bc = [self._add_suffix(f, self.bias_suffix) for f in self.ic_all.list_flat]
        flat_list = [f for f in flat_bc if (self.data_dir / f).exists()]
        if not flat_list:
            flat_list = self.ic_all.list_flat
        mf = MasterFlat(work_dir=str(self.work_dir), list_flat=flat_list)
        # Build without saving, then store under data/MasterFrames
        self.master_flat = mf.build(save=False, method=method)
        if save:
            outdir = self.work_dir / 'MasterFrames'
            outdir.mkdir(parents=True, exist_ok=True)
            self.master_flat.write(outdir / 'master_flat.fits', overwrite=True)
        return self.master_flat

    def normalize_flat(self, degree: int = 25, gauss_sigma: int = 20, save: bool = True) -> CCDData:
        assert self.master_flat is not None
        # Save normalized outputs under MasterFrames at repo root
        fnorm = FlatNormalizer(master_flat=self.master_flat, disp_axis=self.params.disp_axis, work_dir=str(self.work_dir))
        fnorm.response_correct(degree=degree)
        fnorm.illumination_correct(gauss_sigma=gauss_sigma, save=save)
        self.normalized_flat = fnorm.normalized_flat
        return self.normalized_flat

    # Apply flat correction
    def apply_flat_correction(self) -> None:
        assert self.ic_all is not None and self.normalized_flat is not None
        for fname in self.ic_all.list_sci + self.ic_all.list_cal:
            # Apply flat to bias-corrected, trimmed frames (if present)
            inname = self._add_suffix(fname, self.bias_suffix)
            inpath = Path(self.data_dir) / inname
            if not inpath.exists():
                # Fallback to original
                inpath = Path(self.data_dir) / fname
            ccd = CCDData.read(inpath, unit='adu')
            outpath = Path(self.data_dir) / ('f' + inpath.name)
            if outpath.exists():
                continue
            if 'subbias' in ccd.header.get('SUBTRACT_BIAS', '').lower() and 'FLAT_CORRECT' not in ccd.header:
                flat_corrected = ccdp.flat_correct(ccd, self.normalized_flat)
                flat_corrected.header['FLAT_CORRECT'] = 'flatcor'
                flat_corrected.write(outpath, overwrite=True)

    # Cosmic ray cleaning
    def cosmic_ray_clean(self, sigclip: float = 7.0) -> None:
        assert self.ic_all is not None
        # telescope-based defaults if header is missing
        tel = self.ic_all.telescope
        defaults = {'XLT': (2.2, 7.8), 'LJT': (1.1, 6.6), 'HCT': (1.22, 4.8)}
        d_gain, d_rd = defaults.get(tel, (1.0, 5.0))
        for fname in self._filter_base(self.ic_all.list_sci):
            # work from flat-corrected of bias-corrected trimmed inputs
            inname = 'f' + self._add_suffix(fname, self.bias_suffix)
            fpath = Path(self.data_dir) / inname
            if not fpath.exists():
                fpath = Path(self.data_dir) / ('f' + fname)
            ccd = CCDData.read(fpath, unit='adu')
            hdr = ccd.header
            # Prefer values from headers (case-insensitive, multiple aliases)
            gain = self._first_present(hdr, ['GAIN', 'gain']) or d_gain
            rd = self._first_present(hdr, ['RDNOISE', 'rdnoise', 'RON', 'ron', 'READNOISE', 'readnoise']) or d_rd
            new_ccd = ccdp.cosmicray_lacosmic(ccd, readnoise=rd, sigclip=sigclip, verbose=True)
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
            outpath = Path(self.data_dir) / ('cr' + fpath.name)
            if outpath.exists():
                continue
            new_ccd.write(outpath, overwrite=True)

    # 1D extraction
    def extract_1d(self, guess: float | int | None = None) -> None:
        assert self.ic_all is not None
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
        if base.startswith(('f', 'cr', 'a', 'af', 'w', 'c')):
            return True
        if Path(base).stem.endswith(self.bias_suffix):
            return True
        return False

    def _filter_base(self, names: List[str]) -> List[str]:
        return [n for n in names if not self._is_derived_name(n)]

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
