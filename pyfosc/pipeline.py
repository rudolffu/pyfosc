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
from astropy.nddata import CCDData
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
        self.master_bias = mb.build(save=save, method=method)
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

        ic = ImageFileCollection(location=str(self.data_dir), filenames=self.ic_all.list_allbutbias)
        for ccd, fname in ic.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            if params.overscan and params.biassec:
                ccd = ccdp.subtract_overscan(ccd, fits_section=params.biassec, overscan_axis=1)
            if params.trimsec:
                ccd = ccdp.trim_image(ccd, fits_section=params.trimsec)
            ccd = ccdp.subtract_bias(ccd, bias_ref)
            ccd.header['SUBTRACT_BIAS'] = 'subbias'
            if params.trimsec:
                ccd.header['TRIMSEC'] = params.trimsec
            ccd.write(fname, overwrite=True)

    # Master Flat and normalization
    def build_master_flat(self, method: str = 'average', save: bool = True) -> CCDData:
        assert self.ic_all is not None
        mf = MasterFlat(work_dir=str(self.work_dir), list_flat=self.ic_all.list_flat)
        self.master_flat = mf.build(save=save, method=method)
        return self.master_flat

    def normalize_flat(self, degree: int = 25, gauss_sigma: int = 20, save: bool = True) -> CCDData:
        assert self.master_flat is not None
        fnorm = FlatNormalizer(master_flat=self.master_flat, disp_axis=self.params.disp_axis, work_dir=str(self.work_dir))
        fnorm.response_correct(degree=degree)
        fnorm.illumination_correct(gauss_sigma=gauss_sigma, save=save)
        self.normalized_flat = fnorm.normalized_flat
        return self.normalized_flat

    # Apply flat correction
    def apply_flat_correction(self) -> None:
        assert self.ic_all is not None and self.normalized_flat is not None
        for fname in self.ic_all.list_sci + self.ic_all.list_cal:
            ccd = CCDData.read(fname, unit='adu')
            if 'subbias' in ccd.header.get('SUBTRACT_BIAS', '').lower() and 'FLAT_CORRECT' not in ccd.header:
                flat_corrected = ccdp.flat_correct(ccd, self.normalized_flat)
                flat_corrected.header['FLAT_CORRECT'] = 'flatcor'
                flat_corrected.write(f'f{fname}', overwrite=True)

    # Cosmic ray cleaning
    def cosmic_ray_clean(self, sigclip: float = 7.0) -> None:
        assert self.ic_all is not None
        # telescope-based defaults if header is missing
        tel = self.ic_all.telescope
        defaults = {'XLT': (2.2, 7.8), 'LJT': (1.1, 6.6), 'HCT': (1.22, 4.8)}
        d_gain, d_rd = defaults.get(tel, (1.0, 5.0))
        for fname in self.ic_all.list_sci:
            fpath = f'f{fname}'
            ccd = CCDData.read(fpath, unit='adu')
            hdr = ccd.header
            gain = hdr.get('GAIN', d_gain)
            rd = hdr.get('RDNOISE', d_rd)
            new_ccd = ccdp.cosmicray_lacosmic(ccd, readnoise=rd, sigclip=sigclip, verbose=True)
            new_ccd.write('cr' + fpath, overwrite=True)

    # 1D extraction
    def extract_1d(self) -> None:
        assert self.ic_all is not None
        # Objects
        crf_list = [f'crf{x}' for x in self.ic_all.list_sci]
        ic = ImageFileCollection(location=str(self.data_dir), filenames=crf_list)
        extractor = Extract1dSpec(ic_2dspec=ic, disp_axis=self.params.disp_axis)
        extractor.trace_and_extract()
        # Lamps
        flamp_list = [f'f{x}' for x in self.ic_all.list_cal]
        ic_lamp = ImageFileCollection(location=str(self.data_dir), filenames=flamp_list)
        extractor.extract_lamp(ic_lamp=ic_lamp)
        # Rename lamp outputs from a<lamp> -> af<lamp>
        for lf in flamp_list:
            aout = Path('a' + lf)
            if aout.exists():
                aout.rename(Path('af' + lf[1:]))
