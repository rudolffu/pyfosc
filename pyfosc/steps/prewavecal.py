#!/usr/bin/env python
import os
from pathlib import Path
from typing import Dict, List
from astropy.nddata import CCDData
from astropy.io import fits

from pyfosc.pipeline import PreWaveCal


def _read_lines(path: Path) -> List[str]:
    if not path.exists():
        return []
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def _params_from_settings(settings: Dict) -> Dict:
    tel = settings['mysettings']['telescope']
    grism = settings['mysettings'].get('Grism', '')
    # Defaults are conservative; override for known setups
    params = {
        'overscan': False,
        'biassec': None,
        'trimsec': None,
        'disp_axis': 0,
    }
    if tel == 'HCT':
        params.update({'trimsec': '[26:250,165:2800]', 'disp_axis': 0, 'overscan': False})
    elif tel == 'LJT':
        # LJT-G3 from notebooks
        params.update({'trimsec': '[201:1400,2301:4130]',
                       'biassec': '[10:40,1:4612]',
                       'disp_axis': 0, 'overscan': True})
    elif tel == 'XLT':
        # XLT-NewG4 from notebooks
        params.update({'trimsec': '[51:1750,681:1350]', 'disp_axis': 1, 'overscan': False})
    return params


def build_lists() -> None:
    # reuse the existing list generator
    import pyfosc.steps.gen_imglists as gen  # noqa: F401 (exec side-effect)


def bias_and_trim(params: Dict) -> CCDData:
    data_dir = Path('.')
    zero_list = _read_lines(data_dir / 'zero.list')
    if not zero_list:
        raise RuntimeError('zero.list is empty. Run pyfosc steps.gen_imglists first or ensure files are present.')
    mb = MasterBias(work_dir='.', list_bias=zero_list)
    master_bias = mb.build(save=True)
    # Overscan and trim master bias if needed
    if params.get('overscan') and params.get('biassec'):
        master_bias = ccdp.subtract_overscan(master_bias, fits_section=params['biassec'], overscan_axis=1)
    if params.get('trimsec'):
        master_bias = ccdp.trim_image(master_bias, fits_section=params['trimsec'])
    # save final master bias
    outdir = Path('MasterFrames'); outdir.mkdir(exist_ok=True)
    master_bias.write(outdir / 'master_bias.fits', overwrite=True)

    flatnall = _read_lines(data_dir / 'flatnall.list')
    if not flatnall:
        return master_bias
    ic = ImageFileCollection(location='.', filenames=flatnall)
    for ccd, fname in ic.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
        if params.get('overscan') and params.get('biassec'):
            ccd = ccdp.subtract_overscan(ccd, fits_section=params['biassec'], overscan_axis=1)
        if params.get('trimsec'):
            ccd = ccdp.trim_image(ccd, fits_section=params['trimsec'])
        ccd = ccdp.subtract_bias(ccd, master_bias)
        ccd.header['SUBTRACT_BIAS'] = 'subbias'
        if params.get('trimsec'):
            ccd.header['TRIMSEC'] = params['trimsec']
        ccd.write(fname, overwrite=True)
    return master_bias


def build_and_normalize_flat(params: Dict) -> CCDData:
    flat_list = _read_lines(Path('flat.list'))
    mf = MasterFlat(work_dir='.', list_flat=flat_list)
    master_flat = mf.build(save=True)
    fnorm = FlatNormalizer(master_flat=master_flat, disp_axis=params['disp_axis'], work_dir='.')
    fnorm.response_correct(degree=25)
    fnorm.illumination_correct(gauss_sigma=20, save=True)
    normalized = fnorm.normalized_flat
    return normalized


def flat_correct(params: Dict, normalized_flat: CCDData) -> None:
    specall = _read_lines(Path('specall.list'))
    for fname in specall:
        ccd = CCDData.read(fname, unit='adu')
        if 'subbias' in ccd.header.get('SUBTRACT_BIAS', '').lower() and 'FLAT_CORRECT' not in ccd.header:
            flat_corrected = ccdp.flat_correct(ccd, normalized_flat)
            flat_corrected.header['FLAT_CORRECT'] = 'flatcor'
            flat_corrected.write(f'f{fname}', overwrite=True)


def cosmic_ray_clean(settings: Dict) -> None:
    # Determine gain/readnoise heuristics (fallback per telescope)
    tel = settings['mysettings']['telescope']
    if tel == 'XLT':
        gain, rdnoise = 2.2, 7.8
    elif tel == 'LJT':
        gain, rdnoise = 1.1, 6.6
    elif tel == 'HCT':
        gain, rdnoise = 1.22, 4.8
    else:
        gain, rdnoise = 1.0, 5.0
    obj_list = _read_lines(Path('objall.list'))
    for fname in obj_list:
        fpath = f'f{fname}'
        ccd = CCDData.read(fpath, unit='adu')
        # override from header if present
        hdr = ccd.header
        gain = hdr.get('GAIN', gain)
        rdnoise = hdr.get('RDNOISE', rdnoise)
        new_ccd = ccdp.cosmicray_lacosmic(ccd, readnoise=rdnoise, sigclip=7, verbose=True)
        new_ccd.write('cr' + fpath, overwrite=True)


def extract_1d(params: Dict) -> None:
    # Objects
    obj_list = _read_lines(Path('objall.list'))
    crf_list = [f'crf{x}' for x in obj_list]
    ic = ImageFileCollection(location='.', filenames=crf_list)
    extractor = Extract1dSpec(ic_2dspec=ic, disp_axis=params['disp_axis'])
    extractor.trace_and_extract()
    # Lamps
    lamp_list = _read_lines(Path('lampall.list'))
    flamp_list = [f'f{x}' for x in lamp_list]
    ic_lamp = ImageFileCollection(location='.', filenames=flamp_list)
    extractor.extract_lamp(ic_lamp=ic_lamp)
    # Rename lamp outputs from a<lamp> -> af<lamp> for downstream compatibility
    for lf in flamp_list:
        aout = Path('a' + lf)
        if aout.exists():
            aout.rename(Path('af' + lf[1:]))


def main(argv: list[str] | None = None) -> int:
    pw = PreWaveCal('.')
    pw.discover()
    pw.build_master_bias()
    pw.calibrate_bias_and_trim()
    pw.build_master_flat()
    pw.normalize_flat()
    pw.apply_flat_correction()
    pw.cosmic_ray_clean()
    pw.extract_1d()
    print('Pre-wavecal reductions complete.')
    return 0


if __name__ == '__main__':  # pragma: no cover
    raise SystemExit(main())
