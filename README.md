# PyFOSC

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10967240.svg)](https://doi.org/10.5281/zenodo.10967240)

PyFOSC reduces long-slit spectra from FOSC-type instruments, from raw CCD frames to
flux-calibrated, telluric-corrected 1D spectra. It is written for observers who want a
night reduced with two commands and plots to check the result, without installing IRAF.

```bash
cd /path/to/night            # the folder holding the raw FITS files of one night
pyfosc init --telescope XLT --slit slit1.8 --grism G4 --copy-raw
pyfosc run
```

The final spectra are `data/<OBJECT>_<tel><n>_<star>_tel.fits`; the plots to check are in
`QAplots/`.

> PyFOSC is under active development (version 2.0.0-alpha). The Python calibration steps
> have been checked against IRAF reductions on Xinglong 2.16 m nights (grism G4). They have
> not yet been checked on real Lijiang, HCT or P200 data. Always look at the QA plots.

## What it does

- Bias, trimming, flat field, cosmic-ray removal and 1D extraction (optimal and aperture).
- Wavelength calibration by matching each arc lamp to packaged reference lamps, with a
  manual line-marking window as fallback.
- Flux calibration from a standard star, recognised from the `OBJECT` header keyword.
- Telluric correction from the standard star, with a scale measured per band and per night.
- QA plots for every step, and re-runs that skip what is already done.
- No IRAF or PyRAF: it runs on numpy, scipy, astropy, ccdproc and specreduce.

## Supported instruments

| Telescope | Instrument | `--telescope` | Reference lamps packaged for automatic wavelength calibration |
|---|---|---|---|
| Xinglong 2.16 m | BFOSC | `XLT` | Fe-Ar, grism G4, slits 1.8″ and 2.3″ |
| Lijiang 2.4 m | YFOSC | `LJT` | He-Ne, grisms G3, G8, G14, slits 1.8″ and 2.5″ |
| HCT 2 m | HFOSC | `HCT` | Fe-Ar, grism Gr7, slit s167 |
| Palomar 200-inch | DBSP | `P200` | none: mark the lamp lines by hand |

Other grisms and slits work once you have identified a lamp by hand (see
[Calibration steps](docs/calibration-steps.md#lamp-identification)).

## Installation

```bash
git clone https://github.com/rudolffu/pyfosc.git
cd pyfosc
python -m pip install -e .
```

Python 3.9 or later is needed. See the [installation guide](docs/installation.md) for a
clean environment, the optional GUI extras and the legacy IRAF scripts.

## Documentation

| If you want to... | Read |
|---|---|
| reduce your first night, step by step | [Tutorial](docs/tutorial.md) |
| judge whether the reduction is good | [Checking your reduction](docs/checking-your-reduction.md) |
| know what each step does and which options it has | [Calibration steps](docs/calibration-steps.md) |
| find a file or a header keyword | [Files and headers](docs/files.md) |
| fix a run that stopped or looks wrong | [Troubleshooting](docs/troubleshooting.md) |

In a notebook, the same steps are methods of two classes:

```python
from pyfosc.pipeline import PreWaveCal, PostExtraction

pw = PreWaveCal('.')        # bias, flat, cosmic rays, extraction
pe = PostExtraction('.')    # wavelength, flux, telluric
pe.run()
```

`pyfosc notebook --name prewavecal` and `pyfosc notebook --name postcal` write starter
notebooks with every step in its own cell.

## Contributing

- Source code: <https://github.com/rudolffu/pyfosc>
- Bug reports and questions: <https://github.com/rudolffu/pyfosc/issues>
- Run the tests with `pip install pytest` and `pytest -q`. They use synthetic spectra and
  the packaged reference lamps, so no observing data are needed.
- Coding conventions are in [AGENTS.md](AGENTS.md); changes are listed in
  [CHANGES.txt](CHANGES.txt).

## How to cite

Yuming Fu. (2024, April 13). PyFOSC: a pipeline toolbox for BFOSC/YFOSC long-slit
spectroscopy data reduction (Version v1.1.0). Zenodo. <https://doi.org/10.5281/zenodo.10967240>

```bibtex
@software{yuming_fu_2024_10967240,
  author       = {Yuming Fu},
  title        = {{PyFOSC: a pipeline toolbox for BFOSC/YFOSC long-slit spectroscopy data reduction}},
  month        = apr,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v1.1.0},
  doi          = {10.5281/zenodo.10967240},
  url          = {https://doi.org/10.5281/zenodo.10967240}
}
```

Cosmic rays are removed with `ccdproc.cosmicray_lacosmic`, an implementation of L.A.Cosmic;
if you use PyFOSC please also cite van Dokkum (2001, PASP, 113, 1420).

The instruments are described in
[Fan et al. 2016](https://ui.adsabs.harvard.edu/abs/2016PASP..128k5005F/abstract) and
[Zhao et al. 2018](https://ui.adsabs.harvard.edu/abs/2018RAA....18..110Z/abstract)
(Xinglong 2.16 m, BFOSC) and
[Wang et al. 2019](https://ui.adsabs.harvard.edu/abs/2019RAA....19..149W/abstract)
(Lijiang 2.4 m, YFOSC).

## License

PyFOSC is licensed under the BSD 3-Clause License; see [COPYRIGHT](COPYRIGHT).
The standard-star and extinction tables under `iraf_data/onedstds` come from IRAF; their
licence is in [doc/LICENSES](doc/LICENSES).
