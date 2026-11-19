# Files and headers

Reference for what PyFOSC reads and writes. Paths are relative to the night folder.

## Folder layout

| Path | Content |
|---|---|
| `myfosc.json` | telescope, slit, grism (written by `pyfosc init`) |
| `raw/` | untouched backup of the raw frames |
| `data/` | renamed frames and every intermediate and final product |
| `data/imglist.csv` | table of frames and their types |
| `data/database/` | wavelength solutions (`wavesol_<lamp>.json`) |
| `MasterFrames/` | master bias and master flat |
| `QAplots/` | all diagnostic plots and `telluric_summary.csv` |

## File names

Each processing step adds a prefix to the file name, so a name reads right to left as the
history of the file.

| Name | Step that wrote it | Content |
|---|---|---|
| `<name>_bc` | bias and trim | bias-subtracted, trimmed 2D frame |
| `f<name>_bc` | flat field | flat-fielded 2D frame |
| `crf<name>_bc` | cosmic rays | cleaned 2D frame |
| `acrf<name>_bc` | extraction | 1D spectrum of a target, pixel scale, counts |
| `af<name>_bc` | extraction | 1D spectrum of an arc lamp |
| `wacrf<name>_bc`, `waf<name>_bc` | dispersion correction | on a linear wavelength grid |
| `cwacrf<name>_bc` | flux calibration | erg cm⁻² s⁻¹ Å⁻¹ |
| `<OBJECT>_<tel><n>_<star>_tel.fits` | telluric correction | final spectrum |
| `sens.fits` | flux calibration | sensitivity function |
| `telluric_<star>.fits` | telluric correction | transmission template and its error |

In the final name, `<tel>` is `XL`, `LJ` or `HCT`, `<n>` numbers the spectra of the night,
and `<star>` is the standard used for the telluric correction.

## Spectrum format

Extracted spectra are FITS images of shape `(4, 1, N)`, the IRAF multispec layout:

| Band | Content |
|---|---|
| 1 | optimally extracted spectrum (use this one) |
| 2 | simple aperture sum |
| 3 | sky spectrum that was subtracted |
| 4 | uncertainty |

The wavelength of pixel `i` (starting at 1) is `CRVAL1 + (i - CRPIX1) * CD1_1`.

```python
from pyfosc.specio import read_spec

spec = read_spec('data/cwacrfspecltarget0036_J0001+4621_bc.fit')
spec.wave, spec.flux, spec.sigma   # arrays of length N
spec.bands                         # list of the four bands
spec.header                        # FITS header
```

> **Known problem with band 4.** The uncertainty written by the extraction step is 15 to
> 40 times larger than the pixel-to-pixel scatter of the spectra (measured on a Xinglong
> night). Do not use it for signal-to-noise estimates or as fit weights until this is
> fixed; estimate the noise from the scatter of the spectrum instead.

## Header keywords added by PyFOSC

| Keyword | Written by | Meaning |
|---|---|---|
| `REFSPEC1` | dispersion correction | lamp whose solution was used |
| `WAVERMS` | dispersion correction | rms of the wavelength fit (Å) |
| `WAVENLIN` | dispersion correction | number of lamp lines in the fit |
| `DC-FLAG` | dispersion correction | 0: linear wavelength scale |
| `STDNAME` | flux calibration | standard star |
| `FCAIRM` | flux calibration | airmass used for the extinction correction |
| `BUNIT` | flux calibration | `erg/cm2/s/A` |
| `TELSTD` | telluric correction | standard used for the template |
| `TELSHIFT` | telluric correction | wavelength shift applied to the template (Å) |
| `TELSCO2`, `TELSCH2O` | telluric correction | scale of the best-measured O2 and H2O band |
| `TELEO2`, `TELEH2O` | telluric correction | their errors |
| `TELSRC` | telluric correction | origin of those scales (`fit`, `night+fit`, `night`, `airmass`) |
| `TELBND1`, `TELBND2`, ... | telluric correction | one entry per band: name, scale, error, origin |

`sens.fits` also carries `STDFILE` (the standard's file), `CALDIR` (table directory),
`STDAIRM` (its airmass) and `SENSRMS` (rms of the fit in magnitudes).

## Calibration data shipped with PyFOSC

These live in the PyFOSC folder, not in your night folder.

| Path | Content |
|---|---|
| `database/lamp_*.fits`, `database/idlamp_*` | reference arc lamps and their identified lines |
| `iraf_data/newlinelists/fear_bfosc_atlas.txt` | Fe-Ar lines from the BFOSC lamp atlas |
| `iraf_data/newlinelists/henear_scidoc2215.txt` | He-Ne-Ar line list |
| `iraf_data/onedstds/<dir>/<star>.dat` | standard-star tables (wavelength, AB magnitude, bandwidth) |
| `extinction/baoextinct.dat` | extinction curve, Xinglong |
| `extinction/LJextinct.dat` | extinction curve, Lijiang (also used for HCT) |
| `iraf_data/onedstds/kpnoextinct.dat` | extinction curve used for P200 |

## Related

- [Calibration steps](calibration-steps.md): how the files are produced
- [Checking your reduction](checking-your-reduction.md): the plots in `QAplots/`
