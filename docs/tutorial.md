# Tutorial: reduce one night

This walk-through takes one night of long-slit data from raw frames to final spectra.
It uses a Xinglong 2.16 m night (BFOSC, grism G4, 1.8″ slit) as the example; for another
telescope only the `pyfosc init` line changes.

## Before you start

- PyFOSC is installed ([installation guide](installation.md)).
- One folder holds the raw FITS files of the night: bias frames, flat fields, arc lamps,
  at least one spectrophotometric standard star, and your targets.
- The standard star's `OBJECT` header keyword is a name PyFOSC knows, for example
  `BD+284211`, `Feige34`, `G191-B2B`, `HZ44`. Spaces, `+`, `-` and letter case do not
  matter. The known stars are the files in `iraf_data/onedstds/*/`.

## 1. Set up the night

```bash
cd /path/to/night
```

```bash
pyfosc init --telescope XLT --slit slit1.8 --grism G4 --copy-raw
```

This writes `myfosc.json` (telescope, slit and grism), makes `raw/` and `data/`, backs the
raw frames up into `raw/`, and puts renamed copies in `data/` together with `imglist.csv`,
a table of what each frame is. Open `imglist.csv` and check that biases, flats, lamps and
targets are labelled correctly: everything later relies on it.

Values for `--telescope` are `XLT`, `LJT`, `HCT` and `P200`. The slit and grism names
(`slit1.8`, `slit2.3`, `slit2.5`, `G3`, `G4`, `G8`, `G14`, `Gr7`, ...) select the reference
lamp for the wavelength calibration.

## 2. Run the pipeline

```bash
pyfosc run
```

The steps run in this order:

| Step | What happens | Output in `data/` |
|---|---|---|
| Bias and trim | master bias subtracted, overscan trimmed | `*_bc` |
| Flat field | master flat built, normalised, applied | `f*` |
| Cosmic rays | L.A.Cosmic cleaning | `crf*` |
| Extraction | trace and extract 1D spectra | `acrf*` (targets), `af*` (lamps) |
| Lamp identification | wavelength solution for each lamp | `database/wavesol_*.json` |
| Dispersion correction | spectra put on a linear wavelength grid | `wacrf*` |
| Flux calibration | sensitivity function from the standard, applied | `sens.fits`, `cwacrf*` |
| Telluric correction | O2 and H2O bands removed | `<OBJECT>_<tel><n>_<star>_tel.fits` |

A run can be repeated: finished steps are skipped. To redo the wavelength, flux and
telluric steps, add `--overwrite`.

## 3. Look at the plots

Open the `QAplots/` folder. At minimum, check:

1. `wavesol_<lamp>.pdf`: residuals of the lamp lines should scatter by a few tenths of a
   pixel with no trend.
2. `sensfunc.pdf`: the sensitivity curve should be smooth.
3. `telluric_<object>_tel.pdf`: the A band near 7600 Å should be gone without leaving a
   spike or a dip.

[Checking your reduction](checking-your-reduction.md) explains each plot and what to do
when one looks wrong. Do this before using the spectra.

## 4. Use the spectra

Each final file holds four spectra on the same wavelength grid:

```python
from pyfosc.specio import read_spec

spec = read_spec('data/J0859+4846_XL14_bd284211_tel.fits')
wave = spec.wave            # Angstrom
flux = spec.flux            # erg / cm2 / s / A, optimal extraction
sky = spec.bands[2]         # sky spectrum that was subtracted
sigma = spec.sigma          # uncertainty (see the note in "Files and headers")
```

The files follow the IRAF multispec convention, so `splot`, `specutils` and most other
tools can read them as well.

## If the automatic choices are not what you want

- **Choose the standard star yourself:**

  ```bash
  pyfosc run --std-file wacrfspecltarget0034_BD+284211_bc.fit --star bd284211 --overwrite
  ```

- **Choose the standard-star table:** add `--caldir subaru` (or `oke1990`, `spec50cal`, ...).
  Tables of the same star differ by a few percent, so use the same one when comparing with
  an older reduction.
- **Run one step at a time**, from the night folder:

  ```bash
  python -m pyfosc.steps.reidentlamp2m
  ```

  ```bash
  python -m pyfosc.steps.wavecal2m
  ```

  ```bash
  python -m pyfosc.steps.telluric_base2m
  ```

  Each accepts `--help`. (`pyfosc run --help` itself lists only the sub-command; its
  options are `--std-file`, `--star`, `--caldir`, `--overwrite`, `--stop-on-error`, `--cwd`.)
- **Work in a notebook:** `pyfosc notebook --name prewavecal` and
  `pyfosc notebook --name postcal` write notebooks with one cell per step.

## Next

- [Checking your reduction](checking-your-reduction.md)
- [Calibration steps](calibration-steps.md) for how each step works and its options
