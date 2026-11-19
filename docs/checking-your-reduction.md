# Checking your reduction

A pipeline that finishes without errors has not necessarily produced good spectra. This
page lists what to look at after `pyfosc run`, in the order the steps ran, and what each
problem looks like. All plots are in `QAplots/`.

## Before you start

- `pyfosc run` has finished ([tutorial](tutorial.md)).
- Keep the terminal output: warnings printed during the run point at most problems below.

## 1. Extraction: `trace_*.pdf`

One plot per target with the trace of the object along the dispersion axis.

- **Good:** the trace follows the bright ridge of the object smoothly.
- **Bad:** the trace jumps, or follows a neighbouring object or a cosmic ray. Re-extract
  with a starting row: `PreWaveCal('.').extract_1d(guess=<row>)`.

## 2. Wavelength: `wavesol_<lamp>.pdf`

The top panel is the lamp with the identified lines labelled; the bottom panel is the fit
residual of each line. The log line for each lamp looks like

```text
afspecllamp0032_FeAr_bc.fit: 32/45 lines, degree 4, rms 0.359 A (0.13 pix), ...
```

- **Good:** rms of a few tenths of a pixel; residuals scattered around zero with no slope
  or curve; lines spread over the whole wavelength range.
- **Lines found but not used** (32 of 45 above) are normal: at a wide slit, neighbouring
  lines merge and are rejected as outliers.
- **Warning "no lines between X and Y A":** the solution is unconstrained in that range.
  Do not trust wavelengths there; identify more lines by hand (see
  [Calibration steps](calibration-steps.md#lamp-identification)).
- **A faint lamp** gives fewer lines and a worse fit. If you took several lamps, the one
  nearest in time is used for each target; delete or move a bad lamp file out of `data/`
  to exclude it, then re-run with `--overwrite`.

### Check the zero point with a sky line

The lamp and the sky do not illuminate the slit in exactly the same way as a star, and the
instrument flexes during the night. Expect an offset of about a pixel between the lamp
scale and real spectra. Measure it on the [O I] 5577.34 Å sky line in the sky band:

```python
import numpy as np
from pyfosc.specio import read_spec

spec = read_spec('data/wacrfspecltarget0036_J0001+4621_bc.fit')
sky = spec.bands[2]
near = (spec.wave > 5560) & (spec.wave < 5595)
print(spec.wave[near][np.argmax(sky[near])])   # should be close to 5577.3
```

On the Xinglong test night the sky line sat 4 Å (1.4 pixels) from its true wavelength in
every exposure. PyFOSC does not correct this offset; apply it yourself if your science
needs wavelengths better than a pixel.

## 3. Flux calibration: `sensfunc.pdf`

The points are the standard star's bandpasses; the curve is the fitted sensitivity.

- **Good:** a smooth curve with points scattered by a few hundredths of a magnitude
  (the log prints `rms 0.019 mag`). Points inside the telluric bands are not fitted.
- **Bad:** a steep or wavy curve usually means the frame is not the star its header names.
  With two or more standards PyFOSC compares their curves and drops one that disagrees,
  with the warning `labelled <star>, but its sensitivity curve does not match`.
- **Absolute fluxes** depend on slit losses and clouds. Only the shape of the spectrum is
  reliable unless the slit was wide and the night photometric.
- **Red end of blue standards.** With a long-pass order-sorting filter (BFOSC G4 + 385LP),
  second-order blue light of a hot standard reaches the detector redward of about 7600 Å.
  The sensitivity function, and therefore the flux of red targets, is less reliable there.

## 4. Telluric correction

### `telluric_template_<star>.pdf`

The standard star with the continuum used in each band (red), the clean windows that
define it (green bars), and the resulting transmission.

- **Good:** the red line continues the stellar continuum smoothly across each band.
- **Bad:** the red line cuts through a stellar absorption line, or a band is labelled
  `not significant` although it is clearly visible.

### `telluric_<object>_tel.pdf`

Observed and corrected spectrum, the transmission that was applied, and zooms on the O2 B
(6870 Å) and A (7600 Å) bands. The title lists, for each band, the scale applied and where
it came from:

| Origin | Meaning |
|---|---|
| `fit` | measured on this target |
| `night+fit` | this target's measurement combined with the night's value |
| `night` | the night's value (the target is too faint to measure the band) |
| `airmass` | airmass ratio only (nothing better was available) |

- **Good:** the bands are gone and the corrected spectrum joins the continuum on both sides.
- **A narrow spike or dip at the sharp edge of the A band (7594 Å)** is a small
  wavelength offset or resolution difference between standard and target. It is confined
  to one or two pixels; mask them if they matter.
- **An emission line of your target inside a band** (for example Hα of a z ≈ 0.05 galaxy
  in the B band): that band is not fitted for that target and gets the night's value.
  Check the zoom.
- **Warning "the targets need 1.3x the depth seen in <star>"** for an O2 band: expected
  with a very blue standard, whose A band is partly filled in by second-order light. The
  night-level value corrects for it; no action is needed.

`QAplots/telluric_summary.csv` lists every band of every target with its scale, error and
origin. Sort it by `band` to see at a glance whether one target stands out.

## 5. The final spectrum

Plot the final file together with its sky band:

```python
import matplotlib.pyplot as plt
from pyfosc.specio import read_spec

spec = read_spec('data/J0859+4846_XL14_bd284211_tel.fits')
plt.plot(spec.wave, spec.flux)
plt.plot(spec.wave, spec.bands[2] * spec.flux.max() / spec.bands[2].max(), alpha=0.4)
plt.show()
```

Narrow features that line up with strong sky lines (5577, 5890, 6300 Å and the OH lines
redward of 7300 Å) are sky-subtraction residuals, not features of your object.

## Next

- [Troubleshooting](troubleshooting.md) when something above looks wrong
- [Calibration steps](calibration-steps.md) for the options that change each step
