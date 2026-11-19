# Troubleshooting

Problems are grouped by the step where they show up. Each entry gives the message or
symptom, the cause, and what to do. If your problem is not here, open an issue at
<https://github.com/rudolffu/pyfosc/issues> with the terminal output and the relevant plot
from `QAplots/`.

## Wavelength calibration

### "lamp not calibrated, no packaged template for this setup"

There is no reference lamp for your telescope, lamp type and grism (for example Fe-Ar
lamps at Lijiang, or any P200 data). Mark the lines by hand:

```bash
python -m pyfosc.steps.identlamp2m
```

The keys are listed under
[Marking lines by hand](calibration-steps.md#marking-lines-by-hand). One lamp is enough:
its solution is transferred to the other lamps of the same type.

### "lamp not calibrated, template matching failed"

A reference lamp exists but your lamp does not resemble it. Check that `--slit` and
`--grism` in `myfosc.json` describe the observation, and look at the lamp spectrum: a
saturated or nearly empty lamp cannot be matched. Otherwise mark the lines by hand as
above.

### "no lines between X and Y A; the wavelength scale is poorly constrained there"

The fit has no lines over more than a fifth of the detector. Wavelengths in that range can
be wrong by several Å. Add lines there with `identlamp2m` (`l` adds lines from the list
after a first fit).

### Sky lines are not at their wavelengths

An offset of about a pixel that is the same along the spectrum is expected (lamp light
and starlight fill the slit differently). Measure it on [O I] 5577.34 Å as shown in
[Checking your reduction](checking-your-reduction.md#check-the-zero-point-with-a-sky-line)
and apply it yourself. An offset that changes along the spectrum means the solution is
wrong: look at `wavesol_<lamp>.pdf`.

## Flux calibration

### "No wavelength-calibrated spectrum matches a tabulated standard star"

The standard's `OBJECT` keyword is not a known star name. Name the file and the star:

```bash
python -m pyfosc.steps.wavecal2m --std-file wacrfspecltarget0034_BD+284211_bc.fit --star bd284211
```

Star names are the file names in `iraf_data/onedstds/*/` (`bd284211`, `feige34`,
`g191b2b`, `hz44`, ...). Without a standard the spectra are still wavelength calibrated
(`wacrf*`), just not flux calibrated.

### "labelled <star>, but its sensitivity curve does not match the other standard(s)"

One frame named as a standard does not behave like one, usually because another object
was in the slit. It is ignored. If it is the frame you wanted, pass it with `--std-file`
and check `sensfunc.pdf` carefully.

### Fluxes differ from an older IRAF reduction by a few percent

Standard-star tables in different directories differ by that much. Use the same one:

```bash
python -m pyfosc.steps.wavecal2m --caldir subaru --overwrite
```

### "No usable airmass for <object>"

Neither an `AIRMASS` keyword nor usable coordinates and time were found in the header,
so airmass 1.0 was assumed. Add `AIRMASS` to the header of the affected frames and re-run with `--overwrite`.

## Telluric correction

### The A band (7600 Å) is over- or under-corrected

Look at the title of `telluric_<object>_tel.pdf`. If the origin of `O2_7590` is `airmass`,
neither the target nor the night gave a measurement: this happens when the night has
fewer than two targets bright enough to measure the band. There is no better information
in the data; mask the band.

### A spike at 7594 Å or 6867 Å

A one- or two-pixel residual at the sharp edge of the O2 bands, caused by a small
difference in wavelength or resolution between standard and target. Mask those pixels.

### A band was not corrected at all

`telluric_template_<star>.pdf` labels each band. `not significant` means the band was too
weak in the standard to be measured; `no continuum` means the spectrum ends inside the
band or its clean windows; `blue of limit` is the 6270–6330 Å band, which is never
corrected.

### "the targets need 1.3x the depth seen in <star>"

Information, not an error: see
[Checking your reduction](checking-your-reduction.md#4-telluric-correction).

## Running and re-running

### A step was skipped although I changed something

Existing outputs are kept. Add `--overwrite` to the step, or to `pyfosc run`. For the
steps before extraction, delete the output files you want rebuilt.

### `pyfosc run --help` shows no options

The options exist but are not listed there: `--std-file`, `--star`, `--caldir`,
`--overwrite`, `--stop-on-error`, `--cwd`. Each step command
(`python -m pyfosc.steps.wavecal2m --help`) lists its own.

### "pyraf" import errors

Only the legacy scripts under `src/` and `pyfosc.steps.ccdotz`, `makereflat2m`, `doapall`
need PyRAF. Use `pyfosc run` or the step commands named in this documentation instead.

## Related

- [Checking your reduction](checking-your-reduction.md)
- [Calibration steps](calibration-steps.md)
