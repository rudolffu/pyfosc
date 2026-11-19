# Calibration steps

How each step after extraction works, which options it has, and its known limits. Read
this when a default does not suit your data, or when you need to describe the reduction
in a paper. For a first run, use the [tutorial](tutorial.md) instead.

Every step can be run from the night folder as a command or from Python:

| Step | Command | Python (`pe = PostExtraction('.')`) |
|---|---|---|
| Lamp identification | `python -m pyfosc.steps.reidentlamp2m` | `pe.identify_lamps()` |
| Manual identification | `python -m pyfosc.steps.identlamp2m` | `pe.identify_lamps(interactive='always')` |
| Dispersion correction and flux calibration | `python -m pyfosc.steps.wavecal2m` | `pe.apply_wavelength()`, `pe.build_sensfunc()`, `pe.flux_calibrate()` |
| Telluric correction | `python -m pyfosc.steps.telluric_base2m` | `pe.telluric_correct()` |

All commands accept `--cwd <folder>` and `--overwrite`. Without `--overwrite`, outputs
that already exist are kept.

## Lamp identification

**Input** `af*` (extracted lamps). **Output** `data/database/wavesol_<lamp>.json`,
`QAplots/wavesol_<lamp>.pdf`.

1. The lamp is cross-correlated with every packaged reference lamp of the same telescope,
   lamp type and grism (`database/lamp_*.fits`), allowing a shift, a small stretch and a
   flipped dispersion direction. The best match wins.
2. The lines identified in that reference lamp are located in your lamp. Lamp lines are
   images of the slit and flat-topped, so each is centred at the middle of its
   half-maximum points, not with a Gaussian fit.
3. Wavelength against pixel is fitted with a Legendre polynomial (degree 4), clipping
   outliers.
4. **XLT Fe-Ar lamps only:** the solution is then refitted on the lines of the BFOSC lamp
   atlas (`iraf_data/newlinelists/fear_bfosc_atlas.txt`). The packaged reference lamps have
   almost no identified lines between 4900 and 6700 Å, and without this step the solution
   drifts by several Å there.

| Option | Default | Effect |
|---|---|---|
| `--interactive auto\|never\|always` | `auto` | when to open the manual window: only if matching fails, never, or always |
| `--no-refine` | off | keep the reference lamp's own lines (skip step 4) |
| `--refine` | off | also add lines from the He-Ne-Ar list for other lamps; check the plot, misidentifications are possible |
| `--fwidth` | 10 | line width in pixels (about the slit width on the detector) |
| `--degree` | 4 | degree of the polynomial |

### Marking lines by hand

Needed when no reference lamp exists for your grism or slit (all P200 data, LJT Fe-Ar
lamps). `identlamp2m` opens a matplotlib window showing the lamp:

| Key | Action |
|---|---|
| `m` | mark the line under the cursor, then type its wavelength in the terminal |
| `d` | delete the nearest marked line |
| `f` | fit and show the residuals |
| `l` | add further lines from the line list, using the current fit |
| `q` | accept and close |

Mark four or five lines spread over the spectrum, press `f`, then `l`, then `f` again.
The solution of the first lamp is reused as the reference for the other lamps of the same
type in that run. Lamp atlases: the BFOSC and YFOSC calibration-lamp manuals of the
observatories.

**Limits.** The solution describes the lamp. Stars and sky are offset from it by up to a
pixel or two (slit illumination, flexure); see
[Checking your reduction](checking-your-reduction.md#check-the-zero-point-with-a-sky-line).

## Dispersion correction

**Input** `acrf*`, `af*`. **Output** `wacrf*`, `waf*`.

Each spectrum takes the solution of the lamp nearest in time (recorded in the `REFSPEC1`
header keyword) and is rebinned onto a linear wavelength grid with the same number of
pixels. Counts are conserved; the uncertainty band is rebinned in quadrature. All four
bands are kept.

## Flux calibration

**Input** `wacrf*`. **Output** `sens.fits`, `QAplots/sensfunc.pdf`, `cwacrf*`.

1. **Finding the standard.** The `OBJECT` keyword of every spectrum is compared with the
   star names in `iraf_data/onedstds/*/`. With several standards, each one's sensitivity
   curve is computed and a frame whose curve has a different shape from the others is
   dropped (wrong object in the slit). The first remaining one is used.
2. **Sensitivity function.** Counts are summed in the bandpasses of the star's table,
   divided by exposure time, bandwidth and tabulated flux, and corrected for extinction at
   the star's airmass. A spline is fitted through the points, leaving out the telluric
   bands. Magnitudes are converted to flux as IRAF `standard` does (3.68×10⁻²⁰ erg cm⁻²
   s⁻¹ Hz⁻¹ at magnitude 0), so results are directly comparable with IRAF.
3. **Applying it.** Each spectrum is divided by exposure time, dispersion and sensitivity,
   and corrected for extinction at its own airmass. The result is in erg cm⁻² s⁻¹ Å⁻¹.

| Option | Default | Effect |
|---|---|---|
| `--std-file` | automatic | the `wacrf*` file to use as the standard |
| `--star` | from `OBJECT` | star name as in `iraf_data/onedstds` |
| `--caldir` | first match in a fixed order starting with `spec50cal`, `oke1990` | table directory (`subaru`, `irscal`, ...) |
| `--order` | 6 | number of spline pieces of the sensitivity fit |

**Limits.** Airmass is read from the header (`AIRMASS`) or computed from coordinates and
time; if neither works, 1.0 is used with a warning. Extinction is the site's mean curve
(`extinction/baoextinct.dat` for XLT, `LJextinct.dat` for LJT and HCT), not that of your
night. Tables of the same star in different directories differ by a few percent.

## Telluric correction

**Input** `cwacrf*` (or `wacrf*` if no flux calibration was done). **Output**
`<OBJECT>_<tel><n>_<star>_tel.fits`, `telluric_<star>.fits`,
`QAplots/telluric_template_<star>.pdf`, `QAplots/telluric_<object>_tel.pdf`,
`QAplots/telluric_summary.csv`.

The atmosphere's absorption is measured on the standard star and divided out of every
target.

1. **Template.** For each band, a continuum is fitted through clean windows on both sides
   of the band in the standard's spectrum; the transmission is the standard divided by it
   inside the band and exactly 1 elsewhere. A band is kept only if it is detected at more
   than 5σ. Stellar lines of the standard inside a band (Paschen lines of DA white dwarfs,
   He II lines of sdO stars) are estimated from lines of the same series outside the bands
   and removed first.

   | Band | Range (Å) |
   |---|---|
   | O2 B | 6860–6960 |
   | H2O | 7150–7350 |
   | O2 A | 7590–7700 |
   | H2O | 8100–8400 |
   | H2O | 8900–9800 |

2. **Scaling to each target.** A band is applied as `T ** scale`. The starting value is the
   airmass ratio between target and standard. The scale of each band and one wavelength
   shift are then measured on the target by fitting `T ** scale` times a smooth continuum,
   with an error bar.
3. **Night-level values.** The measurements of all targets are averaged per band. Each
   target's own measurement is combined with that average according to their errors: a
   bright target follows its own fit, a faint one gets the night's value.

Why per band and per night: O2 does not change during a night, so all targets should
agree, and they do. But a very blue standard has its A band partly filled in by
second-order light, so targets need about 1.3 times the depth the standard shows, while
the B band needs none. Water does change, which is why it is measured separately.

With two or more usable standards, each target uses the one nearest in time, and the
water bands are interpolated in time between the two standards on either side of it.

| Option | Default | Effect |
|---|---|---|
| `--std-file`, `--star` | automatic | choose the standard |
| `--no-fit` | off | airmass ratio only, no measurement on the targets |
| `--min-snr` | 5 | significance a band needs to be used when the night gives no average |
| `--o2-exponent` | 1 | O2 bands start from `airmass_ratio ** exponent` |
| `--red-limit` | none | do not correct bands redward of this wavelength |

**Limits.**

- The standard must be a hot, nearly featureless star observed with the same setup.
- A resolution difference between standard and target (seeing narrower than the slit)
  leaves a one- or two-pixel residual at the sharp edge of the A band.
- Paschen-line removal has been tested on synthetic spectra only.
- The band at 6270–6330 Å is not corrected (it is weak and overlaps a stellar line of sdO
  standards).

## Related

- [Files and headers](files.md): what each output contains
- [Checking your reduction](checking-your-reduction.md): how to judge the result
