#!/usr/bin/env bash
# Run the standard reduction from a project directory (the one with myfosc.json).
# All steps are pure Python; IRAF/PyRAF is not required.
set -e
python -m pyfosc.steps.prewavecal
python -m pyfosc.steps.reidentlamp2m
# python -m pyfosc.steps.identlamp2m   # identify lamp lines by hand instead
python -m pyfosc.steps.wavecal2m "$@"
python -m pyfosc.steps.telluric_base2m "$@"
