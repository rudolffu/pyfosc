#!/usr/bin/env python
"""Wavelength and flux calibration of extracted spectra (no IRAF).

``acrf*`` -> ``wacrf*`` (linear wavelength grid) -> ``cwacrf*`` (erg/cm2/s/A),
using the lamp solutions from ``reidentlamp2m``/``identlamp2m`` and a
standard star matched from the OBJECT keyword.
"""
from __future__ import annotations

from pyfosc.steps._postcal_cli import (PostExtraction, StandardNotFound, add_standard_args,
                                       ask_for_standard, base_parser)


def main(argv: list[str] | None = None) -> int:
    parser = base_parser("pyfosc.steps.wavecal2m", __doc__.splitlines()[0])
    add_standard_args(parser)
    parser.add_argument("--order", type=int, default=6,
                        help="Number of spline pieces of the sensitivity fit")
    args = parser.parse_args(argv)
    pe = PostExtraction(args.cwd, overwrite=args.overwrite)
    try:
        print("Dispersion correction...")
        pe.apply_wavelength()
        print("Make sensfunc...")
        try:
            pe.build_sensfunc(std_file=args.std_file, star=args.star, caldir=args.caldir,
                              order=args.order)
        except StandardNotFound as exc:
            answer = ask_for_standard(pe)
            if answer is None:
                print(f"Warning: {exc}\nWavelength calibration is done; "
                      "flux calibration was skipped.")
                return 0
            pe.build_sensfunc(std_file=answer[0], star=answer[1], caldir=args.caldir,
                              order=args.order)
        print("Flux calibration...")
        pe.flux_calibrate()
    except (FileNotFoundError, RuntimeError, ValueError, KeyError) as exc:
        print(f"Error: {exc}")
        return 1
    print("---DONE---")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
