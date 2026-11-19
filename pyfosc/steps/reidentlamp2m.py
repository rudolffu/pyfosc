#!/usr/bin/env python
"""Wavelength solutions for arc lamps from the packaged templates (no IRAF).

Each ``af*.fits`` lamp is matched to the reference lamps under ``database/``;
solutions are saved as ``database/wavesol_<lamp>.json`` with QA plots.
"""
from __future__ import annotations

from pyfosc.steps._postcal_cli import PostExtraction, base_parser


def main(argv: list[str] | None = None) -> int:
    parser = base_parser("pyfosc.steps.reidentlamp2m", __doc__.splitlines()[0])
    parser.add_argument("--interactive", choices=["auto", "never", "always"], default="auto",
                        help="When to open the matplotlib identifier (default: only if "
                             "template matching fails and a terminal is attached)")
    parser.add_argument("--refine", action="store_true", default=None,
                        help="Refit on the line list after template matching (default: "
                             "only for XLT Fe-Ar lamps, using the BFOSC atlas lines)")
    parser.add_argument("--no-refine", dest="refine", action="store_false",
                        help="Keep the lines of the matched template only")
    parser.add_argument("--fwidth", type=float, default=10.0, help="Line width in pixels")
    parser.add_argument("--degree", type=int, default=4, help="Degree of the Legendre fit")
    args = parser.parse_args(argv)
    pe = PostExtraction(args.cwd, overwrite=args.overwrite)
    try:
        pe.identify_lamps(interactive=args.interactive, fwidth=args.fwidth, degree=args.degree,
                          refine=args.refine)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"Error: {exc}")
        return 1
    print("---DONE---")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
