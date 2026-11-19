#!/usr/bin/env python
"""Telluric correction using the standard star (no IRAF).

``cwacrf*`` -> ``<OBJECT>_<tel><n>_<star>_tel.fits``.
"""
from __future__ import annotations

from pyfosc.steps._postcal_cli import (PostExtraction, StandardNotFound, add_standard_args,
                                       ask_for_standard, base_parser)


def main(argv: list[str] | None = None) -> int:
    parser = base_parser("pyfosc.steps.telluric_base2m", __doc__.splitlines()[0])
    add_standard_args(parser)
    parser.add_argument("--no-fit", action="store_true",
                        help="Scale the transmission by the airmass ratio only")
    args = parser.parse_args(argv)
    pe = PostExtraction(args.cwd, overwrite=args.overwrite)
    try:
        try:
            pe.telluric_correct(std_file=args.std_file, star=args.star, fit=not args.no_fit)
        except StandardNotFound as exc:
            answer = ask_for_standard(pe)
            if answer is None:
                print(f"Warning: {exc}\nTelluric correction was skipped.")
                return 0
            pe.telluric_correct(std_file=answer[0], star=answer[1], fit=not args.no_fit)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1
    print("---DONE---")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
