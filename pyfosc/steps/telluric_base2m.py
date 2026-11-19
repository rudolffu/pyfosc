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
                        help="Scale the transmission by the airmass law only")
    parser.add_argument("--min-snr", type=float, default=5.0,
                        help="Significance a fitted band scale needs when the night gives "
                             "no value to fall back on (default 5)")
    parser.add_argument("--o2-exponent", type=float, default=1.0,
                        help="O2 bands scale as airmass_ratio**EXPONENT (default 1)")
    parser.add_argument("--red-limit", type=float, default=None,
                        help="Do not correct bands redward of this wavelength (A)")
    args = parser.parse_args(argv)
    options = dict(fit=not args.no_fit, min_snr=args.min_snr, o2_exponent=args.o2_exponent,
                   red_limit=args.red_limit)
    pe = PostExtraction(args.cwd, overwrite=args.overwrite)
    try:
        try:
            pe.telluric_correct(std_file=args.std_file, star=args.star, **options)
        except StandardNotFound as exc:
            answer = ask_for_standard(pe)
            if answer is None:
                print(f"Warning: {exc}\nTelluric correction was skipped.")
                return 0
            pe.telluric_correct(std_file=answer[0], star=answer[1], **options)
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"Error: {exc}")
        return 1
    print("---DONE---")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
