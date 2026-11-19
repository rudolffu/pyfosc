#!/usr/bin/env python
"""Identify arc-lamp lines by hand in a matplotlib window (no IRAF).

The first lamp of each kind is marked interactively; its solution is then
transferred to the remaining lamps. Use ``reidentlamp2m`` for the fully
automatic route based on the packaged templates.
"""
from __future__ import annotations

from pyfosc.steps._postcal_cli import PostExtraction, base_parser


def main(argv: list[str] | None = None) -> int:
    parser = base_parser("pyfosc.steps.identlamp2m", __doc__.splitlines()[0])
    parser.add_argument("--fwidth", type=float, default=10.0, help="Line width in pixels")
    parser.add_argument("--degree", type=int, default=4, help="Degree of the Legendre fit")
    args = parser.parse_args(argv)
    pe = PostExtraction(args.cwd, overwrite=True)
    try:
        pe.identify_lamps(interactive="always", fwidth=args.fwidth, degree=args.degree)
    except (FileNotFoundError, RuntimeError) as exc:
        print(f"Error: {exc}")
        return 1
    print("---DONE---")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
