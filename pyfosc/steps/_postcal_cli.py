"""Shared command-line plumbing for the post-extraction step modules."""
from __future__ import annotations

import argparse
import sys

from pyfosc.postcal import PostExtraction, StandardNotFound


def base_parser(prog: str, description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=prog, description=description)
    parser.add_argument("--cwd", default=".",
                        help="Project directory (with ./data) or the data directory itself")
    parser.add_argument("--overwrite", action="store_true",
                        help="Redo outputs that already exist")
    return parser


def add_standard_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--std-file", help="Wavelength-calibrated standard-star spectrum "
                                           "(wacrf*.fits); default: matched from OBJECT")
    parser.add_argument("--star", help="Standard-star name as in iraf_data/onedstds")
    parser.add_argument("--caldir", help="onedstds directory to take the star from")


def ask_for_standard(pe: PostExtraction) -> tuple[str, str] | None:
    """Prompt for the standard star when none was matched and a terminal is attached."""
    if not sys.stdin or not sys.stdin.isatty():
        return None
    from pyfosc.fits_headers import print_fits_keyword_table

    print_fits_keyword_table(pe.wavecal_science, keyword="OBJECT")
    std_file = input("Enter filename of the standard star spectrum: ").strip()
    star = input("Enter name of the standard star: ").strip()
    return (std_file, star) if std_file and star else None


__all__ = ["PostExtraction", "StandardNotFound", "base_parser", "add_standard_args",
           "ask_for_standard"]
