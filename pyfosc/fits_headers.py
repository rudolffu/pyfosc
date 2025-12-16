from __future__ import annotations

from pathlib import Path
from typing import Iterable

from astropy.io import fits


def print_fits_keyword_table(files: Iterable[str], keyword: str = "OBJECT") -> None:
    """Print a simple two-column table: filename + FITS header keyword value."""
    filenames = [str(f) for f in files]
    if not filenames:
        print(f'No files matched for keyword listing ("{keyword}").')
        return

    rows: list[tuple[str, str]] = []
    for filename in sorted(filenames):
        display_name = Path(filename).name
        value = ""
        try:
            header = fits.getheader(filename, 0)
            value = header.get(keyword, "")
        except Exception as exc:
            value = f"<error: {exc.__class__.__name__}>"
        rows.append((display_name, str(value).strip()))

    name_width = max(len("Filename"), *(len(r[0]) for r in rows))
    value_width = max(len(keyword), *(len(r[1]) for r in rows))

    print(f'{"Filename":<{name_width}}  {keyword:<{value_width}}')
    print(f'{"-" * name_width}  {"-" * value_width}')
    for display_name, value in rows:
        print(f"{display_name:<{name_width}}  {value:<{value_width}}")

