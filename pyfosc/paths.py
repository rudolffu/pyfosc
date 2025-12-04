"""Helpers for locating pyfosc reference data on disk."""
from __future__ import annotations

import os
from pathlib import Path


def find_data_root(module_path: Path | None = None, cwd: Path | None = None) -> Path:
    """Return the directory that contains the reference ``database`` folder.

    Resolution order:
    1) ``PYFOSC_DATA_DIR`` environment variable (must contain ``database``)
    2) A few locations relative to the calling module (useful when running from source)
    3) The current working directory (useful when the user copies the data locally)
    """
    env_root = os.environ.get("PYFOSC_DATA_DIR")
    if env_root:
        env_path = Path(env_root).expanduser().resolve()
        if (env_path / "database").is_dir():
            return env_path

    if module_path:
        here = Path(module_path).resolve()
        candidates = [here.parent.parent.parent, here.parent.parent]
    else:
        candidates = []

    candidates.append(Path(cwd).resolve() if cwd else Path.cwd())

    for base in candidates:
        if (base / "database").is_dir():
            return base

    raise FileNotFoundError(
        "Could not locate database directory. Set PYFOSC_DATA_DIR to the data root."
    )
