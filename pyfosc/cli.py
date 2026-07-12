import argparse
import os
import sys
import json
from pathlib import Path
import subprocess
import shutil

from .project import (
    MANIFEST_NAME,
    ProjectPaths,
    atomic_json_write,
    build_manifest,
    file_sha256,
    validate_manifest,
)


TEMPLATE = None  # no longer used; JSON is built via dict to avoid brace-format issues


def _ensure_dirs(base: Path, staged: bool = True) -> None:
    (base / "raw").mkdir(exist_ok=True)
    (base / ("input" if staged else "data")).mkdir(exist_ok=True)
    if not staged:
        (base / "data" / "database").mkdir(exist_ok=True)


def _write_config(base: Path, telescope: str, slit: str, grism: str, side: str) -> Path:
    cfg_path = base / "myfosc.json"
    existing = json.loads(cfg_path.read_text()) if cfg_path.exists() else None
    cfg = {
        "mysettings": {
            "telescope": telescope,
            "slit": slit,
            "Grism": grism,
            "side": side,
        },
    }
    if existing and "layout" in existing:
        cfg["layout"] = existing["layout"]
    elif existing is None:
        cfg["layout"] = {"version": 2, "mode": "staged", "active_group": None}
    if existing and "prewavecal" in existing:
        cfg["prewavecal"] = existing["prewavecal"]
    atomic_json_write(cfg_path, cfg)
    return cfg_path


def cmd_init(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyfosc init", description="Initialize a working directory for pyfosc")
    parser.add_argument("--cwd", default=os.getcwd(), help="Working directory (default: current directory)")
    parser.add_argument("--telescope", choices=["XLT", "LJT", "HCT", "P200"], help="Instrument/telescope code")
    parser.add_argument("--slit", help="Slit identifier (e.g. slit2.5, s167)")
    parser.add_argument("--grism", help="Grism (e.g. G4, G3, Gr7)")
    parser.add_argument("--side", choices=["blue", "red", ""], default="", help="Channel for P200 (blue/red)")
    parser.add_argument("--copy-raw", action="store_true", help="Copy and backup raw files using pipeline step")
    args = parser.parse_args(argv)

    base = Path(args.cwd).resolve()
    existing_config = base / "myfosc.json"
    staged = True
    if existing_config.exists():
        existing = json.loads(existing_config.read_text())
        staged = existing.get("layout", {}).get("mode") == "staged"
    _ensure_dirs(base, staged=staged)

    tel = args.telescope
    slit = args.slit
    grism = args.grism
    side = args.side

    if tel is None:
        print("Select the instrument:")
        opts = ["BFOSC/XLT (XLT)", "YFOSC/LJT (LJT)", "HFOSC/HCT (HCT)", "DBSP/P200 (P200)"]
        for i, s in enumerate(opts, 1):
            print(f"  {i}. {s}")
        while not tel:
            sel = input("Enter number [1-4]: ").strip()
            tel = {"1": "XLT", "2": "LJT", "3": "HCT", "4": "P200"}.get(sel)

    if tel in {"XLT", "LJT"}:
        slit = slit or input("Enter slit (e.g. slit2.5): ").strip()
        grism = grism or input("Enter Grism (e.g. G4/G3): ").strip()
        side = ""
    elif tel == "HCT":
        slit = slit or input("Enter slit (e.g. s167 or s134): ").strip()
        grism = grism or input("Enter Grism (e.g. Gr7): ").strip()
        side = ""
    elif tel == "P200":
        side = side or ("blue" if input("Channel blue[b]/red[r] [b/r]: ").lower().startswith("b") else "red")
        slit = slit or ""
        grism = grism or ""

    cfg_path = _write_config(base, tel, slit or "", grism or "", side or "")
    print(f"Wrote config: {cfg_path}")

    if args.copy_raw:
        try:
            # Run packaged step that expects myfosc.json in CWD
            subprocess.run([sys.executable, "-m", "pyfosc.steps.list_backup_data"], cwd=str(base), check=True)
            print("Raw data is backed up to ./raw and copied to ./input.")
        except Exception as exc:
            print(f"Warning: copy step failed: {exc}\nYou can run it later: python -m pyfosc.steps.list_backup_data", file=sys.stderr)
    else:
        print("Skipping raw copy. You can run later: python -m pyfosc.steps.list_backup_data")

    return 0


def cmd_run(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyfosc run", description="Run the standard reduction pipeline in the current directory")
    parser.add_argument("--cwd", default=os.getcwd(), help="Working directory (default: current directory)")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop on the first failing step")
    parser.add_argument("--group", help="Manifest run-group id (defaults to layout.active_group)")
    args = parser.parse_args(argv)

    prewave = [sys.executable, "-m", "pyfosc.steps.prewavecal"]
    if args.group:
        prewave.extend(["--group", args.group])
    steps = [
        # New ccdproc/astropy-based pre-wavecal pipeline
        prewave,
        # Existing IRAF-based wave calibration and telluric correction
        "pyfosc.steps.identlamp2m",
        "pyfosc.steps.wavecal2m",
        "pyfosc.steps.telluric_base2m",
    ]

    ok = True
    settings = json.loads((Path(args.cwd) / "myfosc.json").read_text())
    paths = ProjectPaths.from_settings(Path(args.cwd), settings, args.group)
    for step in steps:
        mod = step if isinstance(step, str) else step[2]
        print(f"==> Running {mod}")
        try:
            command = [sys.executable, "-m", step] if isinstance(step, str) else step
            cwd = paths.iraf_work if paths.mode == "staged" and isinstance(step, str) else Path(args.cwd)
            subprocess.run(command, cwd=cwd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Step failed: {mod} (exit {e.returncode})", file=sys.stderr)
            ok = False
            if args.stop_on_error:
                break
    return 0 if ok else 1


def cmd_inspect(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyfosc inspect")
    parser.add_argument("--cwd", default=os.getcwd())
    args = parser.parse_args(argv)
    root = Path(args.cwd).resolve()
    manifest = validate_manifest(root, json.loads((root / MANIFEST_NAME).read_text()))
    active = json.loads((root / "myfosc.json").read_text()).get("layout", {}).get("active_group")
    for name, group in sorted(manifest["groups"].items()):
        marker = "*" if name == active else " "
        counts = " ".join(f"{role}={group.get(role, 0)}" for role in ("bias", "flat", "lamp", "science"))
        print(f"{marker} {name}: {counts}")
    return 0


def cmd_migrate_layout(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyfosc migrate-layout")
    parser.add_argument("--cwd", default=os.getcwd())
    parser.add_argument("--apply", action="store_true", help="Apply the migration; default is dry-run")
    args = parser.parse_args(argv)
    root = Path(args.cwd).resolve()
    settings = json.loads((root / "myfosc.json").read_text())
    source = root / "data"
    candidates = []
    ambiguous = []
    for path in sorted(source.glob("*.fit*")):
        stem = path.stem
        if stem.endswith("_bc") or stem.startswith(
            ("f", "cr", "a", "w", "cwa", "calcon")
        ):
            continue
        if stem.startswith(("bias", "lampflat", "cal", "sci", "object", "specl")):
            candidates.append(path)
        else:
            ambiguous.append(path)
    print(f"Migration {'plan' if not args.apply else 'applying'}: {len(candidates)} base input(s)")
    for path in candidates:
        target = root / "input" / path.name
        collision = " [collision]" if target.exists() else ""
        print(f"  {path.relative_to(root)} -> input/{path.name}{collision}")
    for path in ambiguous:
        print(f"  ambiguous (left untouched): {path.relative_to(root)}")
    if not args.apply:
        return 0
    input_dir = root / "input"
    input_dir.mkdir(exist_ok=True)
    for path in candidates:
        target = input_dir / path.name
        if target.exists():
            same = target.stat().st_size == path.stat().st_size
            same = same and file_sha256(target) == file_sha256(path)
            if not same:
                raise FileExistsError(f"Migration collision: {target}")
        if not target.exists():
            shutil.copy2(path, target)
    manifest = build_manifest(root, settings)
    groups = [name for name, data in manifest["groups"].items() if data.get("science")]
    settings["layout"] = {
        "version": 2,
        "mode": "staged",
        "active_group": groups[0] if len(groups) == 1 else None,
    }
    atomic_json_write(root / "myfosc.json", settings)
    print("Migration complete; existing data/ products were left untouched.")
    return 0


def _write_ipynb(path: Path, cells: list[dict]) -> None:
    nb = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"}
        },
        "nbformat": 4,
        "nbformat_minor": 5
    }
    path.write_text(json.dumps(nb, indent=2))


def cmd_notebook(name: str = "prewavecal", cwd: str | None = None, open_nb: bool = False) -> int:
    base = Path(cwd or os.getcwd()).resolve()
    base.mkdir(parents=True, exist_ok=True)
    if name == "prewavecal":
        nb_path = base / "pyfosc_prewavecal.ipynb"
        cells = [
            {"cell_type": "markdown", "metadata": {}, "source": [
                "# pyfosc pre-wavecal reductions\n",
                "This notebook initializes lists, builds master frames, flat-corrects, runs cosmic-ray cleaning, and extracts 1D spectra.\n"
            ]},
            {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": [
                "from pyfosc.pipeline import PreWaveCal\n",
                "pw = PreWaveCal('.')\n",
                "pw.discover()\n",
                "pw.build_master_bias()\n",
                "pw.calibrate_bias_and_trim()\n",
                "pw.build_master_flat()\n",
                "pw.normalize_flat()\n",
                "pw.apply_flat_correction()\n",
                "pw.cosmic_ray_clean()\n",
                "pw.extract_1d()\n"
            ]}
        ]
        _write_ipynb(nb_path, cells)
        print(f"Wrote {nb_path}")
        if open_nb:
            try:
                subprocess.run([sys.executable, "-m", "jupyter", "notebook", str(nb_path)], check=False)
            except Exception as exc:
                print(f"Couldn't open Jupyter automatically: {exc}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyfosc", description="pyfosc command line interface")
    subparsers = parser.add_subparsers(dest="cmd")

    subparsers.add_parser("init", help="Initialize a working directory for pyfosc")
    subparsers.add_parser("run", help="Run the standard reduction pipeline")
    subparsers.add_parser("inspect", help="List manifest run groups and calibration coverage")
    subparsers.add_parser("migrate-layout", help="Dry-run or apply legacy layout migration")
    nb = subparsers.add_parser("notebook", help="Generate a runnable Jupyter notebook for a step")
    nb.add_argument("--name", choices=["prewavecal"], default="prewavecal", help="Notebook template to generate")
    nb.add_argument("--cwd", default=os.getcwd(), help="Directory where to create the notebook")
    nb.add_argument("--open", action="store_true", dest="open_nb", help="Open notebook with Jupyter after creation")

    args, rest = parser.parse_known_args(argv)
    if args.cmd == "init":
        return cmd_init(rest)
    elif args.cmd == "run":
        return cmd_run(rest)
    elif args.cmd == "notebook":
        return cmd_notebook(name=args.name, cwd=args.cwd, open_nb=args.open_nb)
    elif args.cmd == "inspect":
        return cmd_inspect(rest)
    elif args.cmd == "migrate-layout":
        return cmd_migrate_layout(rest)
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
