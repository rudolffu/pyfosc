import argparse
import os
import sys
import json
from pathlib import Path
import subprocess


TEMPLATE = None  # no longer used; JSON is built via dict to avoid brace-format issues


def _ensure_dirs(base: Path) -> None:
    (base / "raw").mkdir(exist_ok=True)
    (base / "data").mkdir(exist_ok=True)
    (base / "data" / "database").mkdir(exist_ok=True)


def _write_config(base: Path, telescope: str, slit: str, grism: str, side: str) -> Path:
    cfg = {
        "mysettings": {
            "telescope": telescope,
            "slit": slit,
            "Grism": grism,
            "side": side,
        }
    }
    cfg_text = json.dumps(cfg, indent=2)
    cfg_path = base / "myfosc.json"
    cfg_path.write_text(cfg_text)
    # also copy to data/
    data_dir = base / "data"
    data_dir.mkdir(exist_ok=True)
    (data_dir / "myfosc.json").write_text(cfg_text)
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
    _ensure_dirs(base)

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
            print("Raw data is backed up to ./raw and copied to ./data.")
        except Exception as exc:
            print(f"Warning: copy step failed: {exc}\nYou can run it later: python -m pyfosc.steps.list_backup_data", file=sys.stderr)
    else:
        print("Skipping raw copy. You can run later: python -m pyfosc.steps.list_backup_data")

    return 0


def cmd_run(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pyfosc run", description="Run the standard reduction pipeline in the current directory")
    parser.add_argument("--cwd", default=os.getcwd(), help="Working directory (default: current directory)")
    parser.add_argument("--stop-on-error", action="store_true", help="Stop on the first failing step")
    args = parser.parse_args(argv)

    steps = [
        # New ccdproc/astropy-based pre-wavecal pipeline
        "pyfosc.steps.prewavecal",
        # Existing IRAF-based wave calibration and telluric correction
        "pyfosc.steps.identlamp2m",
        "pyfosc.steps.wavecal2m",
        "pyfosc.steps.telluric_base2m",
    ]

    ok = True
    for mod in steps:
        print(f"==> Running {mod}")
        try:
            subprocess.run([sys.executable, "-m", mod], cwd=args.cwd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Step failed: {mod} (exit {e.returncode})", file=sys.stderr)
            ok = False
            if args.stop_on_error:
                break
    return 0 if ok else 1


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
    else:
        parser.print_help()
        return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
