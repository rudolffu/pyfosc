"""Versioned project paths, manifests, and stage provenance."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import shutil
from typing import Any, Dict, Iterable

from astropy.io import fits


MANIFEST_NAME = "pyfosc-manifest.json"
MANIFEST_VERSION = 1
STAGES = ("calibrated", "flat_corrected", "cosmic_cleaned", "extracted")


def _slug(value: Any) -> str:
    text = str(value or "none").lower().replace(".", "p")
    return re.sub(r"[^a-z0-9]+", "", text) or "none"


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_role(filename: str, header) -> str:
    """Instrument-adapter boundary for canonical input classification."""
    name = Path(filename).name.lower()
    for prefix, role in (
        ("bias", "bias"),
        ("lampflat", "flat"),
        ("cal", "lamp"),
        ("sci", "science"),
        ("object", "science"),
        ("speclflat", "flat"),
        ("specllamp", "lamp"),
        ("specltar", "science"),
        ("speclfluxref", "science"),
    ):
        if name.startswith(prefix):
            return role
    obstype = str(header.get("OBSTYPE", header.get("IMAGETYP", ""))).upper()
    if "BIAS" in obstype:
        return "bias"
    if "FLAT" in obstype:
        return "flat"
    if "LAMP" in obstype or obstype == "CAL":
        return "lamp"
    if "SCI" in obstype or "TARGET" in obstype or "FLUXREF" in obstype:
        return "science"
    return "other"


def group_id(settings: Dict[str, Any], header) -> str:
    mine = settings["mysettings"]
    binning = header.get("CCDSUM") or (
        f"{header.get('XBINNING', 1)}x{header.get('YBINNING', 1)}"
    )
    binning = str(binning).replace(" ", "x")
    return "-".join(
        (
            _slug(mine.get("telescope")),
            _slug(mine.get("Grism")),
            _slug(mine.get("slit")),
            f"{header.get('NAXIS1')}x{header.get('NAXIS2')}",
            _slug(binning),
        )
    )


def build_manifest(
    root: Path,
    settings: Dict[str, Any],
    original_by_canonical: Dict[str, str] | None = None,
) -> Dict[str, Any]:
    input_dir = root / "input"
    entries = []
    original_by_canonical = original_by_canonical or {}
    for path in sorted(input_dir.glob("*.fit*")):
        header = fits.getheader(path)
        stat = path.stat()
        entry = {
            "canonical_name": path.name,
            "original_name": original_by_canonical.get(path.name, path.name),
            "role": classify_role(path.name, header),
            "group": group_id(settings, header),
            "setup": {
                "telescope": settings["mysettings"].get("telescope"),
                "grism": settings["mysettings"].get("Grism"),
                "slit": settings["mysettings"].get("slit"),
                "side": settings["mysettings"].get("side"),
            },
            "naxis1": int(header.get("NAXIS1", 0)),
            "naxis2": int(header.get("NAXIS2", 0)),
            "binning": str(header.get("CCDSUM", "1 1")),
            "date_obs": header.get("DATE-OBS"),
            "object": header.get("OBJECT"),
            "filter": header.get("FILTER"),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": file_sha256(path),
        }
        entries.append(entry)
    groups = {}
    for entry in entries:
        group = groups.setdefault(entry["group"], {role: 0 for role in ("bias", "flat", "lamp", "science", "other")})
        group[entry["role"]] += 1
        group["naxis1"] = entry["naxis1"]
        group["naxis2"] = entry["naxis2"]
        group["binning"] = entry["binning"]
    manifest = {
        "version": MANIFEST_VERSION,
        "input_dir": "input",
        "entries": entries,
        "groups": groups,
    }
    atomic_json_write(root / MANIFEST_NAME, manifest)
    return manifest


def validate_manifest(root: Path, manifest: Dict[str, Any]) -> Dict[str, Any]:
    if manifest.get("version") != MANIFEST_VERSION:
        raise ValueError(f"Unsupported manifest version: {manifest.get('version')}")
    input_dir = root / manifest.get("input_dir", "input")
    listed = {entry["canonical_name"] for entry in manifest["entries"]}
    actual = {path.name for path in input_dir.glob("*.fit*")}
    if actual - listed:
        raise ValueError(f"Unlisted FITS inputs found: {sorted(actual - listed)}")
    if listed - actual:
        raise FileNotFoundError(f"Manifest inputs are missing: {sorted(listed - actual)}")
    changed = False
    for entry in manifest["entries"]:
        path = input_dir / entry["canonical_name"]
        stat = path.stat()
        if stat.st_size != entry["size"] or stat.st_mtime_ns != entry["mtime_ns"]:
            digest = file_sha256(path)
            if digest != entry["sha256"]:
                raise ValueError(f"Immutable input changed: {path}")
            entry["size"] = stat.st_size
            entry["mtime_ns"] = stat.st_mtime_ns
            changed = True
    if changed:
        atomic_json_write(root / MANIFEST_NAME, manifest)
    return manifest


def atomic_json_write(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.name}.tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    os.replace(temp, path)


def atomic_ccd_write(ccd, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.stem}.tmp{path.suffix}")
    ccd.write(temp, overwrite=True)
    os.replace(temp, path)


@dataclass(frozen=True)
class ProjectPaths:
    root: Path
    mode: str
    group: str | None = None

    @classmethod
    def from_settings(cls, root: Path, settings: Dict[str, Any], group=None):
        layout = settings.get("layout") or {"mode": "legacy", "version": 1}
        mode = layout.get("mode", "legacy")
        if mode not in {"legacy", "staged"}:
            raise ValueError(f"Unknown layout mode: {mode}")
        selected = group or layout.get("active_group")
        if mode == "staged" and not selected:
            raise ValueError("A staged project requires layout.active_group or run_group")
        return cls(Path(root), mode, selected)

    @property
    def input(self):
        return self.root / ("input" if self.mode == "staged" else "data")

    @property
    def run_root(self):
        return self.root / "reduced" / str(self.group) if self.mode == "staged" else self.root

    def stage(self, name: str):
        return self.run_root / name if self.mode == "staged" else self.root / "data"

    @property
    def masters(self):
        return self.run_root / "masters" if self.mode == "staged" else self.root / "MasterFrames"

    @property
    def provenance(self):
        return self.run_root / "provenance"

    @property
    def qa(self):
        return self.run_root / "QAplots" if self.mode == "staged" else self.root / "QAplots"

    @property
    def iraf_work(self):
        return self.run_root / "iraf_work"

    def ensure(self):
        if self.mode == "staged":
            for path in [self.masters, self.provenance, self.qa, self.iraf_work / "database"]:
                path.mkdir(parents=True, exist_ok=True)
            for stage in STAGES:
                self.stage(stage).mkdir(parents=True, exist_ok=True)


class ProvenanceStore:
    def __init__(self, paths: ProjectPaths):
        self.paths = paths

    def fingerprint(self, stage: str, inputs: Iterable[Dict], parameters: Dict, upstream=None):
        payload = {
            "stage": stage,
            "inputs": sorted((item["canonical_name"], item["sha256"]) for item in inputs),
            "parameters": parameters,
            "upstream": upstream,
            "software": "pyfosc-2.0.0-alpha",
        }
        return sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()

    def read(self, stage: str):
        path = self.paths.provenance / f"{stage}.json"
        return json.loads(path.read_text()) if path.exists() else None

    def current(self, stage: str, fingerprint: str, outputs: Iterable[Path]):
        record = self.read(stage)
        outputs = list(outputs)
        if not record or record.get("fingerprint") != fingerprint:
            return False
        recorded = record.get("outputs", [])
        if not recorded or isinstance(recorded[0], str):
            return False
        by_path = {item["path"]: item for item in recorded}
        for path in outputs:
            if not path.exists():
                return False
            key = str(path.relative_to(self.paths.root))
            item = by_path.get(key)
            if item is None:
                return False
            stat = path.stat()
            if stat.st_size != item["size"]:
                return False
            if stat.st_mtime_ns != item["mtime_ns"] and file_sha256(path) != item["sha256"]:
                return False
        return True

    def write(self, stage: str, fingerprint: str, outputs: Iterable[Path], parameters: Dict):
        output_records = []
        for path in outputs:
            stat = path.stat()
            output_records.append({
                "path": str(path.relative_to(self.paths.root)),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "sha256": file_sha256(path),
            })
        record = {
            "stage": stage,
            "fingerprint": fingerprint,
            "parameters": parameters,
            "outputs": output_records,
        }
        atomic_json_write(self.paths.provenance / f"{stage}.json", record)


def prepare_iraf_handoff(paths: ProjectPaths, science: Iterable[str], lamps: Iterable[str]):
    if paths.mode != "staged":
        return None
    handoff_manifest = paths.iraf_work / "handoff-manifest.json"
    if handoff_manifest.exists():
        previous = json.loads(handoff_manifest.read_text())
        for item in previous.get("files", []):
            old = paths.iraf_work / item["compatibility"]
            if old.exists():
                old.unlink()
    mapping = []
    config = paths.root / "myfosc.json"
    if config.exists():
        shutil.copy2(config, paths.iraf_work / config.name)
    try:
        from .paths import find_data_root

        reference_database = find_data_root(
            module_path=Path(__file__), cwd=paths.root
        ) / "database"
        target_database = paths.iraf_work / "database"
        if reference_database.is_dir():
            for source in reference_database.rglob("*"):
                if not source.is_file():
                    continue
                target = target_database / source.relative_to(reference_database)
                target.parent.mkdir(parents=True, exist_ok=True)
                if not target.exists():
                    shutil.copy2(source, target)
    except FileNotFoundError:
        pass
    for name in science:
        source = paths.stage("extracted") / name
        compat = f"a{name}"
        shutil.copy2(source, paths.iraf_work / compat)
        mapping.append({"canonical": str(source.relative_to(paths.root)), "compatibility": compat, "role": "science"})
    for name in lamps:
        source = paths.stage("extracted") / name
        compat = f"af{name}"
        shutil.copy2(source, paths.iraf_work / compat)
        mapping.append({"canonical": str(source.relative_to(paths.root)), "compatibility": compat, "role": "lamp"})
    atomic_json_write(handoff_manifest, {"files": mapping})
    return mapping
