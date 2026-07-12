import json
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from pyfosc.pipeline import PreWaveCal
from pyfosc.cli import _write_config, cmd_migrate_layout
from pyfosc.project import (
    MANIFEST_NAME,
    ProjectPaths,
    ProvenanceStore,
    build_manifest,
    prepare_iraf_handoff,
    validate_manifest,
)


def _write(path, role, value=10.0):
    prefixes = {
        "bias": "bias",
        "flat": "lampflat",
        "lamp": "cal",
        "science": "sci",
    }
    filename = f"{prefixes[role]}001.fits"
    ccd = CCDData(np.full((10, 12), value), unit=u.adu)
    ccd.header["TELESCOP"] = "Lijiang 2m4"
    ccd.header["FILTER"] = "lslit2_51*grism3"
    ccd.header["FILTER1"] = "lslit2_51"
    ccd.header["FILTER3"] = "grism3"
    ccd.header["CCDSUM"] = "1 1"
    ccd.write(path / filename)
    return filename


def _project(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    for role, value in (("bias", 1), ("flat", 20), ("lamp", 11), ("science", 12)):
        _write(input_dir, role, value)
    settings = {
        "mysettings": {
            "telescope": "LJT",
            "slit": "slit2.5",
            "Grism": "G3",
            "side": "",
        },
        "prewavecal": {
            "trimsec": "[3:8,2:9]",
            "overscan": False,
            "biassec": None,
            "disp_axis": 0,
        },
        "layout": {"version": 2, "mode": "staged", "active_group": None},
    }
    manifest = build_manifest(tmp_path, settings)
    group = next(iter(manifest["groups"]))
    settings["layout"]["active_group"] = group
    (tmp_path / "myfosc.json").write_text(json.dumps(settings))
    return group


def test_manifest_is_authoritative_and_reduced_files_do_not_affect_it(tmp_path):
    _project(tmp_path)
    reduced = tmp_path / "reduced" / "anything" / "calibrated"
    reduced.mkdir(parents=True)
    (reduced / "sci999.fits").write_bytes(b"derived")

    manifest = validate_manifest(
        tmp_path, json.loads((tmp_path / MANIFEST_NAME).read_text())
    )

    assert len(manifest["entries"]) == 4


def test_manifest_rejects_unlisted_input(tmp_path):
    _project(tmp_path)
    (tmp_path / "input" / "extra.fits").write_bytes(b"not listed")

    with pytest.raises(ValueError, match="Unlisted FITS inputs"):
        validate_manifest(tmp_path, json.loads((tmp_path / MANIFEST_NAME).read_text()))


def test_staged_discovery_and_calibration_use_stable_names(tmp_path):
    group = _project(tmp_path)
    pipeline = PreWaveCal(tmp_path, run_group=group)
    pipeline.discover()
    pipeline.build_master_bias()
    pipeline.calibrate_bias_and_trim()

    assert pipeline.ic_all.list_bias == ["bias001.fits"]
    assert (tmp_path / "reduced" / group / "calibrated" / "sci001.fits").exists()
    assert not (tmp_path / "input" / "sci001_bc.fits").exists()
    assert (tmp_path / "reduced" / group / "provenance" / "calibrated.json").exists()


def test_changed_geometry_rebuilds_staged_calibrated_outputs(tmp_path):
    group = _project(tmp_path)
    pipeline = PreWaveCal(tmp_path, run_group=group)
    pipeline.discover()
    pipeline.build_master_bias()
    pipeline.calibrate_bias_and_trim()
    output = tmp_path / "reduced" / group / "calibrated" / "sci001.fits"
    assert CCDData.read(output, unit=u.adu).shape == (8, 6)

    settings = json.loads((tmp_path / "myfosc.json").read_text())
    settings["prewavecal"]["trimsec"] = "[4:7,3:8]"
    (tmp_path / "myfosc.json").write_text(json.dumps(settings))
    changed = PreWaveCal(tmp_path, run_group=group)
    changed.discover()
    changed.build_master_bias()
    changed.calibrate_bias_and_trim()

    assert CCDData.read(output, unit=u.adu).shape == (6, 4)


def test_flat_corrected_stage_uses_separate_directory(tmp_path):
    group = _project(tmp_path)
    pipeline = PreWaveCal(tmp_path, run_group=group)
    pipeline.discover()
    pipeline.build_master_bias()
    pipeline.calibrate_bias_and_trim()
    pipeline.build_master_flat()
    pipeline.normalize_flat(degree=2, gauss_sigma=1)
    pipeline.apply_flat_correction()

    stage = tmp_path / "reduced" / group / "flat_corrected"
    assert (stage / "sci001.fits").exists()
    assert (stage / "cal001.fits").exists()
    assert not (tmp_path / "input" / "fsci001.fits").exists()


def test_provenance_fingerprint_changes_with_parameters(tmp_path):
    group = _project(tmp_path)
    settings = json.loads((tmp_path / "myfosc.json").read_text())
    paths = ProjectPaths.from_settings(tmp_path, settings, group)
    store = ProvenanceStore(paths)
    manifest = json.loads((tmp_path / MANIFEST_NAME).read_text())

    first = store.fingerprint("calibrated", manifest["entries"], {"trimsec": "a"})
    second = store.fingerprint("calibrated", manifest["entries"], {"trimsec": "b"})

    assert first != second


def test_legacy_migration_defaults_to_dry_run(tmp_path, capsys):
    data = tmp_path / "data"
    data.mkdir()
    _write(data, "science")
    settings = {
        "mysettings": {
            "telescope": "LJT", "slit": "slit2.5", "Grism": "G3", "side": ""
        }
    }
    (tmp_path / "myfosc.json").write_text(json.dumps(settings))

    assert cmd_migrate_layout(["--cwd", str(tmp_path)]) == 0

    assert not (tmp_path / "input").exists()
    assert "data/sci001.fits -> input/sci001.fits" in capsys.readouterr().out


def test_iraf_handoff_uses_compatibility_copies(tmp_path, monkeypatch):
    group = _project(tmp_path)
    settings = json.loads((tmp_path / "myfosc.json").read_text())
    paths = ProjectPaths.from_settings(tmp_path, settings, group)
    paths.ensure()
    (paths.stage("extracted") / "sci001.fits").write_bytes(b"science")
    (paths.stage("extracted") / "cal001.fits").write_bytes(b"lamp")
    monkeypatch.setattr("pyfosc.paths.find_data_root", lambda **kwargs: tmp_path)

    mapping = prepare_iraf_handoff(paths, ["sci001.fits"], ["cal001.fits"])

    assert (paths.iraf_work / "asci001.fits").read_bytes() == b"science"
    assert (paths.iraf_work / "afcal001.fits").read_bytes() == b"lamp"
    assert len(mapping) == 2
    assert (paths.stage("extracted") / "sci001.fits").read_bytes() == b"science"


def test_reinitializing_legacy_config_does_not_enable_staged_layout(tmp_path):
    legacy = {
        "mysettings": {
            "telescope": "LJT", "slit": "old", "Grism": "G3", "side": ""
        },
        "prewavecal": {"trimsec": "[1:2,1:2]"},
    }
    (tmp_path / "myfosc.json").write_text(json.dumps(legacy))

    _write_config(tmp_path, "LJT", "slit2.5", "G3", "")

    updated = json.loads((tmp_path / "myfosc.json").read_text())
    assert "layout" not in updated
    assert updated["prewavecal"] == legacy["prewavecal"]
