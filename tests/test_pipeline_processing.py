from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from pyfosc.pipeline import PipelineParams, PreWaveCal
from pyfosc.masters import MasterBias
from pyfosc.fosc import FOSCFileCollection


def _write_ccd(path, value=10.0, shape=(10, 12), **headers):
    ccd = CCDData(np.full(shape, value, dtype=float), unit=u.adu)
    for key, value in headers.items():
        ccd.header[key] = value
    ccd.write(path)


def _pipeline(tmp_path, trimsec="[3:8,2:9]"):
    (tmp_path / "data").mkdir()
    (tmp_path / "myfosc.json").write_text('{"mysettings": {}}')
    pipeline = PreWaveCal(tmp_path)
    pipeline.params = PipelineParams(trimsec=trimsec, overscan=False)
    pipeline.ic_all = SimpleNamespace(
        list_allbutbias=["sci001.fits"],
        list_bias=["bias001.fits"],
        list_flat=[],
        list_sci=["sci001.fits"],
        list_cal=[],
    )
    pipeline.master_bias = CCDData(np.ones((10, 12)), unit=u.adu)
    return pipeline


def test_raw_trimsec_header_does_not_skip_processing(tmp_path):
    pipeline = _pipeline(tmp_path)
    _write_ccd(
        tmp_path / "data" / "sci001.fits",
        TRIMSEC="[1:12,1:10]",
    )

    pipeline.calibrate_bias_and_trim()

    output = tmp_path / "data" / "sci001_bc.fits"
    assert output.exists()
    assert fits.getdata(output).shape == (8, 6)
    header = fits.getheader(output)
    assert header["PYFOSC"] == "PREWAVECAL"
    assert header["TRIMSEC"] == "[3:8,2:9]"
    assert header["SUBTRACT_BIAS"] == "subbias"


def test_calibration_lamp_name_is_not_treated_as_derived(tmp_path):
    pipeline = _pipeline(tmp_path)

    assert pipeline._filter_base(["cal0151_G191-B2B.fits"]) == [
        "cal0151_G191-B2B.fits"
    ]
    assert pipeline._filter_base(
        ["ccal0151_G191-B2B.fits", "cal0151_G191-B2B_bc.fits"]
    ) == []


def test_changed_trimsec_regenerates_stale_output(tmp_path):
    pipeline = _pipeline(tmp_path)
    _write_ccd(tmp_path / "data" / "sci001.fits", TRIMSEC="[1:12,1:10]")
    pipeline.calibrate_bias_and_trim()

    pipeline.params.trimsec = "[4:7,3:8]"
    pipeline.calibrate_bias_and_trim()

    output = tmp_path / "data" / "sci001_bc.fits"
    assert fits.getdata(output).shape == (6, 4)
    assert fits.getheader(output)["TRIMSEC"] == "[4:7,3:8]"


def test_geometry_validation_rejects_mixed_selected_inputs(tmp_path):
    pipeline = _pipeline(tmp_path, trimsec="[1:5,1:5]")
    pipeline.ic_all.list_allbutbias = ["sci001.fits"]
    _write_ccd(tmp_path / "data" / "sci001.fits", shape=(10, 12))
    _write_ccd(tmp_path / "data" / "bias001.fits", shape=(8, 12))

    with pytest.raises(ValueError, match="mixed image geometries"):
        pipeline._validate_input_geometry()


def test_master_flat_requires_processed_flats(tmp_path):
    pipeline = _pipeline(tmp_path)
    pipeline.ic_all.list_flat = ["lampflat001.fits"]

    with pytest.raises(RuntimeError, match="Processed flat inputs are missing"):
        pipeline.build_master_flat()


def test_invalid_section_reports_filename(tmp_path):
    pipeline = _pipeline(tmp_path, trimsec="[1:20,1:5]")
    _write_ccd(tmp_path / "data" / "sci001.fits", shape=(10, 12))
    _write_ccd(tmp_path / "data" / "bias001.fits", shape=(10, 12))

    with pytest.raises(ValueError, match="sci001.fits"):
        pipeline._validate_input_geometry()


def test_empty_master_bias_list_does_not_scan_data_directory(tmp_path):
    (tmp_path / "data").mkdir()
    _write_ccd(tmp_path / "data" / "unrelated.fits")
    builder = MasterBias(work_dir=tmp_path, list_bias=[])

    assert builder.build() is None


def test_ljt_derived_only_collection_can_be_used_for_extraction(tmp_path):
    filename = "crfsci001_target_bc.fits"
    _write_ccd(
        tmp_path / filename,
        TELESCOP="Lijiang 2m4",
        FILTER="lslit2_51*grism3",
    )
    collection = FOSCFileCollection(location=tmp_path, filenames=[filename])

    collection.check_groups()

    assert collection.telescope == "LJT"
    assert collection.list_sci == [filename]
    assert collection.list_bias == []
