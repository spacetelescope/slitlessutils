import numpy as np
from pathlib import Path
import pytest
from tempfile import TemporaryDirectory

from astroquery.mast import Observations

from slitlessutils.core.preprocess.crrej.drizzle import (drizzle, group_by_visit,
                                                         group_by_position_angle)
from slitlessutils.core.wfss import WFSSCollection


@pytest.fixture
def mock_files():
    return [
        "file1_visit_01.fits",
        "file2_visit_01.fits",
        "file3_visit_02.fits",
        "file4_visit_03.fits",
        "file5_visit_02.fits",
    ]


@pytest.fixture
def mock_visits():
    return ["01", "01", "02", "03", "02"]


@pytest.fixture
def mock_pas():
    return [42.13, 42.14, 80.0, 179.99, -179.98, 360.0, 0.01, 0.05, 0.09]


def test_wr96_drizzle():
    """
    Runs the drizzle step embedded in the wr96 example
    """
    mast_links = {
        "jdql01jpq_flc.fits": "mast:HST/product/jdql01jpq_flc.fits",
        "jdql01jxq_flc.fits": "mast:HST/product/jdql01jxq_flc.fits",
    }
    # Setup temporary directory
    with TemporaryDirectory() as tempdir:
        rawdata_dir = Path(tempdir) / "raw_data"
        rawdata_dir.mkdir(parents=True)
        mosaic_dir = Path(tempdir) / "mosaic_files"
        mosaic_dir.mkdir(parents=True)
        # Download wr96 visit files
        for filename in mast_links:
            Observations.download_file(
                mast_links[filename], local_path=str(rawdata_dir / filename)
            )
        rawdata_filepaths = [str(filepath) for filepath in rawdata_dir.iterdir()]

        # Check that our temp folder is indeed empty
        assert len(list(mosaic_dir.iterdir())) == 0
        # Actually perform drizzle
        drizzle(rawdata_filepaths, outdir=mosaic_dir)
        # Confirm we have our output mosaics
        assert len(list(mosaic_dir.iterdir())) > 0


@pytest.mark.parametrize("return_unique_visits", [False, True])
def test_group_by_visits(mock_files, mock_visits, return_unique_visits, monkeypatch):
    def mock_from_list(files):
        return WFSSCollection()

    def mock_get_visits(self):
        return mock_visits

    monkeypatch.setattr(WFSSCollection, "from_list", mock_from_list)
    monkeypatch.setattr(WFSSCollection, "get_visits", mock_get_visits)

    result = group_by_visit(mock_files, return_unique_visits=return_unique_visits)

    if return_unique_visits:
        grouped_files = result[0]
        unique_visits = result[1]
    else:
        grouped_files = result

    assert len(grouped_files) == 3
    assert len(grouped_files[0]) == 2
    assert len(grouped_files[1]) == 2
    assert len(grouped_files[2]) == 1

    assert grouped_files[0] == ["file1_visit_01.fits", "file2_visit_01.fits"]
    assert grouped_files[1] == ["file3_visit_02.fits", "file5_visit_02.fits"]
    assert grouped_files[2] == ["file4_visit_03.fits"]

    if return_unique_visits:
        assert len(unique_visits) == 3
        assert np.array_equal(unique_visits, ["01", "02", "03"])


def test_group_by_pas(mock_files, mock_pas, monkeypatch):
    np.random.seed(42)

    def mock_from_list(files):
        return WFSSCollection()

    def mock_get_pas(self):
        return mock_pas

    monkeypatch.setattr(WFSSCollection, "from_list", mock_from_list)
    monkeypatch.setattr(WFSSCollection, "get_pas", mock_get_pas)

    additional_files = ["add_test1", "add_test2", "add_test3", "add_test4"]
    result = group_by_position_angle(mock_files + additional_files)
    print(result)

    assert len(result) == 5
    assert result[0] == ["file4_visit_03.fits", "file5_visit_02.fits"]
    assert result[1] == ["add_test1", "add_test2"]
    assert result[2] == ['add_test3', 'add_test4']
    assert result[3] == ["file1_visit_01.fits", "file2_visit_01.fits"]
    assert result[4] == ["file3_visit_02.fits",]
