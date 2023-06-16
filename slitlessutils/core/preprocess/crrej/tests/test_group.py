from slitlessutils.core.preprocess.crrej.group import group_by_visit
from slitlessutils.core.wfss import WFSSCollection
import pytest
import numpy as np


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
