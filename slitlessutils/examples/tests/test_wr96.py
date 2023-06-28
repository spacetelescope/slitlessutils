import pytest

from slitlessutils.examples import wr96


@pytest.mark.skip(reason="Example not working yet")
def test_wr96():
    wr96.run_all(plot=False)
