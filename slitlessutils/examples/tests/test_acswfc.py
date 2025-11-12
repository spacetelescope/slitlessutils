import pytest

from slitlessutils.examples import acswfc


@pytest.mark.skip(reason="Example not working yet")
def test_acswfc():
    acswfc.run_all(plot=False)
