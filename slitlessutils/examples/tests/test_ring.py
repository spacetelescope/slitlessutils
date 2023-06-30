import pytest

from slitlessutils.examples import ring


@pytest.mark.skip(reason="Example not working yet")
def test_ring():
    ring.run_all()
