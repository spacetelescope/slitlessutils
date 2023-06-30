import pytest

from slitlessutils.examples import starfield


@pytest.mark.skip(reason="Example not working yet")
def test_starfield():
    starfield.run_all()
