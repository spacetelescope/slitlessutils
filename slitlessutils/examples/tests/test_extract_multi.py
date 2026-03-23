import pytest

from slitlessutils.examples import simulate, extract_multi

def test_extract_multi():
    # Need the outputs of simulate to run extract_multi
    simulate.run_all()
    extract_multi.run_all()
