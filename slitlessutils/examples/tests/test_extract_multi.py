import pytest
from astropy.io import fits
from astropy.table import Table
from numpy.testing import assert_almost_equal

from slitlessutils.examples import extract_multi, simulate


@pytest.mark.filterwarnings("ignore:did not parse as fits unit")
def test_extract_multi():
    # Need the outputs of simulate to run extract_multi
    simulate.run_all()
    extract_multi.run_all()

    x1d = fits.open('pointsources_x1d.fits')
    assert len(x1d) == 6
    assert_almost_equal(x1d[1].data['FLUX'][100], 2.2441835)

    # Make sure this reads into a table. Note that this raises a UnitsWarning
    # due to the constant factor in the unit.
    t = Table(x1d[1])
    assert str(t['FLUX'].unit) == '1e-17 erg/(s*cm**2*AA)'
    assert t['WAVELENGTH'][100] == 10000
    assert_almost_equal(t['FLUX'][110], 2.1282458)
