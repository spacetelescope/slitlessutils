import numpy as np
from astropy.io import fits
from numpy.testing import assert_allclose

from slitlessutils.core.preprocess.astrometry.utils import (get_cd, get_crval,
                                                            set_cd, set_crval)


def test_cd_utils():
    fits_header = fits.Header()
    fits_header["CD1_1A"] = 1.
    fits_header["CD1_2A"] = 2.
    fits_header["CD2_1A"] = 3.
    fits_header["CD2_2A"] = 4.

    assert_allclose(get_cd(fits_header, "A"), np.array([[1., 2.], [3., 4.]]))

    new_vals = np.array([[5., 6.], [7., 8.]])
    set_cd(fits_header, new_vals, "A")

    assert fits_header["CD1_1A"] == 5.0
    assert fits_header["CD1_2A"] == 6.0
    assert fits_header["CD2_1A"] == 7.0
    assert fits_header["CD2_2A"] == 8.0


def test_crval_utils():
    fits_header = fits.Header()
    fits_header["CRVAL1A"] = 1.
    fits_header["CRVAL2A"] = 2.

    assert_allclose(get_crval(fits_header, "A"), np.array([1., 2.]))

    new_vals = np.array([3., 4.])
    set_crval(fits_header, new_vals, "A")

    assert fits_header["CRVAL1A"] == 3.0
    assert fits_header["CRVAL2A"] == 4.0
