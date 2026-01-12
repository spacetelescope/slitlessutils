import os

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

from slitlessutils import config

cfg = config.Config()

reffile = cfg.retrieve_reffiles(update=True)


# Copied over from https://github.com/spacetelescope/ci_watson
@pytest.fixture(scope='function')
def _jail(tmp_path):
    """Perform test in a pristine temporary working directory."""
    old_dir = os.getcwd()
    os.chdir(tmp_path)
    try:
        yield str(tmp_path)
    finally:
        os.chdir(old_dir)


@pytest.fixture(scope='session')
def test_data_dir(tmp_path_factory):
    """
    Create a session-scoped temporary directory for test data.

    This is useful for expensive fixtures that can be shared across tests.
    """
    return tmp_path_factory.mktemp("slitlessutils_test_data")


@pytest.fixture
def basic_wcs_header():
    """
    Create a basic WCS header for testing.

    Returns a header with valid WCS for a 100x100 image centered
    at RA=53, Dec=-27 with 0.05 arcsec/pixel scale.
    """
    w = WCS(naxis=2)
    w.wcs.crpix = [50., 50.]
    w.wcs.crval = [53.0, -27.0]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cd = [[-0.05 / 3600., 0.], [0., 0.05 / 3600.]]
    hdr = w.to_header()

    # Standard HST keywords
    hdr['TELESCOP'] = 'HST'
    hdr['INSTRUME'] = 'WFC3'
    hdr['FILTER'] = 'F105W'

    return hdr


@pytest.fixture
def simple_segmentation_and_direct(tmp_path, basic_wcs_header):
    """
    Create simple segmentation and direct images for testing.

    Returns tuple of (segmentation_path, direct_image_path).

    The images contain 2 point-like sources:
    - Source 1: at pixel (30, 30), bright
    - Source 2: at pixel (70, 70), fainter
    """
    npix = 100
    img = np.zeros((npix, npix), dtype=np.float32)
    seg = np.zeros((npix, npix), dtype=np.int32)

    # Create source footprints
    y, x = np.ogrid[:npix, :npix]

    # Source 1
    r1 = np.sqrt((x - 30)**2 + (y - 30)**2)
    img[r1 < 5] = 100.0
    seg[r1 < 3] = 1

    # Source 2
    r2 = np.sqrt((x - 70)**2 + (y - 70)**2)
    img[r2 < 5] = 50.0
    seg[r2 < 3] = 2

    seg_path = str(tmp_path / 'test_seg.fits')
    sci_path = str(tmp_path / 'test_sci.fits')

    fits.writeto(seg_path, seg, header=basic_wcs_header, overwrite=True)
    fits.writeto(sci_path, img, header=basic_wcs_header, overwrite=True)

    return seg_path, sci_path


@pytest.fixture
def sample_wcs_csv_content():
    """
    Return sample CSV content for simulated WFSSCollection tests.

    This CSV defines 3 simulated observations at different orientations.
    """
    return """dataset,ra,dec,orientat,telescope,instrument,disperser,blocking
sim_a,53.162,-27.791,0.0,HST,WFC3IR,G102,
sim_b,53.162,-27.791,60.0,HST,WFC3IR,G102,
sim_c,53.162,-27.791,120.0,HST,WFC3IR,G102,
"""


@pytest.fixture
def sample_wcs_csv(tmp_path, sample_wcs_csv_content):
    """Create a sample WCS CSV file for testing WFSSCollection."""
    csv_file = tmp_path / "test_wcs.csv"
    csv_file.write_text(sample_wcs_csv_content)
    return str(csv_file)


@pytest.fixture
def wfss_detector():
    """Create a simulated WFSS and return the first detector for testing."""
    import slitlessutils as su
    wfss = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test',
        53.0, -27.0, 0.0, 'G102',
        exptime=1000., background=0.5
    )
    for detname, detector in wfss.items():
        return detname, detector


@pytest.fixture
def sample_throughput():
    """Create a simple top-hat throughput for testing."""
    import slitlessutils as su
    wave = np.linspace(10000., 11000., 100)
    trans = np.ones_like(wave)
    trans[:10] = 0.
    trans[-10:] = 0.
    return su.photometry.Throughput(wave, trans, unit='angstroms')


@pytest.fixture
def flat_sed():
    """Create a flat SED for testing with proper flux values."""
    import slitlessutils as su
    wave = np.linspace(8000., 12000., 100)
    flux = np.ones_like(wave) * 1e-16
    return su.photometry.SED(wave, flux)


@pytest.fixture
def slitlessutils_config():
    """Return the Config singleton."""
    return config.Config()
