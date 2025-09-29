import numpy as np
import pytest

from astropy.io import fits
from astropy.utils.data import download_file

from slitlessutils.core.utilities.embedding import embedsub_full_detector


@pytest.mark.remote_data
def test_subarray_embedding(tmp_path):
    ir_file = 'iebo7aqsq_flt.fits'
    ir_box = 'https://stsci.box.com/shared/static/8mcltriruun5u6wvoelpxuv9i8am4gtm.fits'
    uvis_file = 'iels01aaq_flt.fits'
    uvis_box = 'https://stsci.box.com/shared/static/abny5qtrkvayng10rjwd84jocwtkiggd.fits'

    uvis_file = download_file(uvis_box)
    uvis_file
    ir_file = download_file(ir_box)

    # Manually provide x and y size
    embedded_uvis = embedsub_full_detector(uvis_file, y_size=2051, x_size=4096)
    uvis = fits.open(embedded_uvis)
    assert np.mean(uvis[1].data) > 0
    assert uvis[1].data.shape == (2051, 4096)

    # Specify instrument instead
    embedded_uvis = embedsub_full_detector(uvis_file, instrument='UVIS')
    uvis = fits.open(embedded_uvis)
    assert np.mean(uvis[1].data) > 0
    assert uvis[1].data.shape == (2051, 4096)

    embedded_ir = embedsub_full_detector(ir_file, y_size=1014, x_size=1014)
    ir = fits.open(embedded_ir)
    assert np.mean(ir[1].data) > 0
    assert ir[1].data.shape == (1014, 1014)

    with pytest.raises(ValueError, match='Instrument cannot be set if y_size and x_size are set'):
        bad_embed = embedsub_full_detector(ir_file, 'IR', y_size=1014, x_size=1014)  # noqa

    with pytest.raises(ValueError, match='One of instrument or x_size'):
        bad_embed = embedsub_full_detector(ir_file)  # noqa
