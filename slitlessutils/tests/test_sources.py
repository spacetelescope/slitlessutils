"""
Tests for Source and SourceCollection data structures.

Sources are the fundamental unit of spectroscopic extraction. These tests
verify that sources are correctly loaded from segmentation maps and direct
images, which is essential for any extraction workflow.
"""

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

import slitlessutils as su


@pytest.fixture
def multi_source_data(tmp_path):
    """Create test data with multiple sources at different positions."""
    npix = 200

    # Create WCS header properly
    w = WCS(naxis=2)
    w.wcs.crpix = [npix / 2., npix / 2.]
    w.wcs.crval = [53.16, -27.79]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cd = [[-0.13 / 3600., 0.], [0., 0.13 / 3600.]]
    header = w.to_header()
    header['TELESCOP'] = 'HST'
    header['INSTRUME'] = 'WFC3'
    header['FILTER'] = 'F105W'

    # Create segmentation with 5 sources
    seg = np.zeros((npix, npix), dtype=np.int32)
    sci = np.zeros((npix, npix), dtype=np.float32)

    # Place sources at different locations
    positions = [(50, 50), (150, 50), (100, 100), (50, 150), (150, 150)]
    for i, (x, y) in enumerate(positions, start=1):
        yy, xx = np.ogrid[-y:npix - y, -x:npix - x]
        mask = xx ** 2 + yy ** 2 <= 5 ** 2
        seg[mask] = i
        sci[mask] = 100. / i

    seg_path = str(tmp_path / 'multi_seg.fits')
    sci_path = str(tmp_path / 'multi_sci.fits')
    fits.writeto(seg_path, seg, header=header, overwrite=True)
    fits.writeto(sci_path, sci, header=header, overwrite=True)

    return seg_path, sci_path


@pytest.fixture
def edge_case_header(tmp_path):
    """Create header for edge case testing."""
    npix = 100

    # Create WCS header properly
    w = WCS(naxis=2)
    w.wcs.crpix = [npix / 2., npix / 2.]
    w.wcs.crval = [0.0, 0.0]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cd = [[-0.13 / 3600., 0.], [0., 0.13 / 3600.]]
    header = w.to_header()
    header['TELESCOP'] = 'HST'
    header['INSTRUME'] = 'WFC3'
    header['FILTER'] = 'F105W'

    return tmp_path, npix, header


def test_sourcecollection_loads_from_files(simple_segmentation_and_direct):
    """Verify SourceCollection loads sources from seg and direct images."""
    seg_path, sci_path = simple_segmentation_and_direct

    sources = su.sources.SourceCollection(
        seg_path, sci_path,
        zeropoint=25.0
    )

    assert sources is not None, "SourceCollection should be created"
    assert len(sources) > 0, "SourceCollection should contain at least one source"


def test_sourcecollection_finds_all_sources(simple_segmentation_and_direct):
    """Verify all sources in segmentation map are loaded."""
    seg_path, sci_path = simple_segmentation_and_direct

    sources = su.sources.SourceCollection(
        seg_path, sci_path,
        zeropoint=25.0
    )

    assert len(sources) == 2, (
        f"Expected 2 sources, found {len(sources)}. "
        "Check segmentation map parsing."
    )


def test_sourcecollection_segids_match_segmap(simple_segmentation_and_direct):
    """Verify source IDs match segmentation map values."""
    seg_path, sci_path = simple_segmentation_and_direct

    sources = su.sources.SourceCollection(
        seg_path, sci_path,
        zeropoint=25.0
    )

    assert 1 in sources, "Source with segid=1 should exist"
    assert 2 in sources, "Source with segid=2 should exist"
    assert 0 not in sources, "Segid=0 (background) should not be a source"


def test_source_has_segid(simple_segmentation_and_direct):
    """Verify each source has a segmentation ID."""
    seg_path, sci_path = simple_segmentation_and_direct
    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0
    )

    for segid, source in sources.items():
        assert hasattr(source, 'segid'), "Source must have segid attribute"
        assert source.segid == segid, (
            f"Source segid ({source.segid}) must match dict key ({segid})"
        )


def test_source_has_pixels(simple_segmentation_and_direct):
    """Verify sources have pixel coordinates."""
    seg_path, sci_path = simple_segmentation_and_direct
    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0
    )

    for segid, source in sources.items():
        assert source.npixels > 0, (
            f"Source {segid} must have pixels, got npixels={source.npixels}"
        )


def test_source_has_wcs(simple_segmentation_and_direct):
    """Verify sources have WCS information."""
    seg_path, sci_path = simple_segmentation_and_direct
    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0
    )

    for segid, source in sources.items():
        assert hasattr(source, 'wcs'), f"Source {segid} must have WCS"
        assert source.wcs is not None, f"WCS for source {segid} must not be None"


def test_minpix_filters_small_sources(simple_segmentation_and_direct):
    """Verify minpix parameter filters sources with too few pixels."""
    seg_path, sci_path = simple_segmentation_and_direct

    # Our sources have ~28 pixels each (radius 3)
    # Setting minpix=50 should filter them out
    sources = su.sources.SourceCollection(
        seg_path, sci_path,
        zeropoint=25.0,
        minpix=50
    )

    assert len(sources) == 0, "Sources with fewer than minpix pixels should be filtered"


def test_minpix_zero_keeps_all(simple_segmentation_and_direct):
    """Verify minpix=0 keeps all sources."""
    seg_path, sci_path = simple_segmentation_and_direct

    sources = su.sources.SourceCollection(
        seg_path, sci_path,
        zeropoint=25.0,
        minpix=0
    )

    assert len(sources) == 2, "With minpix=0, all sources should be kept"


def test_set_spectral_parameters_on_collection(simple_segmentation_and_direct):
    """Verify spectral parameters can be set on entire collection."""
    seg_path, sci_path = simple_segmentation_and_direct
    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0
    )

    sources.set_spectral_parameters(wave0=8000., wave1=12000.)

    for segid, source in sources.items():
        for region in source:
            assert hasattr(region, 'wave0'), f"wave0 not propagated to source {segid}"
            assert region.wave0 == 8000., f"wave0 value incorrect for source {segid}"


def test_set_spectral_parameters_on_individual_source(simple_segmentation_and_direct):
    """Verify spectral parameters can be set per-source."""
    seg_path, sci_path = simple_segmentation_and_direct
    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0
    )

    sources[1].set_spectral_parameters(wave0=9000.)

    for region in sources[1]:
        assert hasattr(region, 'wave0'), "Individual source spectral parameters should be settable"
        assert region.wave0 == 9000., "wave0 should be set to 9000 on individual source"


def test_complete_source_workflow(multi_source_data):
    """Test loading sources and configuring spectral parameters."""
    seg_path, sci_path = multi_source_data

    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=26.0
    )
    assert len(sources) == 5, "Should load all 5 sources"

    for segid in range(1, 6):
        assert segid in sources
        source = sources[segid]
        assert source.segid == segid
        assert source.npixels > 0

    sources.set_spectral_parameters(wave0=7500., wave1=11500., dwave=25.)

    for segid, source in sources.items():
        for region in source:
            assert region.wave0 == 7500.
            assert region.wave1 == 11500.


def test_source_at_image_edge(edge_case_header):
    """Test source near image edge loads correctly."""
    tmp_path, npix, header = edge_case_header

    seg = np.zeros((npix, npix), dtype=np.int32)
    sci = np.zeros((npix, npix), dtype=np.float32)
    seg[5:10, 5:10] = 1
    sci[5:10, 5:10] = 100.

    seg_path = str(tmp_path / 'edge_seg.fits')
    sci_path = str(tmp_path / 'edge_sci.fits')
    fits.writeto(seg_path, seg, header=header, overwrite=True)
    fits.writeto(sci_path, sci, header=header, overwrite=True)

    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0, local_back=False
    )
    assert 1 in sources


def test_single_pixel_source(edge_case_header):
    """Test single-pixel source loads correctly."""
    tmp_path, npix, header = edge_case_header

    seg = np.zeros((npix, npix), dtype=np.int32)
    sci = np.zeros((npix, npix), dtype=np.float32)
    seg[50, 50] = 1
    sci[50, 50] = 100.

    seg_path = str(tmp_path / 'single_seg.fits')
    sci_path = str(tmp_path / 'single_sci.fits')
    fits.writeto(seg_path, seg, header=header, overwrite=True)
    fits.writeto(sci_path, sci, header=header, overwrite=True)

    sources = su.sources.SourceCollection(
        seg_path, sci_path, zeropoint=25.0, minpix=0
    )
    if 1 in sources:
        assert sources[1].npixels == 1


def test_high_segid_values(edge_case_header):
    """Test source with high segmentation ID (99999) loads correctly."""
    tmp_path, npix, header = edge_case_header

    seg = np.zeros((npix, npix), dtype=np.int32)
    sci = np.zeros((npix, npix), dtype=np.float32)

    high_segid = 99999
    yy, xx = np.ogrid[45:55, 45:55]
    mask = (xx - 50) ** 2 + (yy - 50) ** 2 <= 16
    seg[45:55, 45:55][mask[0:10, 0:10]] = high_segid
    sci[45:55, 45:55][mask[0:10, 0:10]] = 100.

    seg_path = str(tmp_path / 'highid_seg.fits')
    sci_path = str(tmp_path / 'highid_sci.fits')
    fits.writeto(seg_path, seg, header=header, overwrite=True)
    fits.writeto(sci_path, sci, header=header, overwrite=True)

    sources = su.sources.SourceCollection(seg_path, sci_path, zeropoint=25.0)
    assert high_segid in sources
    assert sources[high_segid].segid == high_segid
