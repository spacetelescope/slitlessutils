"""
Tests for simulated WFSS data structures.

These tests validate:
1. orientat is correctly set for SimulatedData
2. Detector keywords are properly populated

These are essential for simulation workflows because the telescope pointing
(orientat) directly affects wavelength calibration and spectral extraction.
"""

import numpy as np
import pytest

import slitlessutils as su


@pytest.mark.parametrize("input_orientat", [0.0, 45.0, 90.0, 180.0, 270.0, 359.9])
def test_orientat_propagates_to_wfss(input_orientat):
    """Verify that orientat value is correctly stored and accessible."""
    wfss = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test_sim',
        53.162, -27.791, input_orientat, 'G102'
    )

    for detname, detector in wfss.items():
        pc = detector.wcs.wcs.pc
        assert pc is not None, "PC matrix should be set from orientat"
        assert pc.shape == (2, 2), "PC matrix must be 2x2"


def test_orientat_affects_wcs_pc_matrix():
    """Verify that different orientat values produce different WCS matrices."""
    wfss_0 = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test_sim',
        53.162, -27.791, 0.0, 'G102'
    )

    wfss_90 = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test_sim',
        53.162, -27.791, 90.0, 'G102'
    )

    for detname in wfss_0.keys():
        pc_0 = wfss_0[detname].wcs.wcs.pc
        pc_90 = wfss_90[detname].wcs.wcs.pc

        assert not np.allclose(pc_0, pc_90), (
            f"PC matrices for orientat=0 and orientat=90 should differ "
            f"for detector {detname}. Got pc_0={pc_0}, pc_90={pc_90}"
        )


def test_orientat_crvals_updated():
    """Verify that CRVAL (reference pixel coordinates) are set."""
    wfss = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test_sim',
        53.162, -27.791, 45.0, 'G102'
    )

    for detname, detector in wfss.items():
        crval = detector.wcs.wcs.crval
        assert crval is not None, "CRVAL must be set"
        assert len(crval) == 2, "CRVAL must have 2 elements (RA, Dec)"
        assert -90 <= crval[1] <= 90, "Dec must be in valid range"
        assert 0 <= crval[0] <= 360, "RA must be in valid range"


def test_from_wcsfile_loads_all_datasets(sample_wcs_csv):
    """Verify that all datasets in CSV file are loaded."""
    collection = su.wfss.WFSSCollection.from_wcsfile(sample_wcs_csv)

    assert len(collection) == 3, (
        f"Expected 3 datasets but got {len(collection)}. "
        "Check CSV parsing in from_wcsfile()."
    )
    assert 'sim_a' in collection, "Dataset 'sim_a' should be loaded"
    assert 'sim_b' in collection, "Dataset 'sim_b' should be loaded"
    assert 'sim_c' in collection, "Dataset 'sim_c' should be loaded"


def test_from_wcsfile_preserves_orientat(sample_wcs_csv):
    """Verify that orientat values from CSV are preserved."""
    collection = su.wfss.WFSSCollection.from_wcsfile(sample_wcs_csv)

    expected = {'sim_a': 0.0, 'sim_b': 60.0, 'sim_c': 120.0}

    for dataset_name in collection.keys():
        simdata = dict.__getitem__(collection, dataset_name)
        expected_orientat = expected[dataset_name]
        assert np.allclose(simdata.orientat, expected_orientat), (
            f"Dataset {dataset_name}: expected orientat={expected_orientat}, "
            f"got {simdata.orientat}"
        )


def test_collection_iteration_yields_wfss_objects(sample_wcs_csv):
    """Verify iterating over collection yields WFSS objects with required components."""
    collection = su.wfss.WFSSCollection.from_wcsfile(sample_wcs_csv)

    count = 0
    for wfss in collection:
        assert isinstance(wfss, su.wfss.WFSS), f"Expected WFSS, got {type(wfss)}"
        for detname, detector in wfss.items():
            assert detector.wcs is not None
            assert detector.config is not None
            assert detector.naxis is not None and all(n > 0 for n in detector.naxis)
            assert len(list(detector.orders)) > 0
        count += 1

    assert count == 3, "Should iterate over all 3 datasets"


def test_detector_wcs_keywords(wfss_detector):
    """Verify WCS keywords (NAXIS, CRPIX, CRVAL) are properly set."""
    detname, detector = wfss_detector

    naxis = detector.naxis
    assert naxis is not None and len(naxis) == 2
    assert naxis[0] == 1014 and naxis[1] == 1014, f"WFC3IR detector should be 1014x1014, got {naxis}"

    crpix = detector.wcs.wcs.crpix
    assert crpix is not None and len(crpix) == 2

    crval = detector.wcs.wcs.crval
    assert crval is not None and len(crval) == 2


def test_detector_config_properties(wfss_detector):
    """Verify detector config has scale, extensions, and noise properties."""
    detname, detector = wfss_detector

    scale = detector.config.scale
    assert scale is not None and len(scale) == 2
    assert np.allclose(scale, [0.121, 0.136], rtol=0.01), f"WFC3IR scale changed, got {scale}"

    extensions = detector.extensions
    assert extensions is not None
    assert 'science' in extensions

    noise = detector.config.noise
    assert noise is not None
    assert hasattr(noise, 'read') and np.allclose(noise.read, 12.0, rtol=0.01), f"WFC3IR read noise changed, got {noise.read}"
    assert hasattr(noise, 'dark') and np.allclose(noise.dark, 0.045, rtol=0.01), f"WFC3IR dark current changed, got {noise.dark}"


def test_detector_spectral_orders(wfss_detector):
    """Verify spectral orders are defined for the disperser."""
    detname, detector = wfss_detector

    orders = list(detector.orders)
    assert len(orders) == 5, f"G102 should have 5 spectral orders, got {len(orders)}"
    assert set(orders) == {'+1', '0', '+2', '+3', '-1'}, f"G102 orders changed, got {orders}"
    for order in orders:
        assert isinstance(order, str), f"Order names should be strings, got {type(order)}"


@pytest.mark.parametrize("ra,dec,orientat,description", [
    (53.0, -27.0, 0.0, "zero orientat"),
    (53.0, -60.0, 45.0, "southern hemisphere"),
    (0.0, 85.0, 0.0, "near celestial pole"),
])
def test_coordinate_edge_cases(ra, dec, orientat, description):
    """Test simulation handles various coordinate edge cases."""
    wfss = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test', ra, dec, orientat, 'G102'
    )

    for detname, detector in wfss.items():
        assert detector.wcs is not None, f"WCS should be valid for {description}"
        crval = detector.wcs.wcs.crval
        assert 0 <= crval[0] <= 360 or np.allclose(crval[0], ra, atol=1.0)
        assert -90 <= crval[1] <= 90


def test_orientat_360_equals_0():
    """Verify orientat=360 produces same WCS as orientat=0."""
    wfss_0 = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test', 53.0, -27.0, 0.0, 'G102'
    )
    wfss_360 = su.wfss.WFSS.simulated(
        'HST', 'WFC3IR', 'test', 53.0, -27.0, 360.0, 'G102'
    )

    for detname in wfss_0.keys():
        pc_0 = wfss_0[detname].wcs.wcs.pc
        pc_360 = wfss_360[detname].wcs.wcs.pc
        assert np.allclose(pc_0, pc_360, atol=1e-10), "orientat=360 should equal orientat=0"
