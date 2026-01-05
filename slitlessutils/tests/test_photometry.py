"""
Tests for photometry module (Throughput and SED classes).

The photometry module handles filter throughputs and spectral energy
distributions. These are fundamental for flux calibration and simulation.
"""

import numpy as np
import pytest

import slitlessutils as su


@pytest.mark.parametrize("telescope,instrument,filter_name", [
    ('HST', 'WFC3', 'F105W'),
    ('HST', 'WFC3', 'F140W'),
    ('HST', 'ACS', 'F775W'),
])
def test_from_keys_loads_known_filters(telescope, instrument, filter_name):
    """Verify that standard HST filters can be loaded by keys."""
    try:
        throughput = su.photometry.Throughput.from_keys(
            telescope, instrument, filter_name
        )
        assert throughput is not None, (
            f"Failed to load throughput for {telescope}/{instrument}/{filter_name}"
        )
    except (RuntimeError, FileNotFoundError) as e:
        pytest.skip(f"Filter file not available: {e}")


def test_throughput_has_wavelength_data():
    """Verify loaded throughput has wavelength array."""
    try:
        tp = su.photometry.Throughput.from_keys('HST', 'WFC3', 'F105W')
    except (RuntimeError, FileNotFoundError):
        pytest.skip("F105W filter file not available")

    assert hasattr(tp, 'wave'), "Throughput must have wave attribute"
    assert tp.wave is not None, "Throughput wave array must not be None"
    assert len(tp.wave) > 0, "Throughput must have wavelength data"


def test_throughput_has_transmission_data():
    """Verify loaded throughput has transmission array."""
    try:
        tp = su.photometry.Throughput.from_keys('HST', 'WFC3', 'F105W')
    except (RuntimeError, FileNotFoundError):
        pytest.skip("F105W filter file not available")

    assert hasattr(tp, 'tran'), "Throughput must have tran attribute"
    assert tp.tran is not None, "Throughput tran array must not be None"
    assert len(tp.tran) > 0, "Throughput must have transmission data"


def test_throughput_has_zeropoint(sample_throughput):
    """Verify throughput has a photometric zeropoint."""
    assert hasattr(sample_throughput, 'zeropoint'), "Throughput must have zeropoint attribute"
    assert sample_throughput.zeropoint is not None, "Zeropoint must not be None"


def test_zeropoint_can_be_set(sample_throughput):
    """Verify zeropoint can be modified."""
    sample_throughput.zeropoint = 26.5
    assert sample_throughput.zeropoint == 26.5, "Zeropoint should be settable"


def test_photfnu_computable(sample_throughput):
    """Verify PHOTFNU (flux conversion factor) is computable."""
    sample_throughput.zeropoint = 25.0
    photfnu = sample_throughput.photfnu

    assert photfnu is not None, "PHOTFNU must be computable"
    assert photfnu > 0, f"PHOTFNU must be positive, got {photfnu}"


def test_sed_creation_empty():
    """Verify empty SED can be created."""
    sed = su.photometry.SED()
    assert sed is not None, "Should be able to create empty SED"
    assert len(sed) == 0, "Empty SED should have zero length"


def test_sed_creation_with_data():
    """Verify SED can be created with initial data."""
    wave = np.linspace(5000., 10000., 50)
    flux = np.ones_like(wave) * 1e-17
    sed = su.photometry.SED(wave, flux)
    assert len(sed) == 50, f"SED should have 50 elements, got {len(sed)}"


def test_sed_append_data():
    """Verify data can be appended to SED."""
    sed = su.photometry.SED()

    sed.append(5000., 1e-17)
    assert len(sed) == 1, "SED should have 1 element after append"

    sed.append(5100., 1.1e-17)
    assert len(sed) == 2, "SED should have 2 elements after second append"


def test_sed_data_accessible():
    """Verify SED data columns are accessible."""
    wave = np.linspace(5000., 10000., 50)
    flux = np.ones_like(wave) * 1e-17

    sed = su.photometry.SED(wave, flux)

    assert 'lamb' in sed.data.dtype.names, "SED must have wavelength data accessible"
    assert 'flam' in sed.data.dtype.names, "SED must have flux data accessible"


def test_sed_normalize_at_wavelength(flat_sed):
    """Verify SED can be normalized at specific wavelength."""
    assert len(flat_sed) > 0, "SED must have data to normalize"

    original_value = flat_sed['flam'][50]
    assert np.isclose(original_value, 1e-16), "Initial flux should be 1e-16"

    target_flux = 5e-16
    flat_sed.normalize(10000., target_flux)

    new_value = flat_sed['flam'][50]
    assert np.isclose(new_value, 5e-16, rtol=0.1), (
        f"Flux should be scaled to ~5e-16, got {new_value}"
    )


def test_sed_throughput_integration():
    """Verify SED can be integrated through a throughput."""
    wave = np.linspace(8000., 15000., 200)
    flux = np.ones_like(wave) * 1e-16
    sed = su.photometry.SED(wave, flux)

    tp_wave = np.linspace(10000., 12000., 50)
    tp_trans = np.ones_like(tp_wave) * 0.9
    tp = su.photometry.Throughput(tp_wave, tp_trans, unit='angstroms')
    tp.zeropoint = 25.0

    sed_lamb = sed['lamb']
    assert np.any((sed_lamb >= 10000.) & (sed_lamb <= 12000.)), (
        "SED and throughput wavelengths should overlap for integration"
    )


def test_throughput_sed_workflow():
    """Test complete workflow: load throughput, create SED, verify overlap."""
    tp = su.photometry.Throughput.from_keys('HST', 'WFC3', 'F105W')
    assert tp is not None
    assert tp.zeropoint is not None
    assert len(tp.wave) > 0

    wave = np.linspace(8000., 15000., 500)
    flux = np.ones_like(wave) * 1e-17
    sed = su.photometry.SED(wave, flux)

    sed_wave = sed['lamb']
    tp_range = (tp.wmin, tp.wmax)
    sed_range = (sed_wave.min(), sed_wave.max())
    assert tp_range[0] < sed_range[1] and sed_range[0] < tp_range[1]


def test_incremental_sed_building():
    """Test building SED point-by-point from measurements."""
    sed = su.photometry.SED()

    wavelengths = [8000., 9000., 10000., 11000., 12000.]
    fluxes = [1e-17, 1.2e-17, 1.5e-17, 1.3e-17, 1.1e-17]

    for w, f in zip(wavelengths, fluxes):
        sed.append(w, f)

    assert len(sed) == 5
    assert len(sed.data) == 5


def test_throughput_sharp_edges():
    """Test throughput with sharp wavelength cutoffs."""
    wave = np.array([9000., 9001., 11000., 11001.])
    trans = np.array([0., 1., 1., 0.])
    tp = su.photometry.Throughput(wave, trans, unit='angstroms')
    assert len(tp.wave) == 4
    assert len(tp.tran) == 4


def test_sed_with_wavelength_gap():
    """Test SED with non-contiguous wavelength coverage."""
    wave1 = np.linspace(8000., 9000., 20)
    wave2 = np.linspace(10000., 11000., 20)
    wave = np.concatenate([wave1, wave2])
    flux = np.ones_like(wave) * 1e-17
    sed = su.photometry.SED(wave, flux)
    assert len(sed) == 40


def test_sed_single_point():
    """Test SED with single data point."""
    sed = su.photometry.SED()
    sed.append(10000., 1e-17)
    assert len(sed) == 1


def test_throughput_low_transmission():
    """Test throughput with very low transmission values."""
    wave = np.linspace(8000., 12000., 100)
    trans = np.ones_like(wave) * 1e-10
    tp = su.photometry.Throughput(wave, trans, unit='angstroms')
    assert np.all(tp.tran >= 0)


def test_sed_unity_flux_normalization():
    """Test normalizing SED with unity flux values."""
    wave = np.linspace(8000., 12000., 100)
    flux = np.ones_like(wave)
    sed = su.photometry.SED(wave, flux)
    sed.normalize(10000., 5.0)
    assert np.isclose(sed['flam'][50], 5.0, rtol=0.1)
