"""
Tests for the slitlessutils configuration module.

The Config class is a singleton that manages global settings and reference
files. These tests verify that configuration is correctly loaded and
accessible, which is essential for all downstream operations.
"""

import os

from slitlessutils import config


def test_config_singleton_returns_same_instance():
    """Verify that Config() always returns the same instance."""
    cfg1 = config.Config()
    cfg2 = config.Config()

    assert cfg1 is cfg2, (
        "Config() should return same instance (singleton pattern). "
        "Got different instances."
    )


def test_config_is_dict_like():
    """Verify Config behaves like a dictionary."""
    cfg = config.Config()

    assert hasattr(cfg, '__getitem__'), "Config must support item access"
    assert hasattr(cfg, '__setitem__'), "Config must support item assignment"
    assert hasattr(cfg, 'keys'), "Config must support keys()"


def test_fluxscale_exists_and_positive(slitlessutils_config):
    """Verify fluxscale parameter exists and is positive."""
    assert 'fluxscale' in slitlessutils_config or hasattr(slitlessutils_config, 'fluxscale'), (
        "fluxscale parameter must exist in Config"
    )

    fluxscale = slitlessutils_config.fluxscale
    assert fluxscale is not None, "fluxscale must not be None"
    assert fluxscale > 0, f"fluxscale must be positive, got {fluxscale}"


def test_fluxunits_exists(slitlessutils_config):
    """Verify fluxunits parameter exists."""
    assert 'fluxunits' in slitlessutils_config or hasattr(slitlessutils_config, 'fluxunits'), (
        "fluxunits parameter must exist in Config"
    )

    fluxunits = slitlessutils_config.fluxunits
    assert fluxunits is not None, "fluxunits must not be None"
    assert isinstance(fluxunits, str), f"fluxunits must be a string, got {type(fluxunits)}"


def test_compression_settings_exist(slitlessutils_config):
    """Verify HDF5 compression settings exist."""
    assert hasattr(slitlessutils_config, 'compression'), "compression setting must exist"
    assert hasattr(slitlessutils_config, 'compression_opts'), "compression_opts setting must exist"


def test_refpath_exists(slitlessutils_config):
    """Verify reference file path is set."""
    assert hasattr(slitlessutils_config, 'refpath'), "refpath attribute must exist"
    refpath = slitlessutils_config.refpath
    assert refpath is not None, "refpath must not be None"
    assert isinstance(refpath, str), f"refpath must be string, got {type(refpath)}"


def test_refpath_directory_exists(slitlessutils_config):
    """Verify reference file directory exists on disk."""
    refpath = slitlessutils_config.refpath
    assert os.path.isdir(refpath), f"Reference path directory must exist: {refpath}"


def test_instruments_directory_exists(slitlessutils_config):
    """Verify instruments subdirectory exists."""
    instruments_path = os.path.join(slitlessutils_config.refpath, 'instruments')
    assert os.path.isdir(instruments_path), f"Instruments directory must exist: {instruments_path}"


def test_suffixes_defined():
    """Verify standard file suffixes are defined."""
    from slitlessutils.config import SUFFIXES

    expected_keys = ['1d spectra', '2d spectra', 'L-curve', 'group']

    for key in expected_keys:
        assert key in SUFFIXES, f"SUFFIXES must include '{key}'"


def test_suffixes_are_strings():
    """Verify all suffixes are non-empty strings."""
    from slitlessutils.config import SUFFIXES

    for key, suffix in SUFFIXES.items():
        assert isinstance(suffix, str), f"Suffix for '{key}' must be string, got {type(suffix)}"
        assert len(suffix) > 0, f"Suffix for '{key}' must not be empty"
