import numpy as np
import pandas as pd

from matchmaker.measures import (
    erg_s_hz_to_w_hz,
    luminosity_distance,
    powerlaw_scale,
    w_hz_to_erg_s_hz,
)


def test_powerlaw_scale():
    assert powerlaw_scale(1.0, 2.0, 10.0, 0.0) == 10.0
    assert powerlaw_scale(1.0, 2.0, 10.0, -1.0) == 5.0


def test_unit_conversion_roundtrip():
    # 1 W/Hz = 1e7 erg/s/Hz
    assert np.isclose(w_hz_to_erg_s_hz(1.0), 1e7)
    assert np.isclose(erg_s_hz_to_w_hz(w_hz_to_erg_s_hz(3.5)), 3.5)


def test_luminosity_distance_scalar():
    lum = luminosity_distance(0.01, 1.0)
    assert np.isfinite(lum)
    assert lum > 0


def test_luminosity_distance_ndarray():
    # Regression: type check used np.array (a function) instead of np.ndarray,
    # so ndarray input silently returned None
    lum = luminosity_distance(np.array([0.01, 0.02]), np.array([1.0, 2.0]))
    assert lum is not None
    assert len(lum) == 2
    assert np.all(np.isfinite(lum))


def test_luminosity_distance_series():
    lum = luminosity_distance(pd.Series([0.01, 0.02]), pd.Series([1.0, 2.0]))
    assert lum is not None
    assert len(lum) == 2


def test_luminosity_distance_increases_with_redshift():
    low = luminosity_distance(0.01, 1.0)
    high = luminosity_distance(0.05, 1.0)
    assert high > low
