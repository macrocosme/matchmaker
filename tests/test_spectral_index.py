"""Radio spectral indices."""
import numpy as np

from matchmaker.spectral_index import fit_spectral_index, spectral_indices


def test_spectral_index_recovers_power_law():
    f = np.array([144, 887.5, 1400, 3000.0])
    s = 10 * (f / 144) ** -0.7
    si = spectral_indices(f, s, 0.05 * s)
    assert np.isclose(si['alpha'], -0.7) and np.isclose(si['alpha_high'], -0.7)


def test_spectral_index_turnover_shows_as_alpha_above_alpha_high():
    f = np.array([144, 887.5, 1400, 3000.0])
    s = 10 * (f / 144) ** -0.7
    s[0] *= 0.1                                     # absorbed at 144 MHz
    si = spectral_indices(f, s, 0.05 * s)
    assert si['alpha'] > si['alpha_high']


def test_spectral_index_undefined_is_nan_not_zero():
    a, e = fit_spectral_index([144], [5.0], [0.5])
    assert np.isnan(a) and np.isnan(e)


# --- reference PRS (model/prs.py) ---------------------------------------------------------
