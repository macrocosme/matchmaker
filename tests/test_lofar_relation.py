"""The radio luminosity - star formation relation: the broken power law used to be unusable."""
import numpy as np

from matchmaker.model.lofar import broken_power_law, measured_power_law, power_law


def test_power_law_matches_gurkan_normalisation():
    assert np.isclose(np.log10(power_law(1.0)), 22.06)


def test_power_law_accepts_arrays():
    out = power_law(np.array([0.1, 1.0, 10.0]))
    assert out.shape == (3,)
    assert np.all(np.diff(out) > 0)


def test_broken_power_law_is_a_luminosity_not_a_factor_of_22():
    """It used to multiply by lc = 22.02 instead of 10**22.02, so it returned ~1e-2 W/Hz."""
    value = broken_power_law(1.0, 1e10)
    assert 1e21 < value < 1e23


def test_broken_power_law_is_continuous_at_the_break():
    below = broken_power_law(0.01 * (1 - 1e-9), 1e9)
    above = broken_power_law(0.01 * (1 + 1e-9), 1e9)
    assert np.isclose(below, above, rtol=1e-6)


def test_broken_power_law_mass_term_lowers_dwarfs():
    """A 1e9 Msun galaxy gets (0.1)^0.44 of a 1e10 one, a factor 2.8 down."""
    ratio = broken_power_law(1.0, 1e9) / broken_power_law(1.0, 1e10)
    assert np.isclose(ratio, 0.1 ** 0.44)


def test_broken_power_law_accepts_arrays():
    out = broken_power_law(np.array([0.001, 0.1, 1.0]), np.array([1e8, 1e9, 1e10]))
    assert out.shape == (3,)
    assert np.all(np.isfinite(out))


def test_measured_power_law_reduces_to_the_assumed_one():
    assert np.isclose(measured_power_law(0.3, 22.06, 1.07), power_law(0.3))
