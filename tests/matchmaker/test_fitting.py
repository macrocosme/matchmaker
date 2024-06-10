import numpy as np
import pytest

def test_fit_power_law():
    # Test case 1: Simple power law fit
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([1, 4, 9, 16, 25])
    a, k, d_a, d_k = fit_power_law(x, y)
    assert np.isclose(a, 1.0, rtol=1e-3)
    assert np.isclose(k, 2.0, rtol=1e-3)
    assert np.isclose(d_a, 0.0, rtol=1e-3)
    assert np.isclose(d_k, 0.0, rtol=1e-3)

    # Test case 2: Power law fit with uncertainties
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([1, 4, 9, 16, 25])
    y_err = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    a, k, d_a, d_k = fit_power_law(x, y, y_err=y_err)
    assert np.isclose(a, 1.0, rtol=1e-3)
    assert np.isclose(k, 2.0, rtol=1e-3)
    assert np.isclose(d_a, 0.0, rtol=1e-3)
    assert np.isclose(d_k, 0.0, rtol=1e-3)

    # Test case 3: Power law fit with relative uncertainties
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([1, 4, 9, 16, 25])
    y_err = np.array([0.1, 0.2, 0.3, 0.4, 0.5]) / y
    a, k, d_a, d_k = fit_power_law(x, y, y_err=y_err, absolute_sigma=False)
    assert np.isclose(a, 1.0, rtol=1e-3)
    assert np.isclose(k, 2.0, rtol=1e-3)
    assert np.isclose(d_a, 0.0, rtol=1e-3)
    assert np.isclose(d_k, 0.0, rtol=1e-3)

    # Test case 4: Power law fit with non-zero uncertainties
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([1, 4, 9, 16, 25])
    y_err = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    a, k, d_a, d_k = fit_power_law(x, y, y_err=y_err)
    assert np.isclose(a, 1.0, rtol=1e-3)
    assert np.isclose(k, 2.0, rtol=1e-3)
    assert np.isclose(d_a, 0.0, rtol=1e-3)
    assert np.isclose(d_k, 0.0, rtol=1e-3)

    # Test case 5: Power law fit with negative values
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([-1, -4, -9, -16, -25])
    a, k, d_a, d_k = fit_power_law(x, y)
    assert np.isclose(a, -1.0, rtol=1e-3)
    assert np.isclose(k, 2.0, rtol=1e-3)
    assert np.isclose(d_a, 0.0, rtol=1e-3)
    assert np.isclose(d_k, 0.0, rtol=1e-3)

    # Test case 6: Power law fit with zero values
    x = np.array([1, 2, 3, 4, 5])
    y = np.array([0, 0, 0, 0, 0])
    a, k, d_a, d_k = fit_power_law(x, y)
    assert np.isclose(a, 0.0, rtol=1e-3)
    assert np.isclose(k, 0.0, rtol=1e-3)
    assert np.isclose(d_a, 0.0, rtol=1e-3)
    assert np.isclose(d_k, 0.0, rtol=1e-3)

# Run the tests
test_fit_power_law()