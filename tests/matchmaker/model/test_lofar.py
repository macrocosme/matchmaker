import pytest
from matchmaker.model.lofar import broken_power_law

def test_broken_power_law():
    # Test case 1: psi <= psi_break
    psi = 0.005
    mstar = 1e10
    expected_lum = 22.02 * psi**0.52 * (mstar / 1e10)**0.44
    assert broken_power_law(psi, mstar) == pytest.approx(expected_lum)

    # Test case 2: psi > psi_break
    psi = 0.02
    mstar = 1e10
    expected_lum = 22.02 * psi**1.01 * (mstar / 1e10)**0.44 * 0.01**(0.52 - 1.01)
    assert broken_power_law(psi, mstar) == pytest.approx(expected_lum)

    # Test case 3: psi == psi_break
    psi = 0.01
    mstar = 1e10
    expected_lum = 22.02 * psi**1.01 * (mstar / 1e10)**0.44
    assert broken_power_law(psi, mstar) == pytest.approx(expected_lum)

    # Test case 4: mstar = 0
    psi = 0.01
    mstar = 0
    expected_lum = 0
    assert broken_power_law(psi, mstar) == pytest.approx(expected_lum)

# Run the tests
test_broken_power_law()