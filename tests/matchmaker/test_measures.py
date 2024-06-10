import pandas as pd
import numpy as np
from astropy.cosmology import FlatLambdaCDM

def test_luminosity_distance():
    # Create sample data
    redshift = 0.5
    flux = 1.0
    alpha = -0.7
    output_units = 'W_Hz'
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Call the luminosity_distance function
    result = luminosity_distance(redshift, flux, alpha, output_units, cosmo)

    # Check if the result is correct
    assert np.isclose(result, 1.0869565217391304e+26)

def test_luminosity_distance_with_series():
    # Create sample data
    redshift = pd.Series([0.5, 0.6, 0.7])
    flux = pd.Series([1.0, 2.0, 3.0])
    alpha = -0.7
    output_units = 'W_Hz'
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Call the luminosity_distance function
    result = luminosity_distance(redshift, flux, alpha, output_units, cosmo)

    # Check if the result is correct
    expected_result = pd.Series([1.0869565217391304e+26, 2.1739130434782607e+26, 3.260869565217391e+26])
    assert result.equals(expected_result)

def test_luminosity_distance_with_array():
    # Create sample data
    redshift = np.array([0.5, 0.6, 0.7])
    flux = np.array([1.0, 2.0, 3.0])
    alpha = -0.7
    output_units = 'W_Hz'
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Call the luminosity_distance function
    result = luminosity_distance(redshift, flux, alpha, output_units, cosmo)

    # Check if the result is correct
    expected_result = np.array([1.0869565217391304e+26, 2.1739130434782607e+26, 3.260869565217391e+26])
    assert np.allclose(result, expected_result)

# Run the tests
test_luminosity_distance()
test_luminosity_distance_with_series()
test_luminosity_distance_with_array()