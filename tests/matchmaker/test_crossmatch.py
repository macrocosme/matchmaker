import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from matchmaker.crossmatch import crossmatch, Ellipse

def test_crossmatch():
    # Create sample catalogs
    obj1 = Catalog(pd.DataFrame({
        'ra': np.random.uniform(1, 5, 5),
        'dec': np.random.uniform(10, 50, 5)
    }), name='obj1')

    obj2 = Catalog(pd.DataFrame({
        'ra': np.random.uniform(2, 6, 5),
        'dec': np.random.uniform(20, 60, 5)
    }), name='obj2')

    # Call the crossmatch function
    crossmatch(obj1, obj2)

    # Check if the matches are correct
    assert obj1.get_matches('obj2') == [0, 1, 2, 3]
    assert obj2.get_matches('obj1') == [1, 2, 3, 4]

def test_crossmatch_with_ellipse_filter():
    # Create sample catalogs
    obj1 = Catalog(pd.DataFrame({
        'ra': np.random.uniform(1, 5, 5),
        'dec': np.random.uniform(10, 50, 5)
    }), name='obj1')

    obj2 = Catalog(pd.DataFrame({
        'ra': np.random.uniform(2, 6, 5),
        'dec': np.random.uniform(20, 60, 5)
    }), name='obj2')

    # Call the crossmatch function with ellipse filter
    crossmatch(obj1, obj2, ellipse_filter=True)

    # Check if the matches are correct
    assert obj1.get_matches('obj2') == [0, 1, 2, 3]
    assert obj2.get_matches('obj1') == [1, 2, 3, 4]

def test_crossmatch_with_mask():
    # Create sample catalogs
    obj1 = Catalog(pd.DataFrame({
        'ra': np.random.uniform(1, 5, 5),
        'dec': np.random.uniform(10, 50, 5)
    }), name='obj1')

    obj2 = Catalog(pd.DataFrame({
        'ra': np.random.uniform(2, 6, 5),
        'dec': np.random.uniform(20, 60, 5)
    }), name='obj2')

    # Create masks
    mask1 = np.array([True, False, True, False, True])
    mask2 = np.array([False, True, False, True, False])

    # Call the crossmatch function with masks
    crossmatch(obj1, obj2, mask1=mask1, mask2=mask2)

    # Check if the matches are correct
    assert obj1.get_matches('obj2') == [0, 2, 4]
    assert obj2.get_matches('obj1') == [1, 3]

# Run the tests
test_crossmatch()
test_crossmatch_with_ellipse_filter()
test_crossmatch_with_mask()