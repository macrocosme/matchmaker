import numpy as np
import astropy.units as u
from .dataset import Catalog

def test_chance_association_ab():
    # Create sample data
    obj1 = Catalog()
    obj1_ids = np.array([0, 1, 2])
    n_in_area = 700
    area = 1 * u.deg**2
    catalogue_name = 'lotss'
    verbose = False

    # Call the chance_association_ab function
    chance_association_ab(obj1, obj1_ids, n_in_area, area, catalogue_name, verbose)

    # Check if the result is correct
    expected_total = 3
    expected_expected = 0.0
    assert obj1.matches[catalogue_name].total == expected_total
    assert obj1.matches[catalogue_name].expected == expected_expected

# Run the test
test_chance_association_ab()