import math

import astropy.units
import pandas as pd
import numpy as np
import astropy.units as u
from .dataset import Catalog

chance_association_probability = lambda n, theta, area: n * (math.pi * theta**2) / area
chance_association_probability_ab = lambda n, a, b, area: n * (math.pi * a * b) / area
expected_background_source = lambda chance_association, n: n * chance_association

# TODO: chance association when not using 'a' in clu

def chance_association(obj1:Catalog,
                       obj1_ids:np.ndarray,
                       theta:astropy.units.Quantity=6 * u.arcmin,
                       n_in_area=700,
                       area=1 * u.deg**2,
                       catalogue_name='lotss'):
    """
    Calculate the chance association for a given set of objects.

    Parameters:
    - obj1 (Catalog): The catalog object containing the objects.
    - obj1_ids (np.ndarray): The array of object IDs.
    - theta (astropy.units.Quantity): The angular separation threshold.
    - n_in_area (int): The number of objects in the search area.
    - area (astropy.units.Quantity): The search area.
    - catalogue_name (str): The name of the catalogue.

    Returns:
    None
    """
    total = obj1_ids.size

    expected = chance_association_probability(
        n=n_in_area,
        theta=theta.to(u.deg),
        area=area
    ).value

    obj1.matches[catalogue_name].set_chance(total, total * expected)

def chance_association_ab(obj1:Catalog,
                          obj1_ids:np.ndarray,
                          n_in_area=700,
                          area=1 * u.deg**2,
                          catalogue_name='lotss',
                          verbose=False):
    """
    Calculate the chance association between objects in obj1 and a given catalog.

    Parameters:
    - obj1 (Catalog): The catalog containing the objects.
    - obj1_ids (np.ndarray): The array of object IDs in obj1.
    - n_in_area (int, optional): The number of objects in the area. Default is 700.
    - area (Quantity, optional): The area of the region. Default is 1 square degree.
    - catalogue_name (str, optional): The name of the catalog. Default is 'lotss'.
    - verbose (bool, optional): Whether to print verbose output. Default is False.

    Returns:
    None
    """
    expected = 0
    total = 0

    for obj1_id in obj1_ids:
        partial = chance_association_probability_ab(
            n=n_in_area,
            a=obj1.semi_major(with_unit=True, to_unit=u.deg)[obj1_id],
            b=obj1.semi_minor(with_unit=True, to_unit=u.deg)[obj1_id],
            area=area
        ).value

        if not pd.isna(partial):
            total += 1
            expected += partial

        if verbose:
            print ("""Total expected false positives from obj1_ids so far: {} 
Expected false association(s) for region: {:0.2f} 
Region:
    semi-major: {:0.2f}
    semi-minor: {:0.2f}
""".format(np.round(expected).astype(int),
           partial,
           obj1.semi_major(with_unit=True)[obj1_id],
           obj1.semi_minor(with_unit=True)[obj1_id])
            )

    obj1.matches[catalogue_name].set_chance(total, expected)
