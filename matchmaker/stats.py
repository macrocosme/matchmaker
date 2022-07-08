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
