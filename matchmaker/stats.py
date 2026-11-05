import math

import astropy.units
import astropy.units as u
import numpy as np

from .dataset import Catalog

chance_association_probability = lambda n, theta, area: n * (math.pi * theta**2) / area
chance_association_probability_ab = lambda n, a, b, area: n * (math.pi * a * b) / area
expected_background_source = lambda chance_association, n: n * chance_association

# TODO: chance association when not using 'a' in clu

def _n_searched(obj1, n_searched):
    if n_searched is not None:
        return n_searched
    if obj1.df is None:
        raise ValueError('n_searched is required when the source catalogue has no DataFrame loaded')
    return len(obj1.df)


def chance_association(obj1:Catalog,
                       obj1_ids:np.ndarray,
                       theta:astropy.units.Quantity=6 * u.arcmin,
                       n_in_area=700,
                       area=1 * u.deg**2,
                       catalogue_name='lotss',
                       n_searched=None):
    """Expected number of chance associations within a circle of radius theta.

    Every source-catalogue object that was searched has the same probability
    p = n_in_area * pi * theta^2 / area of an unrelated target falling inside theta,
    whether or not it ended up matched. The expected number of chance associations
    is therefore p * n_searched.

    Stores Chance(total=obj1_ids.size, expected=p * n_searched) on
    obj1.matches[catalogue_name], so `total` (matches found) can be compared with
    `expected` (matches expected by chance).

    Parameters
    ----------
    obj1_ids: indices of the matched source-catalogue objects (sets `total` only)
    n_searched: number of source-catalogue objects searched; defaults to len(obj1.df)

    Note: before 2026-09 this multiplied p by obj1_ids.size (the number of matches)
    instead of the number of objects searched. Results from that version are not comparable.
    """
    p = chance_association_probability(
        n=n_in_area,
        theta=theta.to(u.deg),
        area=area
    ).value

    obj1.matches[catalogue_name].set_chance(obj1_ids.size, p * _n_searched(obj1, n_searched))


def chance_association_ab(obj1:Catalog,
                          obj1_ids:np.ndarray,
                          n_in_area=700,
                          area=1 * u.deg**2,
                          catalogue_name='lotss',
                          searched_ids=None,
                          verbose=False):
    """As chance_association, but each searched object contributes its own ellipse area.

    expected = sum over searched objects of n_in_area * pi * a * b / area, skipping objects
    without a size. searched_ids defaults to every row of obj1.df. Before 2026-09 the sum
    ran over the matched objects only (see chance_association).
    """
    a = np.asarray(obj1.semi_major(with_unit=True, to_unit=u.deg).value)
    b = np.asarray(obj1.semi_minor(with_unit=True, to_unit=u.deg).value)
    ids = np.arange(len(a)) if searched_ids is None else np.asarray(searched_ids)
    partial = chance_association_probability_ab(n=n_in_area, a=a[ids] * u.deg, b=b[ids] * u.deg, area=area).value
    expected = float(np.nansum(partial))

    if verbose:
        print(f'{np.isfinite(partial).sum()} searched objects with a size; '
              f'expected chance associations: {expected:.2f}; matched: {obj1_ids.size}')

    obj1.matches[catalogue_name].set_chance(obj1_ids.size, expected)
