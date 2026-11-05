"""chance_association counts every object searched, not only the matches."""
import math

import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from matchmaker.dataset import Catalog, Column
from matchmaker.stats import chance_association, chance_association_ab


def _catalog(n):
    cat = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    cat.name = 'src'
    cat.df = pd.DataFrame({'ra': np.zeros(n), 'dec': np.zeros(n)})
    cat.add_match('lotss')
    return cat


def test_chance_association_scales_with_searched_not_matched():
    cat = _catalog(1000)
    chance_association(cat, np.array([0, 1, 2]), theta=6 * u.arcsec, n_in_area=700, area=1 * u.deg**2)
    p = 700 * math.pi * (6 / 3600) ** 2
    assert cat.matches['lotss'].chance.total == 3
    assert np.isclose(cat.matches['lotss'].chance.expected, 1000 * p)


def test_chance_association_dr2_example():
    cat = _catalog(18159)
    chance_association(cat, np.arange(708), theta=6 * u.arcsec, n_in_area=2051534, area=5635 * u.deg**2)
    assert round(cat.matches['lotss'].chance.expected) == 58


def test_chance_association_needs_n_searched_without_dataframe():
    cat = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    cat.add_match('lotss')
    with pytest.raises(ValueError):
        chance_association(cat, np.array([0]))
    chance_association(cat, np.array([0]), n_searched=10, theta=6 * u.arcsec, n_in_area=700, area=1 * u.deg**2)
    assert np.isclose(cat.matches['lotss'].chance.expected, 10 * 700 * math.pi * (6 / 3600) ** 2)


def test_chance_association_ab_sums_over_all_searched():
    class Ell(Catalog):
        def semi_major(self, mask=None, with_unit=False, to_unit=None):
            return (self.df['a'].values * u.arcsec).to(to_unit)

        def semi_minor(self, mask=None, with_unit=False, to_unit=None):
            return (self.df['b'].values * u.arcsec).to(to_unit)

    cat = Ell(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    cat.df = pd.DataFrame({'ra': [0, 0, 0], 'dec': [0, 0, 0], 'a': [6.0, 6.0, np.nan], 'b': [3.0, 3.0, 1.0]})
    cat.add_match('lotss')
    chance_association_ab(cat, np.array([0]), n_in_area=700, area=1 * u.deg**2)
    per = 700 * math.pi * (6 / 3600) * (3 / 3600)
    assert np.isclose(cat.matches['lotss'].chance.expected, 2 * per)     # NaN-size object skipped
    assert cat.matches['lotss'].chance.total == 1


# --- CLU, WISE, sources -------------------------------------------------------------------
