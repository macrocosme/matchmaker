import math

import astropy.units as u
import numpy as np

from matchmaker.dataset import Catalog, Column
from matchmaker.stats import (
    chance_association,
    chance_association_probability,
    chance_association_probability_ab,
    expected_background_source,
)


def test_chance_association_probability_circle():
    # n * pi * theta^2 / area
    assert np.isclose(chance_association_probability(700, 0.1, 1), 700 * math.pi * 0.01)


def test_chance_association_probability_ellipse():
    assert np.isclose(chance_association_probability_ab(700, 0.1, 0.05, 1), 700 * math.pi * 0.005)


def test_expected_background_source():
    assert expected_background_source(0.1, 50) == 5.0


def test_chance_association_sets_chance_on_match():
    cat = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    cat.name = 'src'
    cat.add_match('lotss')

    ids = np.array([0, 1, 2])
    chance_association(cat, ids, theta=6 * u.arcmin, n_in_area=700, area=1 * u.deg**2,
                       catalogue_name='lotss', n_searched=3)

    chance = cat.matches['lotss'].chance
    per_source = 700 * math.pi * 0.1**2  # 6 arcmin = 0.1 deg
    assert chance.total == 3
    assert np.isclose(chance.expected, 3 * per_source)
