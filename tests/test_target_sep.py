"""crossmatch with a per-target acceptance radius (target_sep)."""
import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from matchmaker.crossmatch import crossmatch
from matchmaker.dataset import Catalog, Column


def _cat(name, ra, dec):
    c = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    c.name = name
    c.df = pd.DataFrame({'ra': ra, 'dec': dec})
    return c


def test_radius_follows_the_matched_target():
    # two sources, each 7.5" from its own target; target 0 has a 6" beam, target 1 a 9" beam
    off = 7.5 / 3600
    src = _cat('src', [10.0, 50.0], [0.0, 0.0])
    tgt = _cat('tgt', [10.0, 50.0], [off, off])
    crossmatch(src, tgt, unit_sep=u.arcsec, n_sigma=0, target_sep=np.array([6.0, 9.0]))
    assert list(src.matches['tgt'].mask) == [False, True]


def test_target_sep_must_match_target_rows():
    src, tgt = _cat('src', [10.0], [0.0]), _cat('tgt', [10.0, 20.0], [0.0, 0.0])
    with pytest.raises(ValueError):
        crossmatch(src, tgt, unit_sep=u.arcsec, target_sep=np.array([6.0]))


def test_target_sep_not_with_legacy():
    src, tgt = _cat('src', [10.0], [0.0]), _cat('tgt', [10.0], [0.0])
    with pytest.raises(ValueError):
        crossmatch(src, tgt, unit_sep=u.arcsec, target_sep=np.array([6.0]), legacy=True)
