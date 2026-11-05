"""Survey footprints from source positions (matchmaker.footprint)."""
import astropy.units as u
import numpy as np
import pandas as pd
import pytest

pytest.importorskip('astropy_healpix')

from matchmaker.dataset import Catalog, Column  # noqa: E402
from matchmaker.footprint import Footprint  # noqa: E402


def _patch(n=200_000, ra=(150, 170), dec=(20, 40), seed=0):
    rng = np.random.default_rng(seed)
    sd = rng.uniform(np.sin(np.radians(dec[0])), np.sin(np.radians(dec[1])), n)
    return rng.uniform(*ra, n), np.degrees(np.arcsin(sd))


def test_area_of_a_uniform_patch():
    ra, dec = _patch()
    fp = Footprint.from_positions(ra, dec, nside=64)
    true = np.radians(20) * (np.sin(np.radians(40)) - np.sin(np.radians(20))) * (180 / np.pi) ** 2
    # edge cells count in full, so the area exceeds the truth by roughly perimeter x cell width
    assert true < fp.area.to_value(u.deg**2) < 1.15 * true


def test_contains():
    ra, dec = _patch()
    fp = Footprint.from_positions(ra, dec, nside=64)
    assert list(fp.contains([160, 160, 200], [30, 60, 30])) == [True, False, False]


def test_values_are_the_per_cell_majority():
    ra, dec = _patch()
    res = np.where(dec < 30, 9.0, 6.0)
    fp = Footprint.from_positions(ra, dec, nside=64, values=res)
    assert list(fp.value_at([160, 160, 200], [25, 35, 30])) == [9.0, 6.0, pytest.approx(np.nan, nan_ok=True)]
    total = fp.area_where(6.0) + fp.area_where(9.0)
    assert np.isclose(total.value, fp.area.value)


def test_catalog_filter_footprint_keeps_inside_and_original_index():
    ra, dec = _patch()
    fp = Footprint.from_positions(ra, dec, nside=64)
    cat = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    cat.name = 'clu'
    cat.df = pd.DataFrame({'ra': [160.0, 10.0, 165.0], 'dec': [30.0, 30.0, 25.0]})
    cat.filter_footprint(fp)
    assert list(cat.df['index_df_clu']) == [0, 2]
    assert list(cat.df.index) == [0, 1]
