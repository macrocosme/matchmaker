"""RACS-low loader."""
import astropy.units as u

from matchmaker.dataset.racs import Racs


def test_racs_is_racs_low_frequency():
    assert Racs(load_data=False).central_frequency == 887.5 * u.MHz


def test_racs_flux_columns_match_racs_low_dr1_file():
    r = Racs(load_data=False)
    assert r.cols.measure.total_flux.label == 'total_flux_source'
    assert r.cols.measure.total_flux_err.label == 'e_total_flux_source'


def test_racs_rejects_a_file_that_is_not_racs_low(tmp_path):
    """RACS-mid reaches declination +49 and holds ~3.1 million sources at 1367.5 MHz."""
    import numpy as np
    import pandas as pd
    import pytest

    r = Racs(load_data=False)
    r.df = pd.DataFrame({'ra': [10.0, 20.0], 'dec': [45.0, 48.0], 'e_ra': [1.0, 1.0],
                         'total_flux_source': [1.0, 2.0]})
    r.n_source = 3_100_000
    with pytest.raises(ValueError, match='RACS-low DR1'):
        r.verify_release()
    assert np.isclose(r.central_frequency.value, 887.5)


def test_racs_accepts_the_published_file():
    import pandas as pd

    r = Racs(load_data=False)
    r.df = pd.DataFrame({'ra': [10.0], 'dec': [29.1], 'e_ra': [1.0], 'total_flux_source': [1.0]})
    r.n_source = 2123638
    r.verify_release()          # must not raise


def test_racs_says_where_the_file_is_when_it_is_missing():
    import pytest

    r = Racs(load_data=False)
    r.file_location = '/nonexistent/racs.fits'
    with pytest.raises(FileNotFoundError, match='reproduction_package'):
        r.load_data()
