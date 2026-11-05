"""Positional-error units: RACS in arcsec; NVSS RA in seconds of time."""
import astropy.units as u
import numpy as np
import pandas as pd
from astropy.table import Table

from matchmaker.crossmatch import crossmatch
from matchmaker.dataset import Catalog, Column
from matchmaker.dataset.nvss import Nvss
from matchmaker.dataset.racs import Racs


def test_racs_errors_in_arcsec():
    r = Racs(load_data=False)
    assert r.cols.ra_err.unit == u.arcsec and r.cols.dec_err.unit == u.arcsec


def _nvss_file(tmp_path, ra, dec, e_ra_s, e_dec_arcsec):
    t = Table({'RAJ2000': ra, 'DEJ2000': dec, 'e_RAJ2000': e_ra_s, 'e_DEJ2000': e_dec_arcsec,
               'S1_4': np.full(len(ra), 5.0), 'e_S1_4': np.full(len(ra), 0.5), 'resFlux': np.zeros(len(ra))})
    p = tmp_path / 'nvss.fits'
    t.write(p)
    return str(p)


def test_nvss_ra_error_converted_from_seconds_of_time(tmp_path):
    n = Nvss(load_data=False)
    n.file_location = _nvss_file(tmp_path, [10.0, 20.0], [0.0, 60.0], [0.3, 0.3], [3.2, 3.2])
    n.load_data()
    ra_err = n.prop_to_unit('ra_err', u.arcsec)
    assert np.allclose(ra_err, [0.3 * 15, 0.3 * 15 * 0.5])       # 15 cos(dec) arcsec per second
    assert np.allclose(n.prop_to_unit('dec_err', u.arcsec), [3.2, 3.2])


def test_nvss_as_source_no_longer_matches_at_degree_scale(tmp_path):
    """With errors read as degrees, a 30" offset was accepted at sep=6" (acceptance ~ degrees)."""
    n = Nvss(load_data=False)
    n.file_location = _nvss_file(tmp_path, [150.0], [30.0], [0.3], [3.2])
    n.load_data()
    target = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    target.name = 'target'
    target.df = pd.DataFrame({'ra': [150.0], 'dec': [30.0 + 30 / 3600]})
    crossmatch(n, target, sep=6, unit_sep=u.arcsec, n_sigma=1.0)
    assert not n.matches['target'].mask[0]
