"""Chance-association Monte Carlo (monte_carlo.mc)."""
import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from monte_carlo import mc


def test_offsets_are_true_angles():
    rng = np.random.default_rng(0)
    dec = np.full(20000, 60.0)
    d_ra, d_dec = mc.random_offsets(len(dec), 10.0, dec, rng)
    on_sky_ra = d_ra * np.cos(np.radians(dec))
    assert np.abs(on_sky_ra).max() <= 10 / 60 + 1e-9        # true angle respects the box
    assert np.abs(d_ra).max() > 10 / 60                      # coordinate offset is larger at dec 60
    assert np.abs(d_dec).max() <= 10 / 60 + 1e-9


def test_shift_wraps_ra_and_reflects_at_the_pole():
    ra, dec = mc.shift_positions(np.array([359.9, 10.0]), np.array([0.0, 89.9]),
                                 np.array([0.2, 0.0]), np.array([0.0, 0.2]))
    assert ra[0] == pytest.approx(0.1)
    assert dec[1] == pytest.approx(89.9) and ra[1] == pytest.approx(190.0)


def _setup(n_gal=500, n_radio=20000, seed=1):
    rng = np.random.default_rng(seed)
    g_ra, g_dec = rng.uniform(150, 152, n_gal), rng.uniform(30, 32, n_gal)
    r_ra, r_dec = rng.uniform(149, 153, n_radio), rng.uniform(29, 33, n_radio)
    radio = SkyCoord(r_ra * u.deg, r_dec * u.deg)
    z = rng.uniform(0.01, 0.04, n_gal)
    sfr = 10 ** rng.uniform(-2, -0.5, n_gal)
    flux = 10 ** rng.uniform(0, 1.5, n_radio)
    return dict(galaxy_ra=g_ra, galaxy_dec=g_dec, galaxy_z=z, galaxy_sfr=sfr,
                radio_coords=radio, radio_flux=flux)


def test_more_matches_with_a_bigger_radius():
    kw = _setup()
    small = mc.run(**kw, radius_arcsec=6.0, n_realisations=20, rng=np.random.default_rng(0))
    big = mc.run(**kw, radius_arcsec=30.0, n_realisations=20, rng=np.random.default_rng(0))
    assert big['n_matches'].mean() > 10 * small['n_matches'].mean()   # area scales as radius^2
    assert (big['n_outliers'] <= big['n_matches']).all()


def test_per_source_radius_is_used():
    kw = _setup()
    n = len(kw['radio_coords'])
    rng = np.random.default_rng(3)
    wide = mc.run(**kw, radius_arcsec=np.full(n, 20.0), n_realisations=15, rng=np.random.default_rng(7))
    mixed = mc.run(**kw, radius_arcsec=rng.choice([6.0, 20.0], n), n_realisations=15,
                   rng=np.random.default_rng(7))
    assert mixed['n_matches'].mean() < wide['n_matches'].mean()


def test_positional_error_widens_the_acceptance():
    kw = _setup()
    n = len(kw['radio_coords'])
    a = mc.run(**kw, radius_arcsec=6.0, n_realisations=15, rng=np.random.default_rng(5))
    b = mc.run(**kw, radius_arcsec=6.0, radio_sigma_arcsec=np.full(n, 10.0), n_sigma_position=1.0,
               n_realisations=15, rng=np.random.default_rng(5))
    assert b['n_matches'].sum() > a['n_matches'].sum()


def test_reproducible_and_region_counts_add_up():
    kw = _setup()
    region = np.where(kw['galaxy_dec'] < 31, '9', '6')
    r1 = mc.run(**kw, radius_arcsec=30.0, region=region, n_realisations=10, rng=np.random.default_rng(2))
    r2 = mc.run(**kw, radius_arcsec=30.0, region=region, n_realisations=10, rng=np.random.default_rng(2))
    assert np.array_equal(r1['n_matches'], r2['n_matches'])
    assert np.array_equal(r1['n_matches_6'] + r1['n_matches_9'], r1['n_matches'])
    assert np.array_equal(r1['n_outliers_6'] + r1['n_outliers_9'], r1['n_outliers'])


def test_footprint_drops_galaxies_that_shift_out():
    pytest.importorskip('astropy_healpix')
    from matchmaker.footprint import Footprint
    kw = _setup()
    fp = Footprint.from_positions(np.random.default_rng(0).uniform(150, 152, 200000),
                                 np.random.default_rng(1).uniform(30, 32, 200000), nside=64)
    r = mc.run(**kw, radius_arcsec=6.0, footprint=fp, shift_arcmin=30.0, n_realisations=10,
               rng=np.random.default_rng(0))
    assert (r['n_galaxies'] < len(kw['galaxy_ra'])).any() and (r['n_galaxies'] > 0).all()


def test_p_value():
    assert mc.p_value([1, 2, 3, 10], 3) == 0.5 and mc.p_value([1, 2], 99) == 0.0
