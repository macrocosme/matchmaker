import astropy.units as u
import numpy as np
import pandas as pd

from matchmaker.crossmatch import Ellipse, crossmatch
from matchmaker.dataset import Catalog, Column

ARCSEC = 1. / 3600.


def make_catalog(name, ra_deg, dec_deg, ra_err_arcsec=None, dec_err_arcsec=None):
    cat = Catalog(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
    cat.name = name
    data = {'ra': ra_deg, 'dec': dec_deg}
    if ra_err_arcsec is not None:
        cat.cols.ra_err = Column('ra_err', u.arcsec)
        cat.cols.dec_err = Column('dec_err', u.arcsec)
        data['ra_err'] = ra_err_arcsec
        data['dec_err'] = dec_err_arcsec
    cat.df = pd.DataFrame(data)
    return cat


class GalaxyCatalog(Catalog):
    """Minimal catalog with an angular extent, for the ellipse filter."""

    def __init__(self, name, ra_deg, dec_deg, major_arcsec, minor_arcsec, pa_deg):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
        self.name = name
        self.cols.major = Column('major', u.arcsec)
        self.cols.minor = Column('minor', u.arcsec)
        self.cols.pa = Column('pa', u.deg)
        self.df = pd.DataFrame({'ra': ra_deg, 'dec': dec_deg,
                                'major': major_arcsec, 'minor': minor_arcsec,
                                'pa': pa_deg})

    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        return self.prop_to_unit('major', to_unit or self.cols.major.unit,
                                 mask=mask, with_unit=with_unit) / 2

    def semi_minor(self, mask=None, with_unit=False, to_unit=None):
        return self.prop_to_unit('minor', to_unit or self.cols.minor.unit,
                                 mask=mask, with_unit=with_unit) / 2


class TestEllipse:
    # Convention: pa measured East of North, so pa=0 puts the major axis
    # along +dec and pa=90 along +ra

    def test_center_is_inside(self):
        e = Ellipse(x=10., y=20., a=1., b=0.5, pa=0.)
        assert e.is_point_inside(10., 20.)

    def test_pa0_major_axis_points_north(self):
        e = Ellipse(x=0., y=0., a=1., b=0.5, pa=0.)
        assert e.is_point_inside(0., 0.99)       # along major (North)
        assert not e.is_point_inside(0.6, 0.)    # beyond b to the East
        assert e.is_point_inside(0.4, 0.)        # within b to the East

    def test_pa90_major_axis_points_east(self):
        e = Ellipse(x=0., y=0., a=2., b=0.5, pa=90.)
        assert e.is_point_inside(1.5, 0.)        # East, within a
        assert not e.is_point_inside(0., 1.5)    # North, beyond b

    def test_rotated_major_axis(self):
        # Point at r=0.9 along the pa=45deg direction (NE) must be inside a=1
        pa = 45.
        e = Ellipse(x=0., y=0., a=1., b=0.1, pa=pa)
        r = 0.9
        x = r * np.sin(np.radians(pa))
        y = r * np.cos(np.radians(pa))
        assert e.is_point_inside(x, y)
        # Same distance perpendicular to the major axis: outside (b=0.1)
        assert not e.is_point_inside(-r * np.cos(np.radians(pa)),
                                     r * np.sin(np.radians(pa)))

    def test_ra_offsets_scaled_by_cos_dec(self):
        # 1-degree circle at dec=60: 1.9deg of RA is only 0.95deg on sky
        e = Ellipse(x=10., y=60., a=1., b=1., pa=0.)
        assert e.is_point_inside(11.9, 60.)
        assert not e.is_point_inside(12.2, 60.)  # 1.1deg on sky

    def test_ra_wraparound(self):
        e = Ellipse(x=359.9, y=0., a=0.2, b=0.2, pa=0.)
        assert e.is_point_inside(0.05, 0.)

    def test_vectorized_input(self):
        e = Ellipse(x=0., y=0., a=1., b=1., pa=0.)
        inside = e.is_point_inside(np.array([0., 2.]), np.array([0., 0.]))
        assert inside.tolist() == [True, False]


class TestCrossmatch:
    def test_simple_match_within_separation(self):
        # source: 3 objects; target: close match for the first only
        source = make_catalog('src', [10.0, 50.0, 200.0], [45.0, -10.0, 30.0])
        target = make_catalog('tgt', [10.0001, 120.0], [45.0001, 60.0])

        crossmatch(source, target, sep=2, unit_sep=u.arcsec)

        mask = source.matches['tgt'].mask
        assert mask.tolist() == [True, False, False]

        match = target.matches['src']
        assert match.idx[0] == 0  # nearest target of source 0
        assert match.filtered_idx.tolist() == [0]
        assert match.mask_idx.tolist() == [0]

    def test_no_match_when_separation_too_small(self):
        source = make_catalog('src', [10.0], [45.0])
        target = make_catalog('tgt', [10.1], [45.1])  # ~510 arcsec away

        crossmatch(source, target, sep=1, unit_sep=u.arcsec)

        assert source.matches['tgt'].mask.tolist() == [False]

    def test_crossmatch_does_not_mutate_target_coordinates(self):
        source = make_catalog('src', [10.0], [45.0])
        target = make_catalog('tgt', [359.9, 0.1], [89.5, -89.5])
        before = target.df.copy()

        crossmatch(source, target, sep=1, unit_sep=u.arcsec)

        pd.testing.assert_frame_equal(target.df, before)

    def test_match_recorded_on_both_catalogs(self):
        source = make_catalog('src', [10.0], [45.0])
        target = make_catalog('tgt', [10.0], [45.0])

        crossmatch(source, target, sep=1, unit_sep=u.arcsec)

        assert 'tgt' in source.matches
        assert 'src' in target.matches
        assert source.matches['tgt'].from_to == 0
        assert target.matches['src'].from_to == 1


class TestErrorAwareRadius:
    def test_positional_errors_grow_acceptance_radius(self):
        # Target 3" away with 2.5" errors per axis: sigma ~3.54",
        # threshold = 0.5" + 3.54" -> matched
        source = make_catalog('src', [10.0], [0.0])
        target = make_catalog('tgt', [10.0 + 3 * ARCSEC], [0.0],
                              ra_err_arcsec=[2.5], dec_err_arcsec=[2.5])

        crossmatch(source, target, sep=0.5, unit_sep=u.arcsec)
        assert source.matches['tgt'].mask.tolist() == [True]

    def test_n_sigma_zero_is_pure_radius_match(self):
        source = make_catalog('src', [10.0], [0.0])
        target = make_catalog('tgt', [10.0 + 3 * ARCSEC], [0.0],
                              ra_err_arcsec=[2.5], dec_err_arcsec=[2.5])

        crossmatch(source, target, sep=0.5, unit_sep=u.arcsec, n_sigma=0)
        assert source.matches['tgt'].mask.tolist() == [False]

    def test_acceptance_is_monotonic_in_distance(self):
        # Regression for the shift-grid holes: with err=2" and sep=0.5",
        # a target at 1" must match if one at 2" does
        for d_arcsec in (0.4, 1.0, 1.6, 2.0):
            source = make_catalog('src', [10.0], [0.0])
            target = make_catalog('tgt', [10.0 + d_arcsec * ARCSEC], [0.0],
                                  ra_err_arcsec=[2.0], dec_err_arcsec=[2.0])
            crossmatch(source, target, sep=0.5, unit_sep=u.arcsec)
            assert source.matches['tgt'].mask.tolist() == [True], f'hole at {d_arcsec}"'

    def test_nan_errors_treated_as_zero(self):
        source = make_catalog('src', [10.0], [0.0])
        target = make_catalog('tgt', [10.0 + 3 * ARCSEC], [0.0],
                              ra_err_arcsec=[np.nan], dec_err_arcsec=[np.nan])

        crossmatch(source, target, sep=0.5, unit_sep=u.arcsec)
        assert source.matches['tgt'].mask.tolist() == [False]


class TestLegacyFlag:
    def test_legacy_reproduces_shift_grid_hole(self):
        # Documented legacy quirk: target at 1" with 2" errors falls in a hole
        # of the union-of-9-circles acceptance region (the new method matches it)
        source = make_catalog('src', [10.0], [0.0])
        target = make_catalog('tgt', [10.0 + 1 * ARCSEC], [0.0],
                              ra_err_arcsec=[2.0], dec_err_arcsec=[2.0])

        crossmatch(source, target, sep=0.5, unit_sep=u.arcsec, legacy=True)
        assert source.matches['tgt'].mask.tolist() == [False]

        crossmatch(source, target, sep=0.5, unit_sep=u.arcsec)
        assert source.matches['tgt'].mask.tolist() == [True]

    def test_legacy_matches_shifted_position(self):
        # Target exactly err away: the shifted pass lands on the source
        source = make_catalog('src', [10.0], [0.0])
        target = make_catalog('tgt', [10.0 + 2 * ARCSEC], [0.0],
                              ra_err_arcsec=[2.0], dec_err_arcsec=[2.0])

        crossmatch(source, target, sep=0.5, unit_sep=u.arcsec, legacy=True)
        assert source.matches['tgt'].mask.tolist() == [True]

    def test_legacy_records_matches_on_both_catalogs(self):
        source = make_catalog('src', [10.0], [45.0])
        target = make_catalog('tgt', [10.0], [45.0])

        crossmatch(source, target, sep=1, unit_sep=u.arcsec, legacy=True)

        assert source.matches['tgt'].mask.tolist() == [True]
        assert target.matches['src'].idx.tolist() == [0]


class TestEllipseFilter:
    # Galaxy at dec=45 with a=60", b=10", major axis East (pa=90)
    def make_galaxy(self):
        return GalaxyCatalog('gal', [10.0], [45.0],
                             major_arcsec=[120.0], minor_arcsec=[20.0],
                             pa_deg=[90.0])

    def test_target_along_major_axis_kept(self):
        galaxy = self.make_galaxy()
        # 30" East on sky (RA offset scaled by cos dec): inside a=60"
        target = make_catalog('tgt', [10.0 + 30 * ARCSEC / np.cos(np.radians(45.))],
                              [45.0])

        crossmatch(galaxy, target, sep=1, unit_sep=u.arcmin, ellipse_filter=True)
        assert galaxy.matches['tgt'].mask.tolist() == [True]

    def test_target_along_minor_axis_rejected(self):
        galaxy = self.make_galaxy()
        # 30" North: outside b=10"
        target = make_catalog('tgt', [10.0], [45.0 + 30 * ARCSEC])

        crossmatch(galaxy, target, sep=1, unit_sep=u.arcmin, ellipse_filter=True)
        assert galaxy.matches['tgt'].mask.tolist() == [False]

    def test_target_positional_error_inflates_ellipse(self):
        galaxy = self.make_galaxy()
        # 30" North with 25" positional error: 30 <= 10 + 25 -> kept
        target = make_catalog('tgt', [10.0], [45.0 + 30 * ARCSEC],
                              ra_err_arcsec=[0.0], dec_err_arcsec=[25.0])

        crossmatch(galaxy, target, sep=1, unit_sep=u.arcmin, ellipse_filter=True)
        assert galaxy.matches['tgt'].mask.tolist() == [True]
