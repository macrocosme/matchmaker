"""NED Local Volume Sample: galaxies within 1000 Mpc with photometry and derived properties.

An alternative parent sample to CLU, which has not been updated since 2019. NED-LVS is
regenerated a few times a year and includes redshifts from DESI and other recent surveys, so it
is more complete, and it carries angular diameters, which CLU effectively does not.

Star formation rates come in two flavours: `SFR_W4` from the 22 micron luminosity alone, and
`SFR_hybrid` from far-ultraviolet plus 22 micron. The hybrid rate is the closer analogue of CLU's
`sfr_fuv` but exists for only about a third of the catalogue, since it needs GALEX coverage;
choose with `sfr_column` and be aware that mixing the two across a sample introduces a systematic
that follows GALEX coverage.
"""
import astropy.units as u
import numpy as np

from ..utils import load_fits_as_dataframe
from . import DATA_BASE_PATH, Catalog, Column

GALAXY_TYPES = ('G', 'GPair', 'GTrpl')


class NedLvs(Catalog):
    file_location = DATA_BASE_PATH + 'nedlvs/NEDLVS_20260424.fits'
    name = 'nedlvs'

    def __init__(self, load_data=False, constrain=True, use_distance=True,
                 sfr_column='SFR_hybrid', constrain_mass=True, constrain_sfr=True,
                 galaxies_only=True):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg),
                         use_distance=use_distance)
        self.cols.z = Column('z')
        self.cols.mstar = Column('Mstar', u.solMass)
        self.cols.distance = Column('DistMpc', u.Mpc)
        self.cols.a = Column('Diam', u.arcsec, 'angular diameter')
        self.cols.b2a = Column('Diam_ba')
        self.cols.pa = Column('Diam_pa', u.deg)

        self.sfr_column = sfr_column
        self.set_survey_specific_columns()
        self.cols.sfr = getattr(self.cols.measure, sfr_column.lower())

        if load_data:
            self.load_data(constrain=constrain, constrain_mass=constrain_mass,
                           constrain_sfr=constrain_sfr, galaxies_only=galaxies_only)

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self, constrain=True, constrain_mass=True, constrain_sfr=True,
                  galaxies_only=True, mstar_min=1e7, mstar_max=3e9):
        self.df = load_fits_as_dataframe(self.file_location)
        if galaxies_only:
            self.df = self.df[self.df['objtype'].isin(GALAXY_TYPES)]
        if constrain:
            keep = np.ones(len(self.df), dtype=bool)
            if constrain_mass:
                mstar = self.df['Mstar'].values
                keep &= np.isfinite(mstar) & (mstar > mstar_min) & (mstar <= mstar_max)
            if constrain_sfr:
                sfr = self.df[self.sfr_column].values
                keep &= np.isfinite(sfr) & (sfr > 0)
            self.df = self.df[keep]
        self.df = self.df.reset_index(drop=True)

    def set_survey_specific_columns(self):
        self.cols.ra = Column('ra', u.deg, 'RA J2000 in degrees')
        self.cols.dec = Column('dec', u.deg, 'Dec J2000 in degrees')
        self.cols.measure = Column()
        self.cols.measure.sfr_hybrid = Column('SFR_hybrid', u.Msun * u.yr**-1,
                                              'FUV + 22 micron star formation rate')
        self.cols.measure.sfr_hybrid_err = Column('SFR_hybrid_unc', u.Msun * u.yr**-1)
        self.cols.measure.sfr_w4 = Column('SFR_W4', u.Msun * u.yr**-1,
                                          '22 micron star formation rate')
        self.cols.measure.sfr_w4_err = Column('SFR_W4_unc', u.Msun * u.yr**-1)
        self.cols.measure.mstar = Column('Mstar', u.solMass, 'stellar mass from W1')
        self.cols.measure.mstar_err = Column('Mstar_unc', u.solMass)
        for band in ('FUV', 'NUV', 'J', 'H', 'Ks', 'W1', 'W2', 'W3', 'W4'):
            setattr(self.cols.measure, f'm_{band.lower()}', Column(f'm_{band}', u.mag))
            setattr(self.cols.measure, f'lum_{band.lower()}', Column(f'Lum_{band}', u.solLum))

    def angular_diameter(self, mask=None, to_unit=None):
        """Angular diameter in arcsec, or `to_unit`. Unlike CLU's, this is populated and sane."""
        values = self.df[self.cols.a.label]
        values = values if mask is None else values.iloc[mask]
        out = np.asarray(values, dtype=float) * self.cols.a.unit
        return out.to(to_unit).value if to_unit else out.value
