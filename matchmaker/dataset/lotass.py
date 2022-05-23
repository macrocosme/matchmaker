import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)


class Lotass(Catalog):
    file_location = DATA_BASE_PATH + 'data/lotass/lotass_tan_combined.csv'
    name = 'lotass'

    def __init__(self, load_data=False, constrain=True):
        super().__init__(ra=Column('RA_d', u.deg), dec=Column('DEC_d', u.deg))
        self.set_survey_specific_columns()

        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.all.Smean_128_167

        if load_data:
            self.load_data()

    def load_data(self):
        df = pd.read_csv(self.file_location)

        self.df = df

    def set_survey_specific_columns(self):
        self.cols.all = Column()
        self.cols.all.psr_j2000 = Column('PSR (J2000)', None, None)
        self.cols.all.psr_disc = Column('PSR (disc.)', None, None)
        self.cols.all.ra_d = Column('RA_d', u.deg, None)
        self.cols.all.dec_d = Column('DEC_d', u.deg, None)
        self.cols.all.ra = Column('RA (J2000)', None, None)
        self.cols.all.dec = Column('DEC (J2000)', None, None)
        self.cols.all.Epoch = Column('Epoch (MJD)', None, None)
        self.cols.all.P = Column('P (s)', u.s, None)
        self.cols.all.Pdot = Column('Pdot (10^-15)', 10**-15, None)
        self.cols.all.DM = Column('DM (pc cm-3)', u.parsec*u.cm**-3, None)
        self.cols.all.Ntoa = Column('Ntoa', None, None)
        self.cols.all.TRES = Column('TRES (us)', u.us, None)
        self.cols.all.X2red = Column('X2red', None, None)
        self.cols.all.Smean_128_167 = Column('Smean_128_167', None, None)
        self.cols.all.Smean_128_167_error = Column('S_error(mean_128_167)', None, None)
        self.cols.all.S_128mhz = Column('S(128 MHz)', None, 'Flux density at 128 MHz')
        self.cols.all.S_128mhz_error = Column('S_error (128 MHz)', None, 'Flux density uncertainty at 128 MHz')
        self.cols.all.S_167mhz = Column('S(167 MHz)', None, 'Flux density at 167 MHz')
        self.cols.all.S_167mhz_error = Column('S_error(167 MHz)', None, None)
        self.cols.all.S_334mhz = Column('S(334 MHz)', None, 'Flux density at 334 MHz')
        self.cols.all.S_334mhz_error = Column('S_error(334 MHz)', None, None)
        self.cols.all.S_1532mhz = Column('S(1532 MHz)', None, None)
        self.cols.all.S_1532mhz_error = Column('S_error(1532 MHz)', None, None)
        self.cols.all.alpha = Column('alpha', None, None)
        self.cols.all.alpha_error = Column('alpha_error', None, None)
        self.cols.all.Span = Column('Span (months)', None, None)
        self.cols.all.nobs_lofar = Column('nobs_lofar', None, None)
        self.cols.all.nobs_lovell = Column('nobs_lovell', None, None)
        self.cols.all.Obs_length_lofar = Column('Obs_length_lofar (min)', None, None)
        self.cols.all.Obs_length_lovell = Column('Obs_length_lovell (min)', None, None)
