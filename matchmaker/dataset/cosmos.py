import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe

class Gsnr(Catalog):
    file_location = DATA_BASE_PATH + 'cosmos2015/cosmos2015v21_small.fits'
    name = 'galactic_snr'

    def __init__(self, load_data=False):
        super().__init__(ra=Column('ALPHA_J2000', u.deg), dec=Column('DELTA_J2000', u.deg), use_distance=use_distance)

        self.set_survey_specific_columns()

        if load_data:
            self.load_data()

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self):
        df = load_fits_as_dataframe(self.file_location)
        self.df = df

    def set_survey_specific_columns(self):
        self.cols.ALPHA_J2000 = Column('ALPHA_J2000', u.deg, None)
        self.cols.DELTA_J2000 = Column('DELTA_J2000', u.deg, None)
        self.cols.ID_XMM = Column('ID_XMM', None, None)
        self.cols.FLUX_XMM_05_2 = Column('FLUX_XMM_05_2', u.erg * u.cm**-2 * u.s**-1, None)
        self.cols.FLUX_XMM_2_10 = Column('FLUX_XMM_2_10', u.erg * u.cm**-2 * u.s**-1, None)
        self.cols.FLUX_XMM_5_10 = Column('FLUX_XMM_5_10', u.erg * u.cm**-2 * u.s**-1, None)
        self.cols.HARDNESS_XMM = Column('HARDNESS_XMM', None, None)
        self.cols.ID_CHANDRA09 = Column('ID_CHANDRA09', None, None)
        self.cols.FLUX_CHANDRA_05_2 = Column('FLUX_CHANDRA_05_2', u.erg * u.cm**-2 * u.s**-1, None)
        self.cols.FLUX_CHANDRA_2_10 = Column('FLUX_CHANDRA_2_10', u.erg * u.cm**-2 * u.s**-1, None)
        self.cols.FLUX_CHANDRA_05_10 = Column('FLUX_CHANDRA_05_10', u.erg * u.cm**-2 * u.s**-1, None)
        self.cols.PHOTOZ = Column('PHOTOZ', None, None)
        self.cols.TYPE = Column('TYPE', None, None)
        self.cols.AGE = Column('AGE', None, None)
        self.cols.M_NUV = Column('M_NUV', u.mag, None)
        self.cols.M_U = Column('M_U', u.mag, None)
        self.cols.M_B = Column('M_B', u.mag, None)
        self.cols.M_V = Column('M_V', u.mag, None)
        self.cols.M_R = Column('M_R', u.mag, None)
        self.cols.M_I = Column('M_I', u.mag, None)
        self.cols.M_Z = Column('M_Z', u.mag, None)
        self.cols.M_Y = Column('M_Y', u.mag, None)
        self.cols.M_J = Column('M_J', u.mag, None)
        self.cols.M_H = Column('M_H', u.mag, None)
        self.cols.M_K = Column('M_K', u.mag, None)
        self.cols.MASS_MED = Column('MASS_MED', None, 'log Stellar mass from BC03 best-fit template. median of the PDF')
        self.cols.MASS_MED_MIN68 = Column('MASS_MED_MIN68', None, 'lower limit, 68% confidence level')
        self.cols.MASS_MED_MAX68 = Column('MASS_MED_MAX68', None, 'upper limit, 68% confidence level')
        self.cols.MASS_BEST = Column('MASS_BEST', None, None)
        self.cols.SFR_MED = Column('SFR_MED', None, 'log SFR from BC03 best-fit template. median of the PDF')
        self.cols.SFR_MED_MIN68 = Column('SFR_MED_MIN68', None, 'lower limit, 68% confidence level')
        self.cols.SFR_MED_MAX68 = Column('SFR_MED_MAX68', None, 'upper limit, 68% confidence level')
        self.cols.SFR_BEST = Column('SFR_BEST', None, 'log SFR from BC03 best-fit template. Taken at the minimum chi2')
        self.cols.SSFR_MED = Column('SSFR_MED', None, None)
        self.cols.SSFR_MED_MIN68 = Column('SSFR_MED_MIN68', None, None)
        self.cols.SSFR_MED_MAX68 = Column('SSFR_MED_MAX68', None, None)
        self.cols.SSFR_BEST = Column('SSFR_BEST', None, None)

