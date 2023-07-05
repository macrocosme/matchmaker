import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe


class Tgss(Catalog):
    file_location = DATA_BASE_PATH + 'TGSS/TGSSADR1_7sigma_catalog.fits'
    name = 'tgss'


    def __init__(self, load_data=False, constrain=False):
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg))
        self.cols.ra_err = Column('E_RA', u.deg)
        self.cols.dec_err = Column('E_DEC', u.deg)
        self.cols.major = Column('Maj', u.arcsec)
        self.cols.minor = Column('Min', u.arcsec)
        self.cols.pa = Column('PA', u.deg)
        # self.set_survey_specific_columns()
        self.cols.measure = Column()
        self.cols.measure.total_flux = Column('Total_flux', u.mJy)
        self.cols.measure.total_flux_err = Column('E_Total_flux', u.mJy)
        self.cols.measure.peak_flux = Column('Peak_flux', u.mJy * u.beam**-1)
        self.cols.measure.peak_flux_err = Column('E_Peak_flux', u.mJy * u.beam**-1)
        self.cols.measure.s_code = Column('Source_code', None)

        # self.beam = 2.5 * u.arcsec
        # self.frequency_band = 'L'
        self.central_frequency = 150 * u.MHz

        # self.area = 10575 * u.deg**2

        if load_data:
            self.load_data(constrain=constrain)

    def load_data(self, constrain=False):
        df = load_fits_as_dataframe(self.file_location)

        # if constrain:
        #     df = df.loc[df[self.cols.measure.total_flux.label] >= 1.5]

        self.df = df

    # Should simply fetch the same code from the Lotss class as it is identical.
    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        a = self.prop_to_unit('major', self.cols.major.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                a = (a * self.cols.major.unit).to(to_unit).value
            else:
                a = a.to(to_unit)

        return a

    # idem
    def semi_minor(self, mask=None, with_unit=False, to_unit=None):
        b = self.prop_to_unit('minor', self.cols.minor.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                b = (b * self.cols.minor.unit).to(to_unit).value
            else:
                b = b.to(to_unit)

        return b

    # def set_survey_specific_columns(self):
    #     self.cols.all = Column()
    #     self.cols.all.ra = Column('RAJ2000', u.deg, None)
    #     self.cols.all.dec = Column('DEJ2000', u.deg, None)
    #
    #     self.cols.measure = Column()
    #     self.cols.measure.fpeak = Column('FPEAK', u.mJy, None)
    #     self.cols.measure.fint = Column('FINT', u.mJy, None)
    #     self.cols.measure.rms = Column('RMS', u.mJy, None)
    #     self.cols.measure.major = Column('MAJOR', u.arcsec, None)
    #     self.cols.measure.minor = Column('MINOR', u.arcsec, None)
    #     self.cols.measure.posang = Column('POSANG', u.deg, None)
    #     self.cols.measure.fitted_major = Column('FITTED_MAJOR', u.arcsec, None)
    #     self.cols.measure.fitted_minor = Column('FITTED_MINOR', u.arcsec, None)
    #     self.cols.measure.fitted_posang = Column('FITTED_POSANG', u.deg, None)
    #
    #     self.cols.source = Column()
    #     self.cols.source.fldname = Column('FLDNAME', None, None)
    #     self.cols.source.nsdss = Column('NSDSS', None, None)
    #
    #     self.cols.crossmatch = Column()
    #     self.cols.crossmatch.sdss_sep = Column('SDSS_SEP', u.arcsec, None)
    #     self.cols.crossmatch.sdss_mag = Column('SDSS_MAG', u.mag, None)
    #     self.cols.crossmatch.sdss_class = Column('SDSS_CLASS', None, None)
    #     self.cols.crossmatch.ntmass = Column('NTMASS', None, None)
    #     self.cols.crossmatch.tmass_sep = Column('TMASS_SEP', u.arcsec, None)
    #     self.cols.crossmatch.tmass_mag = Column('TMASS_MAG', u.mag, None)
    #
    #     self.cols.all.year = Column('YEAR', None, None)
    #     self.cols.all.mjd = Column('MJD', None, None)
    #     self.cols.all.mjdrms = Column('MJDRMS', None, None)
    #     self.cols.all.mjdstart = Column('MJDSTART', None, None)
    #     self.cols.all.mjdstop = Column('MJDSTOP', None, None)
