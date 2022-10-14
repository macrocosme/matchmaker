import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe


class Vitale(Catalog):
    file_location = DATA_BASE_PATH + 'data/first/vitale_et_al_2015_effelsberg.fits'
    name = 'vitale'


    def __init__(self, load_data=False, constrain=False):
        super().__init__(ra=Column('RAJ2000', u.deg), dec=Column('DEJ2000', u.deg))
        self.cols.ra_err = Column()
        self.cols.dec_err = Column()
        self.cols.major = Column('MAJOR', u.arcsec)
        self.cols.minor = Column('MINOR', u.arcsec)
        self.cols.pa = Column('POSANG', u.deg)
        self.set_survey_specific_columns()
        self.cols.measure.total_flux = self.cols.measure.f1_4

        # self.beam = 2.5 * u.arcsec
        # self.frequency_band = 'L'
        # self.central_frequency = 1400 * u.MHz

        # self.area = 10575 * u.deg**2

        if load_data:
            self.load_data(constrain=constrain)

    def load_data(self, constrain=False):
        df = load_fits_as_dataframe(self.file_location)

        # Some RA and DEC uncertainties are marked as -99
        # which makes the code crash later, so I replace it with 0 for now.
        # The following comparison may break at some point...
        # df.loc[df[self.cols.ra_err.label] == -99] = 0.
        # df.loc[df[self.cols.dec_err.label] == -99] = 0.
        if constrain:
            df = df.loc[df[self.cols.measure.total_flux.label] >= 1.5]

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

    def set_survey_specific_columns(self):
        self.cols.all = Column()
        self.cols.all.ra = Column('RAJ2000', u.deg, None)
        self.cols.all.dec = Column('DEJ2000', u.deg, None)
        self.cols.all.recno = Column('recno', None, 'Record number assigned by the VizieR team. Should Not be used for identification.')

        self.cols.measure = Column()
        self.cols.measure.z = Column('z', None, '[0/0.4] Redshift')
        self.cols.measure.f1_4 = Column('F1_4', u.mJy, 'Integrated flux density at 1.4GHz')
        self.cols.measure.f10_45 = Column('f10_45', u.Jy, 'Flux density at 10.45GHz')
        self.cols.measure.f10_45_err = Column('e_f10_45', u.Jy, 'rms uncertainty on F10.45GHz')
        self.cols.measure.f4_85 = Column('f4_85', u.Jy, 'Flux density at 4.85GHz')
        self.cols.measure.f4_85_err = Column('e_f4_85', u.Jy, 'rms uncertainty on F4.85GHz')

        self.cols.measure.mtype = Column('MType', None, 'Optical morphology from SDSS images (1)')
        self.cols.measure.rmorph = Column('RMorph', None, 'Radio morphology from FIRST images (2)')
        self.cols.measure.activity = Column('Activity', None, 'Activity type from NED (3)')
        self.cols.measure.niid = Column('NIId', None, '[LCS FS] Spectral classification according to the [NII]-based diagnostic diagram (4)')

