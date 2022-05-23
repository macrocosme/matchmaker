import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_apj_table
from astropy.coordinates import SkyCoord

class Vlba(Catalog):
    file_location = DATA_BASE_PATH + 'data/vlba_5ghz/vlba_5ghz.txt'
    name = 'vlba'


    def __init__(self, load_data=False, constrain=True):
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg))
        self.frequency_band = None
        self.central_frequency = 5 * u.MHz

        self.set_survey_specific_columns()
        self.cols.measure = Column()

        self.cols.major = self.cols.survey.r
        self.cols.minor = self.cols.survey.r
        self.cols.pa = self.cols.survey.pa_jet

        self.cols.measure.flux_8_5_GHz = self.cols.survey.f8_5
        self.cols.measure.total_flux = self.cols.survey.f5
        self.cols.measure.peak_flux = self.cols.survey.p5

        if load_data:
            self.load_data()

    def load_data(self):
        df = load_apj_table(self.file_location)

        # need to convert ra and dec input to deg.
        df['RA'] = df.apply(lambda row: SkyCoord("{} {} {} {}{} {} {}".format(row['RAh'], row['RAm'], row['RAs'],
                                                               row['DE-'], row['DEd'], row['DEm'], row['DEs']), unit=(u.hourangle, u.deg)).ra.deg, axis=1)
        
        df['DEC'] = df.apply(lambda row: SkyCoord("{} {} {} {}{} {} {}".format(row['RAh'], row['RAm'], row['RAs'],
                                                               row['DE-'], row['DEd'], row['DEm'], row['DEs']), unit=(u.hourangle, u.deg)).dec.deg, axis=1)

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
        self.cols.survey = Column()
        self.cols.survey.num = Column('Num', None, 'VIPS source number')
        self.cols.survey.name = Column('Name', None, 'VIPS source name')
        self.cols.survey.ra_h = Column('RAh', None, None)
        self.cols.survey.ra_m = Column('RAm', None, None)
        self.cols.survey.ra_s = Column('RAs', None, None)
        self.cols.survey.dec_sign = Column('DE-', None, None)
        self.cols.survey.dec_d = Column('DEd', None, None)
        self.cols.survey.dec_m = Column('DEm', None, None)
        self.cols.survey.dec_s = Column('DEs', None, None)
        self.cols.survey.obs_y = Column('Obs.Y', None, 'UT year of the observation')
        self.cols.survey.obs_m = Column('Obs.M', None, 'UT month of the observation')
        self.cols.survey.obs_d = Column('Obs.D', None, 'UT day of the observation')
        self.cols.survey.f8_5 = Column('F8.5', u.mJy, 'CLASS survey 8.5 GHz flux density')
        self.cols.survey.f5 = Column('F5', u.mJy, 'Total cleaned 5 GHz flux density')
        self.cols.survey.p5 = Column('P5', u.mJy, 'Peak 5 GHz flux density')
        self.cols.survey.rms5 = Column('rms5', u.mJy, 'The rms noise in the 5 GHz VLBA map')
        self.cols.survey.ngc = Column('NGC', None, 'Number of dominant Gaussian components fit')
        self.cols.survey.ta = Column('Ta', None, 'Automatic source type')
        self.cols.survey.te = Column('Te', None, 'Visual inspection source type')
        self.cols.survey.r = Column('R', 'u.mas', 'Mean radius dominant Gaussian components')
        self.cols.survey.d_max = Column('Dmax', u.mas, 'Dominant Gaussian components maximum separation')
        self.cols.survey.pa_jet = Column('PAjet', u.deg, 'Jet position angle')

