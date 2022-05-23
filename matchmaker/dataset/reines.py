import numpy as np
import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)

class Reines(Catalog):
    file_location = DATA_BASE_PATH + 'data/reines/Reines_table3.csv'
    name = 'reines'

    def __init__(self, load_data=False, constrain=True, use_a=True, use_distance=False):
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg), use_distance=use_distance)

        self.set_survey_specific_columns()
        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.survey.s9
        self.cols.measure.luminosity = self.cols.survey.l9
        self.cols.measure.alpha = self.cols.survey.alpha
        self.cols.measure.offset = self.cols.survey.offset

        if load_data:
            self.load_data(constrain=constrain, use_a=use_a)

    def load_data(self, constrain=False, use_a=True):
        df = pd.read_csv(self.file_location, index_col='ID')
        self.df = df

    def set_survey_specific_columns(self):
        self.cols.survey = Column()
        self.cols.survey.id = Column('ID', None, None)
        self.cols.survey.name = Column('Name', None, None)
        self.cols.survey.ra = Column('RA', u.deg, None)
        self.cols.survey.dec = Column('DEC', u.deg, None)
        self.cols.survey.offset = Column('Offset', u.arcsec, 'Offset from optical center of the galaxy in units of arcseconds')
        self.cols.survey.s9 =  Column('S9 GHz', u.uJy, r'Flux density at 9 GHz in units of μJy')
        self.cols.survey.s9_err =  Column('S9 err', u.uJy, r'Flux density at 9 GHz in units of μJy')
        self.cols.survey.l9 = Column('log L9 GHz', u.W * u.Hz**-1, 'log spectral luminosity at 9 GHz in units of W Hz−1')
        self.cols.survey.alpha = Column('alpha', None, r'spectral index (Sν ∝ ν α) determined from flux densities at 9.00 and 10.65 GHz')
        self.cols.survey.alpha_err = Column('alpha err', None, r'spectral index (Sν ∝ ν α) determined from flux densities at 9.00 and 10.65 GHz')
        self.cols.survey.point = Column('Point source', None, 'Designation as a point source (or not), determined from IMFIT')
        self.cols.survey.note = Column('Note', None, 'indication of whether the radio source is in Sample A, Sample B, or is a background source. “AGN” indicates radio sources that are too luminous to be consistent with star formation. “SF” indicates radio sources consistent with star formation')
