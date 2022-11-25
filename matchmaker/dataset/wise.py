import numpy as np
import pandas as pd
from ..utils import load_fits_as_dataframe
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..model import wise

class Wise(Catalog):
    file_location = DATA_BASE_PATH + 'data/wise/allwise.fits'
    name = 'wise'

    def __init__(self, load_data=False, constrain=True, use_distance=False):
        super().__init__(ra=Column('RAJ2000', u.deg), dec=Column('DEJ2000', u.deg), use_distance=use_distance)
        self.set_survey_specific_columns()

        if load_data:
            self.load_data(constrain=constrain)

    @property
    def n_source(self):
        return len(self.df)

    @property
    def w1w2(self):
        return self.df[self.cols.source.w1.label] - self.df[self.cols.source.w2.label]

    @property
    def w2w3(self):
        return self.df[self.cols.source.w2.label] - self.df[self.cols.source.w3.label]

    def load_data(self, constrain=False):
        df = load_fits_as_dataframe(self.file_location)
        df.drop(['recno'], axis=1, inplace=True)
        self.df = df

    def non_agn_mask(self, mask=None, invert_logic=False):
        from matplotlib.path import Path
        w1w2 = self.w1w2
        w2w3 = self.w2w3
        if type(w2w3) in [np.float, np.float32, np.float64]:
            return ~Path(wise.AGN_color_color_box).contains_point([w2w3, w1w2]) if not invert_logic \
                else Path(wise.AGN_color_color_box).contains_point([w2w3, w1w2])
        else:
            return ~Path(wise.AGN_color_color_box).contains_points(np.array([w2w3, w1w2]).T) if not invert_logic \
                else Path(wise.AGN_color_color_box).contains_points(np.array([w2w3, w1w2]).T)

    def set_survey_specific_columns(self):
        self.cols.source = Column()
        self.cols.source.ra = Column('RAJ2000', u.deg, "RA J2000")
        self.cols.source.dec = Column('DEJ2000', u.deg, "Dec J2000")
        self.cols.source.w1 = Column('W1mag', u.mag, "Uncertainty on magnitude in W1 filter.")
        self.cols.source.w3_err = Column('e_W3mag', u.mag, "Uncertainty on magnitude in W3 filter.")
        self.cols.source.w2 = Column('W2mag', u.mag, "Uncertainty on magnitude in W2 filter.")
        self.cols.source.w2_err = Column('e_W2mag', u.mag, "Uncertainty on magnitude in W2 filter.")
        self.cols.source.w3 = Column('W3mag', u.mag, "Uncertainty on magnitude in W3 filter.")
        self.cols.source.w3_err = Column('e_W3mag', u.mag, "Uncertainty on magnitude in W3 filter.")
        self.cols.source.w4 = Column('W4mag', u.mag, "Uncertainty on magnitude in W4 filter.")
        self.cols.source.w4_err = Column('e_W4mag', u.mag, "Uncertainty on magnitude in W4 filter.")
        self.cols.source.colour_index = Column('W3-W2', None, "W3-W2 colour index")
        self.cols.source.N_2MASS_in_3arcsec = Column('N2', None, "Number of matching 2MASS sources within 3 arcsec")
