import numpy as np
import pandas as pd
from ..utils import load_fits_as_dataframe
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..model import wise

class Gaia(Catalog):
    file_location = DATA_BASE_PATH + 'data/krajwade/gaiadr3_galaxies.fits'
    name = 'gaia'

    def __init__(self, load_data=False, constrain=True, use_distance=False):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg), use_distance=use_distance)
        self.set_survey_specific_columns()

        if load_data:
            self.load_data(constrain=constrain)

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self, constrain=False):
        df = load_fits_as_dataframe(self.file_location)
        # df.drop(['recno'], axis=1, inplace=True)
        self.df = df

    def set_survey_specific_columns(self):
        self.cols.source = Column()
        self.cols.source.ra = Column('ra', u.deg, "RA J2000")
        self.cols.source.dec = Column('dec', u.deg, "Dec J2000")

