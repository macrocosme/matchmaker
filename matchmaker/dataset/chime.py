import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)

class Chime(Catalog):
    file_location = DATA_BASE_PATH + 'data/chimefrbcat1.csv'
    name = 'chime'

    def __init__(self, load_data=False):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
        self.cols.ra_err = Column('ra_err', None) # add unit
        self.cols.dec_err = Column('dec_err', None) # add unit

        if load_data:
            self.load_data()

    def load_data(self):
        df_chime = pd.read_csv(self.file_location)
        self.df = df_chime
