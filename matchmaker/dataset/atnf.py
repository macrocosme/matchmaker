import astropy.units as u
from psrqpy import QueryATNF
from . import (Catalog, Column, DATA_BASE_PATH)

class Atnf(Catalog):
    name = 'atnf'

    def __init__(self, load_data=False):
        super().__init__(ra=Column('RAJD', u.deg), dec=Column('DECJD', u.deg))
        self.cols.ra_err = Column('RAJD_ERR', u.deg)
        self.cols.dec_err = Column('DECJD_ERR', u.deg)

        if load_data:
            self.load_data()

    def load_data(self):
        query = QueryATNF(params=['JNAME', 'BNAME', 'RAJD', 'DECJD'])
        self.df = query.dataframe


