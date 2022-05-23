import pandas as pd
import astropy.units as u
from frbcat import TNS
from . import (Catalog, Column, DATA_BASE_PATH)

class TnsFrb(Catalog):
    name = 'tns'

    def __init__(self, load_data=False):
        super().__init__(ra=Column('ra_frac', u.deg), dec=Column('dec_frac', u.deg))

        # For some reason, these two columns include units directly.
        self.cols.ra_err = Column('ra_err', None) # add unit
        self.cols.dec_err = Column('decl_err', None) # add unit

        if load_data:
            self.load_data()

    def ra_err(self, with_unit=False, to_unit=None):
        if to_unit is not None:
            ra_err = u.Quantity(self.df['ra_err']).to(to_unit)

        if with_unit:
            return ra_err
        else:
            return ra_err.value

    def dec_err(self, with_unit=False, to_unit=None):
        if to_unit is not None:
            ra_err = u.Quantity(self.df['decl_err']).to(to_unit)

        if with_unit:
            return ra_err
        else:
            return ra_err.value

    def load_data(self):
        tns = TNS(tns_name='macrocosme', tns_id=2641)
        self.df = tns.df
        # units = tns.units
        # self.df = df_chime
