import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from pandas import read_csv


class Asiago(Catalog):
    file_location = DATA_BASE_PATH + 'data/asiago-supernova/asiago-sn.csv'
    name = 'asiago'

    def __init__(self, load_data=False, constrain=True):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))
        self.cols.sn_type = Column('sn_type')
        self.cols.unconfirmed_flag = Column('unconfirmed_flag')
        self.cols.galaxy_name = Column('galaxy_name')
        self.cols.max_epoch = Column('max_epoch')
        self.cols._class = Column('class')

        if load_data:
            self.load_data()

    def load_data(self):
        df = read_csv(self.file_location)

        self.df = df
