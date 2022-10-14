import numpy as np
import astropy.units as u
import pandas as pd

from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe, load, get_matched
from ..measures import luminosity_distance as _luminosity_distance
from ..model.lofar import power_law, sigma as lum_sfr_sigma
from ..model.wise import is_AGN

class Califa(Catalog):
    file_location = DATA_BASE_PATH + 'data/califa/califa_dr3.fits'
    name = 'califa'
    boxes = None

    def __init__(self, load_data=True, constrain=False):
        # Required and practical columns
        super().__init__(ra=Column('_RAJ2000', u.deg), dec=Column('_RAJ2000', u.deg))

        self.set_survey_specific_columns()

        if load_data:
            self.load_data(constrain=constrain)

    def load_data(self, constrain=False):
        # LOAD LOTSS MAIN TABLE
        self.df = load_fits_as_dataframe(self.file_location)

        self.n_source = len(self.df)

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols._RAJ2000 = Column('_RAJ2000', u.deg, 'RA J2000 in degrees')
        self.cols._DEJ2000 = Column('_DEJ2000', u.deg, 'DEC J2000 in degrees')
        self.cols.NAME = Column('NAME', None, 'Source name')
        self.cols.ID = Column('ID', int, 'Source ID')
        self.cols.RAJ2000 = Column('RAJ2000', None, 'RA J2000 in hms')
        self.cols.DEJ2000 = Column('DEJ2000', None, 'DEC J2000 in dms')
        self.cols.V500 = Column('V500', None, 'Path to file in V500')
        self.cols.V1200 = Column('V1200', None, 'Path to file in V1200')
        self.cols.COMB = Column('COMB', None, 'Path to file in COMB')
