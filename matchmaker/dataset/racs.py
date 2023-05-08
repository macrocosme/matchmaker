import numpy as np
import astropy.units as u
import pandas as pd

from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe, load, get_matched
from ..measures import luminosity_distance as _luminosity_distance
from ..model.lofar import power_law, sigma as lum_sfr_sigma
from ..model.wise import is_AGN

class Racs(Catalog):
    file_location = DATA_BASE_PATH + 'racs/racs.fits'
    name = 'racs'
    boxes = None

    def __init__(self, load_data=True, constrain=False):
        # Required and practical columns
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))

        self.set_survey_specific_columns()
        self.central_frequency = (1800+700)//2 * u.MHz      # From McConnel et al. (2020) PASA
        self.bandwidth = (1800-700) * u.MHz

        if load_data:
            self.load_data(constrain=constrain)

    def load_data(self, constrain=False):
        # LOAD LOTSS MAIN TABLE
        self.df = load_fits_as_dataframe(self.file_location)

        self.n_source = len(self.df)

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.ra = Column('ra', u.deg, 'RA J2000 in degrees')
        self.cols.dec = Column('dec', u.deg, 'DEC J2000 in degrees')
        self.cols.ra_err = Column('e_ra', u.deg, 'RA J2000 in degrees')
        self.cols.dec_err = Column('e_dec', u.deg, 'DEC J2000 in degrees')
        self.cols.measure = Column()
        self.cols.measure.total_flux = Column('total_flux', u.mJy)
        self.cols.measure.total_flux_err = Column('e_total_flux', u.mJy)
        self.cols.measure.peak_flux = Column('peak_flux', u.mJy * u.beam**-1)
        self.cols.measure.peak_flux_err = Column('e_peak_flux', u.mJy * u.beam**-1)
