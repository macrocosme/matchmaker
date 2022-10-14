import numpy as np
import astropy.units as u
import pandas as pd

from . import (Catalog, Column, DATA_BASE_PATH)

class Arts(Catalog):
    file_location = DATA_BASE_PATH + 'data/arts/ARTS_one_off_FRB_catalog.csv'
    name = 'arts'
    boxes = None

    def __init__(self, load_data=True):
        # Required and practical columns
        super().__init__(ra=Column('ellipse_ra', u.deg), dec=Column('ellipse_dec', u.deg))
        self.cols.major = Column('ellipse_a', u.arcsec)
        self.cols.minor = Column('ellipse_b', u.arcsec)
        self.cols.pa = Column('ellipse_angle', u.deg)

        self.central_frequency = 1370 * u.MHz
        self.frequency_band = 'L'

        self.set_survey_specific_columns()

        if load_data:
            self.load_data()


    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        a = self.prop_to_unit('major', self.cols.major.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                a = (a * self.cols.major.unit).to(to_unit).value
            else:
                a = a.to(to_unit)

        return a

    def semi_minor(self, mask=None, with_unit=False, to_unit=None):
        b = self.prop_to_unit('minor', self.cols.minor.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                b = (b * self.cols.minor.unit).to(to_unit).value
            else:
                b = b.to(to_unit)

        return b

    def load_data(self):
        # LOAD LOTSS MAIN TABLE
        df = pd.read_csv(self.file_location)
        self.df = df
        self.n_source = len(self.df)

    def set_survey_specific_columns(self):
        # Radio columns
         self.cols.name = Column('NAME', None, None)
         self.cols.tns_name = Column('TNS Name', None, None)
         self.cols.dm = Column('DM [pc cm**-3]', None, None)
         self.cols.sn_amber = Column('S/N (AMBER)', None, None)
         self.cols.ra_j2000 = Column('RA', None, None)
         self.cols.dec_j2000 = Column('DEC', None, None)
         self.cols.pointing_name = Column('POINTING NAME', None, None)
         self.cols.cb_found = Column('CB Found', None, None)
         self.cols.sb_found = Column('SB Found', None, None)
         self.cols.toa = Column('ToA', None, None)
         self.cols.sn_opt = Column('S/N (opt)', None, None)
         self.cols.dm_opt = Column('DM (opt)', None, None)
         self.cols.width = Column('Width (ms)', u.ms, None)
         self.cols.iquv = Column('IQUV', None, None)
         self.cols.rm = Column('RM', None, None)
         self.cols.node_directory = Column('Node Directory', None, None)
         self.cols.central_directory = Column('Central Directory (arts041)', None, None)
         self.cols.detection_mjd = Column('detection_mjd', None, None)
         self.cols.bary_mjd = Column('bary_mjd', None, None)
         self.cols.log_Samber_Smin = Column('log(Samber/Smin)', None, None)
         self.cols.ellipse_ra = Column('ellipse_ra', u.deg, None)
         self.cols.ellipse_dec = Column('ellipse_dec', u.deg, None)
         self.cols.ellipse_a = Column('ellipse_a', u.deg, None)
         self.cols.ellipse_b = Column('ellipse_b', u.deg, None)
         self.cols.ellipse_angle = Column('ellipse_angle', u.deg, None)
         self.cols.ellipse_array = Column('ellipse_array', None, None)
