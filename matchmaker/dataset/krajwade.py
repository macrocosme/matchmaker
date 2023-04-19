import numpy as np
import astropy.units as u
import pandas as pd

from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe, load, get_matched
from ..measures import luminosity_distance as _luminosity_distance
from ..model.lofar import power_law, sigma as lum_sfr_sigma
from ..model.wise import is_AGN

# Could possibly simply overload the Lotss class as it's almost all the same fields
class Krajwade(Catalog):
    file_folder = DATA_BASE_PATH + 'data/krajwade/'
    file_location = file_folder + 'GalField1-MFS-image.pybdsm.srl.fits'
    name = 'krajwade'
    boxes = None

    def __init__(self, load_data=True, constrain=False, smin=None, single_only=False, input_filename=None):
        # Required and practical columns
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg))
        self.cols.ra_err = Column('E_RA', u.arcsec)
        self.cols.dec_err = Column('E_DEC', u.arcsec)
        self.cols.major = Column('DC_Maj', u.arcsec)
        self.cols.minor = Column('DC_Min', u.arcsec)
        self.cols.pa = Column('DC_PA', u.deg)

        self.beam = {
            'bmaj': 0.001666666666666667 * u.deg,
            'bmin': 0.001666666666666667 * u.deg,
            'bpa': 90 * u.deg
        }

        self.area_dr2 = (4178 * u.deg**2) + (1457 * u.deg**2)
        self.central_frequency = 144 * u.MHz
        self.frequency_band = None

        self.set_survey_specific_columns()

        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.radio.total_flux
        self.cols.measure.total_flux_err = self.cols.radio.e_total_flux

        self.cols.measure.peak_flux = self.cols.radio.peak_flux
        self.cols.measure.peak_flux_err = self.cols.radio.e_peak_flux

        if input_filename is not None:
            self.file_location = self.file_folder + input_filename

        if load_data:
            self.load_data(constrain=constrain, smin=smin, single_only=single_only)

    def load_data(self, constrain=False, smin=None, single_only=False):
        # LOAD LOTSS MAIN TABLE
        self.df = load_fits_as_dataframe(self.file_location)

        if constrain:
            # _smin = 0.8 if smin is None else smin
            if smin is not None:
                const = (self.df[self.cols.measure.total_flux.label] >= smin) # 90% complete point-sources
            if single_only:
                try:
                    const &= (self.df[self.cols.radio.s_code.label] == 'S')  # single component
                except UnboundLocalError:
                    const = (self.df[self.cols.radio.s_code.label] == 'S')  # single component
            df = self.df.loc[const]
            df = df.reset_index().rename(
                {'index': 'index_df_original'},
                axis = 'columns'
            )
            self.df = df

        self.n_source = len(self.df)

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

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.radio = Column()
        self.cols.radio.source_name = Column('Source_Name')
        self.cols.radio.peak_flux = Column('Peak_flux', u.Jy * u.beam**-1)
        self.cols.radio.e_peak_flux = Column('E_Peak_flux', u.Jy * u.beam**-1)
        self.cols.radio.total_flux = Column('Total_flux', u.Jy)
        self.cols.radio.e_total_flux = Column('E_Total_flux', u.Jy)
        self.cols.radio.maj = Column('Maj', u.arcsec)
        self.cols.radio.e_maj = Column('E_Maj', u.arcsec)
        self.cols.radio.min = Column('Min', u.arcsec)
        self.cols.radio.e_min = Column('E_Min', u.arcsec)
        self.cols.radio.pa = Column('PA', u.deg)
        self.cols.radio.e_pa = Column('E_PA', u.deg)
        self.cols.radio.dc_maj = Column('DC_Maj', u.arcsec)
        self.cols.radio.e_dc_maj = Column('E_DC_Maj', u.arcsec)
        self.cols.radio.dc_min = Column('DC_Min', u.arcsec)
        self.cols.radio.e_dc_min = Column('E_DC_Min', u.arcsec)
        self.cols.radio.dc_pa = Column('DC_PA', u.deg)
        self.cols.radio.isl_rms = Column('Isl_rms', u.mJy * u.beam**-1)
        self.cols.radio.s_code = Column('S_Code')
        self.cols.radio.mosaic_id = Column('Mosaic_ID')
        self.cols.radio.masked_fraction = Column('Masked_Fraction')
