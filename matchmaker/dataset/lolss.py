import numpy as np
import astropy.units as u
import pandas as pd

from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe, load, get_matched
from ..measures import luminosity_distance as _luminosity_distance
from ..model.lofar import power_law, sigma as lum_sfr_sigma
from ..model.wise import is_AGN

class Lolss(Catalog):
    file_location = DATA_BASE_PATH + 'data/LoLSS/LoLSS_DR1_rolling_v1.0.srl.fits'
    name = 'lolss'
    boxes = None

    def __init__(self, load_data=True, constrain=False):
        # Required and practical columns
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg))
        self.cols.ra_err = Column('E_RA', u.arcsec)
        self.cols.dec_err = Column('E_DEC', u.arcsec)
        self.cols.major = Column('DC_Maj', u.arcsec)
        self.cols.minor = Column('DC_Min', u.arcsec)
        self.cols.pa = Column('DC_PA', u.deg)

        # self.beam = {
        #     'bmaj': 0.001666666666666667 * u.deg,
        #     'bmin': 0.001666666666666667 * u.deg,
        #     'bpa': 90 * u.deg
        # }

        self.central_frequency = 54 * u.MHz
        # self.frequency_band = None

        self.set_survey_specific_columns()

        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.radio.total_flux
        self.cols.measure.total_flux_err = self.cols.radio.e_total_flux

        self.cols.measure.peak_flux = self.cols.radio.peak_flux
        self.cols.measure.peak_flux_err = self.cols.radio.e_peak_flux

        if load_data:
            self.load_data(constrain=constrain)

    def load_data(self, constrain=False):
        # LOAD LOTSS MAIN TABLE
        self.df = load_fits_as_dataframe(self.file_location)

        if constrain:
            df = self.df.loc[self.df[self.cols.measure.total_flux.label] >= 0.8]  # 90% complete point-sources
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
    
    def set_distance(self, source:Catalog):
        # Add distance columns information
        self.cols.distance = source.cols.distance
        self.cols.z = source.cols.z

        # Add distances as new dataframe column
        lotss_dist = np.empty(len(self.df))
        mask_idx = source.matches[self.name].mask_idx
        filtered_idx = self.matches[source.name].filtered_idx
        lotss_dist[filtered_idx] = source.df.iloc[mask_idx][source.cols.distance.label].values
        self.df[self.cols.distance.label] = lotss_dist
        lotss_dist[filtered_idx] = source.df.iloc[mask_idx][source.cols.z.label].values
        self.df[self.cols.z.label] = lotss_dist
        self.use_distance = True

        sc_source = source.as_SkyCoord(mask_idx)
        sc_target = self.as_SkyCoord(filtered_idx)
        self.matches[source.name].sep3d = sc_source.separation_3d(sc_target) #.to(u.kpc).value
        # self.cols.separation = Column('separation3d-{}'.format(source.name), u.kpc)
        # self.df[self.cols.separation.label] = lotss_dist

    def luminosity_distance(self, mask=None, output_units='W_Hz', with_unit=False):
        df = self.df if mask is None else self.df.iloc[mask]
        return _luminosity_distance(
            df[self.cols.z.label],
            df[self.cols.measure.total_flux.label],
            alpha=-0.8,
            output_units=output_units,
            with_unit=with_unit
        )

    def filter_s_code(self, mask=None, s_code='S'):
        df = self.df if mask is None else self.df.iloc[mask]
        return df['S_Code'].values == s_code

    def set_dr2_regions_boxes(self):
        # lotss region 1 (left)
        l1_ra_max = self.df.loc[self.df[self.cols.ra.label] < 70, self.cols.ra.label].max()
        l1_ra_min = self.df.loc[self.df[self.cols.ra.label] < 70, self.cols.ra.label].min()
        l1_dec_max = self.df.loc[(self.df[self.cols.ra.label] <= l1_ra_max), self.cols.dec.label].max()
        l1_dec_min = self.df.loc[(self.df[self.cols.ra.label] <= l1_ra_max), self.cols.dec.label].min()
        l1_box = [l1_ra_min, l1_ra_max, l1_dec_min, l1_dec_max]

        # lotss region 2 (central)
        mask = (self.df[self.cols.ra.label] > 100) & (self.df[self.cols.ra.label] < 300)
        l2_ra_max = self.df.loc[mask, self.cols.ra.label].max()
        l2_ra_min = self.df.loc[mask, self.cols.ra.label].min()
        mask = (self.df[self.cols.ra.label] > l2_ra_min) & (self.df[self.cols.ra.label] < l2_ra_max)
        l2_dec_max = self.df.loc[mask, self.cols.dec.label].max()
        l2_dec_min = self.df.loc[mask, self.cols.dec.label].min()
        l2_box = [l2_ra_min, l2_ra_max, l2_dec_min, l2_dec_max]

        # lotss region 3 (right)
        l3_ra_max = self.df.loc[self.df[self.cols.ra.label] > 300, self.cols.ra.label].max()
        l3_ra_min = self.df.loc[self.df[self.cols.ra.label] > 300, self.cols.ra.label].min()
        l3_dec_max = self.df.loc[(self.df[self.cols.ra.label] >= l3_ra_min), self.cols.dec.label].max()
        l3_dec_min = self.df.loc[(self.df[self.cols.ra.label] >= l3_ra_min), self.cols.dec.label].min()
        l3_box = [l3_ra_min, l3_ra_max, l3_dec_min, l3_dec_max]

        self.boxes = [l1_box, l2_box, l3_box]

    def compact_sources_idx(self, mask=None):
        sample = self.df if mask is None else self.df.iloc[mask]
        compact_idx = sample.loc[
            # (our_radios['Maj'] <= lotss_beam['bmaj'] * 3600) |
            (sample[self.cols.major.label] <= self.beam['bmaj'].to(self.cols.major.unit))
        ].index

        return compact_idx

    def mask_idx_of_compact_sources(self, catalog_name, mask=None):
        return [np.where(i == self.matches[catalog_name].idx)[0][0] for i in self.compact_sources_idx(mask=mask)]

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.radio = Column()
        self.cols.radio.source_name = Column('Source_Name')
        self.cols.radio.peak_flux = Column('Peak_flux', u.mJy * u.beam**-1)
        self.cols.radio.e_peak_flux = Column('E_Peak_flux', u.mJy * u.beam**-1)
        self.cols.radio.total_flux = Column('Total_flux', u.mJy)
        self.cols.radio.e_total_flux = Column('E_Total_flux', u.mJy)
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
