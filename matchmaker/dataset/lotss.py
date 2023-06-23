import astropy.units
import numpy as np
import astropy.units as u
import pandas as pd

from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe, load, get_matched
from ..measures import distance_atomic, luminosity_distance as _luminosity_distance, powerlaw_scale
from ..model.lofar import power_law, sigma as lum_sfr_sigma

class Lotss(Catalog):
    """Class describing LoTSS DR2"""

    file_location = DATA_BASE_PATH + 'lotss/LoTSS_DR2_v110_masked.srl.fits'
    name = 'lotss'
    boxes = None

    field_file_location = DATA_BASE_PATH + 'lotss/obslist.csv'

    def __init__(self, load_data=True, constrain=False, smin=None, compact_only=False, single_only=False):
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

        self.hips_dir_high_res = 'https://lofar-surveys.org/public_hips/LoTSS_DR2_high_hips/'
        self.hips_dir_high_res = 'https://lofar-surveys.org/public_hips/LoTSS_DR2_low_hips/'
        self.hips_dir = self.hips_dir_high_res

        self.set_survey_specific_columns()

        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.radio.total_flux
        self.cols.measure.total_flux_err = self.cols.radio.e_total_flux

        self.cols.measure.peak_flux = self.cols.radio.peak_flux
        self.cols.measure.peak_flux_err = self.cols.radio.e_peak_flux

        if load_data:
            self.load_data(constrain=constrain, smin=smin, compact_only=compact_only, single_only=single_only)

    def load_data(self, constrain=False, smin=None, compact_only=False, single_only=False):
        # LOAD LOTSS MAIN TABLE
        df = load_fits_as_dataframe(self.file_location)
        self.df = df

        if constrain:
            # _smin = 0.8 if smin is None else smin
            if smin is not None:
                const = (self.df[self.cols.measure.total_flux.label] >= smin) # 90% complete point-sources

            if compact_only:
                self.df['compact'] = self.r99_sigmoid(return_indices=False)
                try:
                    const &= (self.df['compact'] == True)
                except UnboundLocalError:
                    const = (self.df['compact'] == True)

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

        self.set_dr2_regions_boxes()
        self.n_source = len(self.df)

        # LOAD LOTSS FIELDS LIST
        self.fields = Column()
        try:
            self.fields.df = pd.read_csv(self.field_file_location)
        except FileNotFoundError:
            pass

    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        """Get semi-major axis

        Parameters
        ----------
        mask:numpy.array    (optional)
            Mask describing rows to be returned
        with_unit:bool      (optional)
            Return semi-minor axis with or without a unit (as a quantity)
        to_unit:astropy.units.Unit
            Unit in which the quantity should be returned (can be a conversion)

        Returns
        -------
        semi_major
            Type will depend on input parameters (float, quantity, array, or array of quantities)
        """
        a = self.prop_to_unit('major', self.cols.major.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                a = (a * self.cols.major.unit).to(to_unit).value
            else:
                a = a.to(to_unit)

        return a

    def semi_minor(self, mask:np.array=None, with_unit:bool=False, to_unit:astropy.units.Unit=None):
        """Get semi-minor axis

        Parameters
        ----------
        mask:numpy.array    (optional)
            Mask describing rows to be returned
        with_unit:bool      (optional)
            Return semi-minor axis with or without a unit (as a quantity)
        to_unit:astropy.units.Unit
            Unit in which the quantity should be returned (can be a conversion)

        Returns
        -------
        semi_minor
            Type will depend on input parameters (float, quantity, array, or array of quantities)
        """
        b = self.prop_to_unit('minor', self.cols.minor.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                b = (b * self.cols.minor.unit).to(to_unit).value
            else:
                b = b.to(to_unit)

        return b
    
    def set_distance(self, source:Catalog):
        """Set z and distance columns from another source Catalog

        Parameter
        ---------


        """
        # Add distance columns information

        self.cols.z = source.cols.z
        self.cols.distance = source.cols.distance

        # Add distances as new dataframe column
        lotss_dist = np.empty(len(self.df))
        mask_idx = source.matches[self.name].mask_idx
        filtered_idx = self.matches[source.name].filtered_idx

        distance_atomic(source.df.iloc[mask_idx][source.cols.z.label].values)

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

    def luminosity_distance(self, mask=None, output_units='W_Hz', with_unit=False, with_error=False, log=False):
        df = self.df if mask is None else self.df.iloc[mask]
        l = _luminosity_distance(
                df[self.cols.z.label],
                df[self.cols.measure.total_flux.label],
                alpha=-0.8,
                output_units=output_units,
                with_unit=with_unit
            )
        if not with_error:
            return l
        else:
            l_err_low = _luminosity_distance(
                df[self.cols.z.label],
                df[self.cols.measure.total_flux.label] - df[self.cols.measure.total_flux_err.label],
                alpha=-0.8,
                output_units=output_units,
                with_unit=with_unit
            )
            l_err_high = _luminosity_distance(
                df[self.cols.z.label],
                df[self.cols.measure.total_flux.label] + df[self.cols.measure.total_flux_err.label],
                alpha=-0.8,
                output_units=output_units,
                with_unit=with_unit
            )
            if not log:
                return l, l_err_low, l_err_high
            else:
                return np.log10(l), np.log10(l_err_low), np.log10(l_err_high)

    def radio_loudness(self, obj1, mask_lotss=None, mask1=None, alphas=None, verbose=False):
        matched_obj1, matched_obj2 = get_matched(obj1, self, mask1=mask1, mask2=mask_lotss)
        S = matched_obj2[self.cols.measure.total_flux.label]
        try:
            flux = (S.values * self.cols.measure.total_flux.unit).to(u.Jy).value
            if alphas is not None: #This probably wont work for series -- to be done later. for now, loop.
                scaled_flux = powerlaw_scale(self.central_frequency.to(u.GHz),
                                             1.4 * u.GHz,
                                             flux,
                                             alphas)
                flux = scaled_flux
        except AttributeError:
            flux = (S * self.cols.measure.total_flux.unit).to(u.Jy).value
            scaled_flux = powerlaw_scale(self.central_frequency.to(u.GHz),
                                         1.4 * u.GHz,
                                         flux,
                                         -0.7 if np.isnan(alphas) or alphas is None else alphas)
            flux = scaled_flux

        m = -2.5 * np.log10(flux/(3631 * u.Jy).value) # Should be double checked to scale from 1.4 GHz to 144 MHz

        try:
            g = obj1.df.iloc[mask1][obj1.cols.sdss.modelmag_g.label].values
        except AttributeError:
            g = obj1.df.iloc[mask1][obj1.cols.sdss.modelmag_g.label]
        r_g = 0.4 * (g - m)
        if np.isnan(r_g):
            try:
                g = obj1.df.iloc[mask1][obj1.cols.kron.kronmag_g.label].values
            except AttributeError:
                g = obj1.df.iloc[mask1][obj1.cols.kron.kronmag_g.label]
            r_g = 0.4 * (g - m)
        if verbose:
            print ('s', flux)
            print ('m', m)
            print ('g', g)
            print ('r_g', r_g)
        return r_g

    def distance_to_lum_sfr_relation(self, obj1, mask_lotss=None, mask1=None):
        matched_obj1, matched_obj2 = get_matched(obj1, self, mask1=mask1, mask2=mask_lotss)
        lum_w_hz = self.luminosity_distance(mask=mask_lotss)
        sfr = matched_obj1[obj1.cols.sfr.label]
        if type(mask_lotss) in [int, np.int64, np.int32]:
            return (np.log10(lum_w_hz) - np.log10(power_law(sfr))) / np.log10(lum_sfr_sigma)
        else:
            return (np.log10(lum_w_hz) - np.log10(power_law(sfr).values)) / np.log10(lum_sfr_sigma)
        # for l, pl in zip(lum_w_hz, power_law(sfr)):
        #     print (l, pl, l-pl)
        # return lum_w_hz - power_law(sfr)

    def select_lum_sfr_outliers(self, obj1, mask1=None, mask2=None, n_sigma=1.5, filter_on_sfr=None, store=True, verbose=False):
        matched_obj1, matched_obj2 = get_matched(obj1, self, mask1=mask1, mask2=mask2)
        lum_w_hz = self.luminosity_distance(mask=self.matches[obj1.name].filtered_idx if mask2 is None else mask2)
        sfr = matched_obj1[obj1.cols.sfr.label]
        if verbose:
            nulls = 0
            for l, s in zip(lum_w_hz, sfr):
                if np.isnan(s):
                    nulls += 1
            print ('count(sfr==nan):', nulls)
        # print (lum_w_hz.size, len(sfr))
        cond = np.log10(lum_w_hz) >= np.log10(power_law(sfr)) + (n_sigma*np.log10(lum_sfr_sigma))
        if filter_on_sfr is not None:
            cond &= (sfr > filter_on_sfr)
        if store:
            self.matches[obj1.name].filtered_l_sfr_idx = np.where(cond)[0]
        else:
            return np.where(cond)[0]

    def filter_s_code(self, mask=None, s_code='S'):
        df = self.df if mask is None else self.df.iloc[mask]
        return df['S_Code'].values == s_code

    def si_sp(self, catalog_name=None, mask=None):
        """
        if 'catalog_name' isn't none, it overrides 'mask'
        """
        if catalog_name is not None:
            mask = self.matches[catalog_name].filtered_idx
        matches = self.df.iloc[mask] if mask is not None else self.df
        si_sp = np.log(matches['Total_flux']/matches['Peak_flux'])
        return si_sp

    def sn(self, catalog_name=None, mask=None):
        """
        if 'catalog_name' isn't none, it overrides 'mask'
        """
        if catalog_name is not None:
            mask = self.matches[catalog_name].filtered_idx
        matches = self.df.iloc[mask] if mask is not None else self.df
        s_n = matches['Total_flux']/matches['E_Total_flux']
        return s_n

    def filter_mask_to_r99_sigmoid(self, catalog_name):
        # Subset of source-target cross-section where source size in LoTSS is less or equal to R99.9 sigmoid (Shimwell+2022)
        # matches = self.df.iloc[self.matches[catalog_name].filtered_idx]
        si_sp = self.si_sp(catalog_name)
        s_n = self.sn(catalog_name)
        sigmoid = (0.42 + (1.08 / (1+(s_n/96.57)**2.49)))

        # filtered_compact_idx and be used in conjunction with mask_idx and filtered_idx.
        # E.g. source.df.iloc[self.matches[source.name].mask_idx[self.matches[source.name].filtered_compact_idx]]
        # E.g. self.df.iloc[self.matches[source.name].filtered_idx[self.matches[source.name].filtered_compact_idx]]
        self.matches[catalog_name].filtered_compact_idx = np.where(si_sp <= sigmoid)[0]

    def r99_sigmoid(self, mask=None, return_indices=False):
        df = self.df if mask is None else self.df.iloc[mask]

        si_sp = np.log(df['Total_flux']/df['Peak_flux'])
        s_n = df['Total_flux']/df['E_Total_flux']
        sigmoid = (0.42 + (1.08 / (1+(s_n/96.57)**2.49)))
        if not return_indices:
            return si_sp <= sigmoid
        else:
            return np.where(si_sp <= sigmoid)[0]

    def filter_mask_to_beam_size(self, catalog_name):
        # Subset of source-target cross-section where source size in LoTSS is less or equal to beam size
        matches = self.df.iloc[self.matches[catalog_name].filtered_idx][self.cols.radio.dc_maj.label]
        beam = self.beam['bmaj'].to(self.cols.radio.dc_maj.unit)

        # filtered_compact_idx and be used in conjunction with mask_idx and filtered_idx.
        # E.g. source.df.iloc[self.matches[source.name].mask_idx[self.matches[source.name].filtered_compact_idx]]
        # E.g. self.df.iloc[self.matches[source.name].filtered_idx[self.matches[source.name].filtered_compact_idx]]
        self.matches[catalog_name].filtered_compact_idx = np.where(matches <= beam)[0]

    def get_closest_field(self, ra, dec, min_date=None):
        # Adapted from H. Vedantham's script
        cd = np.cos(0.5 * np.pi/180.0 * (dec+self.fields.df[self.cols.dec.label]))
        dis = (cd**2 * (self.fields.df[self.cols.ra.label]-ra)**2 + (self.fields.df[self.cols.dec.label]-dec)**2)**0.5
        I = np.where(dis < 2)[0]

        if len(I) > 0:
            idx = None
            if min_date is not None:
                # Check if date is before observed time
                for ii in I[np.argsort(dis[I].values)]:
                    if min_date < self.fields.df.iloc[ii]['date_end']:
                        idx = ii
                if idx is not None:
                    # Match found within time frame
                    return self.fields.df.iloc[idx]['field'], self.fields.df.iloc[idx]['date_end'] #, I, dis[I]
                else:
                    # No observation after min_date
                    return None, None
            else:
                # Otherwise return nearest field center
                return self.fields.df.iloc[I[np.argmin(dis[I].values)]]['field'], None
        else:
            # No field found
            return None, None

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

        self.cols.optical_id = Column()
        self.cols.optical_id.id_flag = Column('ID_flag')
        self.cols.optical_id.id_name = Column('ID_name')
        self.cols.optical_id.id_ra = Column('ID_ra', u.deg)
        self.cols.optical_id.id_dec = Column('ID_dec', u.deg)
        self.cols.optical_id.ml_lr = Column('ML_LR')
        self.cols.optical_id.lgz_size = Column('LGZ_Size', u.arcsec)
        self.cols.optical_id.lgz_width = Column('LGZ_Width', u.arcsec)
        self.cols.optical_id.lgz_pa = Column('LGZ_PA', u.deg)
        self.cols.optical_id.lgz_assoc = Column('LGZ_Assoc')
        self.cols.optical_id.lgz_assoc_qual = Column('LGZ_Assoc_Qual')
        self.cols.optical_id.lgz_id_qual = Column('LGZ_ID_Qual')

        self.cols.optical_photometry = Column()
        self.cols.optical_photometry.allwise = Column('AllWISE')
        self.cols.optical_photometry.objid = Column('objid')
        self.cols.optical_photometry.gfapflux = Column('gFApFlux', u.uJy)
        self.cols.optical_photometry.gfapfluxerr = Column('gFApFluxErr', u.uJy)
        self.cols.optical_photometry.gfapmag = Column('gFApMag', u.mag, 'AB mag')
        self.cols.optical_photometry.gfapmagerr = Column('gFApMagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.rfapflux = Column('rFApFlux', u.uJy)
        self.cols.optical_photometry.rfapfluxerr = Column('rFApFluxErr', u.uJy)
        self.cols.optical_photometry.rfapmag = Column('rFApMag', u.mag, 'AB mag')
        self.cols.optical_photometry.rfapmagerr = Column('rFApMagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.ifapflux = Column('iFApFlux', u.uJy)
        self.cols.optical_photometry.ifapfluxerr = Column('iFApFluxErr', u.uJy)
        self.cols.optical_photometry.ifapmag = Column('iFApMag', u.mag, 'AB mag')
        self.cols.optical_photometry.ifapmagerr = Column('iFApMagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.zfapflux = Column('zFApFlux', u.uJy)
        self.cols.optical_photometry.zfapfluxerr = Column('zFApFluxErr', u.uJy)
        self.cols.optical_photometry.zfapmag = Column('zFApMag', u.mag, 'AB mag')
        self.cols.optical_photometry.zfapmagerr = Column('zFApMagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.yfapflux = Column('yFApFlux', u.uJy)
        self.cols.optical_photometry.yfapfluxerr = Column('yFApFluxErr', u.uJy)
        self.cols.optical_photometry.yfapmag = Column('yFApMag', u.mag, 'AB mag')
        self.cols.optical_photometry.yfapmagerr = Column('yFApMagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.w1flux = Column('w1Flux', u.uJy)
        self.cols.optical_photometry.w1fluxerr = Column('w1FluxErr', u.uJy)
        self.cols.optical_photometry.w1mag = Column('w1Mag', u.mag, 'AB mag')
        self.cols.optical_photometry.w1magerr = Column('w1MagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.w2flux = Column('w2Flux', u.uJy)
        self.cols.optical_photometry.w2fluxerr = Column('w2FluxErr', u.uJy)
        self.cols.optical_photometry.w2mag = Column('w2Mag', u.mag, 'AB mag')
        self.cols.optical_photometry.w2magerr = Column('w2MagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.w3flux = Column('w3Flux', u.uJy)
        self.cols.optical_photometry.w3fluxerr = Column('w3FluxErr', u.uJy)
        self.cols.optical_photometry.w3mag = Column('w3Mag', u.mag, 'AB mag')
        self.cols.optical_photometry.w3magerr = Column('w3MagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.w4flux = Column('w4Flux', u.uJy)
        self.cols.optical_photometry.w4fluxerr = Column('w4FluxErr', u.uJy)
        self.cols.optical_photometry.w4magerr = Column('w4MagErr', u.mag, 'AB mag')
        self.cols.optical_photometry.gfkronflux = Column('gFKronFlux', u.uJy)
        self.cols.optical_photometry.gfkronfluxerr = Column('gFKronFluxErr', u.uJy)
        self.cols.optical_photometry.rfkronflux = Column('rFKronFlux', u.uJy)
        self.cols.optical_photometry.rfkronfluxerr = Column('rFKronFluxErr', u.uJy)
        self.cols.optical_photometry.ifkronflux = Column('iFKronFlux', u.uJy)
        self.cols.optical_photometry.ifkronfluxerr = Column('iFKronFluxErr', u.uJy)
        self.cols.optical_photometry.zfkronflux = Column('zFKronFlux', u.uJy)
        self.cols.optical_photometry.zfkronfluxerr = Column('zFKronFluxErr', u.uJy)
        self.cols.optical_photometry.yfkronflux = Column('yFKronFlux', u.uJy)
        self.cols.optical_photometry.yfkronfluxerr = Column('yFKronFluxErr', u.uJy)

        self.cols.photo_z = Column()
        self.cols.photo_z.z_best = Column('z_best')
        self.cols.photo_z.z_best_source = Column('z_best_source')
        self.cols.photo_z.z_spec = Column('z_spec')
        self.cols.photo_z.z_spec_source = Column('z_spec_source')
        self.cols.photo_z.z1_median = Column('z1_median')
        self.cols.photo_z.z1_min = Column('z1_min')
        self.cols.photo_z.z1_max = Column('z1_max')
        self.cols.photo_z.z1_area = Column('z1_area')
        self.cols.photo_z.z2_median = Column('z2_median')
        self.cols.photo_z.z2_min = Column('z2_min')
        self.cols.photo_z.z2_max = Column('z2_max')
        self.cols.photo_z.z2_area = Column('z2_area')
        self.cols.photo_z.specagn = Column('specAGN')
        self.cols.photo_z.mqcagn = Column('mqcAGN')
        self.cols.photo_z.xrayclass = Column('XrayClass')
        self.cols.photo_z._2rxs_id = Column('2RXS_ID')
        self.cols.photo_z.xmmsl2_id = Column('XMMSL2_ID')
        self.cols.photo_z.irclass = Column('IRClass')
        self.cols.photo_z.ebv = Column('EBV')
        self.cols.photo_z.panstarrs_missing = Column('PanSTARRS_Missing')

        self.cols.restframe_magnitudes = Column()
        self.cols.restframe_magnitudes.u_rest = Column('u_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.g_rest = Column('g_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.r_rest = Column('r_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.i_rest = Column('i_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.z_rest = Column('z_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.w1_rest = Column('w1_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.w2_rest = Column('w2_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.w3_rest = Column('w3_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.u_rest = Column('U_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.b_rest = Column('B_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.v_rest = Column('V_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.i_rest = Column('I_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.j_rest = Column('J_rest', u.mag, 'AB mag')
        self.cols.restframe_magnitudes.ks_rest = Column('Ks_rest', u.mag, 'AB mag')
