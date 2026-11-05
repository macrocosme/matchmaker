import astropy.units as u
import numpy as np
import pandas as pd

from ..footprint import Footprint
from ..measures import distance_atomic, powerlaw_scale
from ..measures import luminosity_distance as _luminosity_distance
from ..model.compactness import ENVELOPE_DR2, is_compact
from ..model.lofar import power_law
from ..model.lofar import sigma as lum_sfr_sigma
from ..utils import get_matched, load_fits_as_dataframe
from . import DATA_BASE_PATH, Catalog, Column


class Lotss(Catalog):
    """LoTSS source catalogue, DR2 (default) or DR3.

    release='dr2': LoTSS_DR2_v110_masked.srl.fits, 6" everywhere, footprint as three RA/Dec
    boxes (as in Vohl+2023), area 5,635 deg^2.
    release='dr3': LoTSS_DR3_v1.0.srl.fits (Shimwell+2026). Resolution is 6" or 9" per source
    (Resolution column; not a clean declination cut). The footprint is built from all source
    positions at load time (self.footprint, before any cuts), and self.area comes from it.
    """
    file_location = DATA_BASE_PATH + 'lotss/LoTSS_DR2_v110_masked.srl.fits'
    file_location_dr3 = DATA_BASE_PATH + 'lotss/LoTSS_DR3_v1.0.srl.fits'
    name = 'lotss'
    boxes = None
    footprint = None

    field_file_location = DATA_BASE_PATH + 'lotss/obslist.csv'

    def __init__(self, load_data=True, constrain=False, smin=None, compact_only=False, single_only=False,
                 release='dr2', footprint_nside=256, min_snr=None, compactness_params=None):
        if release not in ('dr2', 'dr3'):
            raise ValueError(f"release must be 'dr2' or 'dr3', not {release!r}")
        self.release = release
        self.footprint_nside = footprint_nside
        # R99.9 envelope coefficients: one tuple, or {resolution_arcsec: tuple}. DR2 uses the
        # published Shimwell+2022 values; DR3 publishes none, so they must be measured first
        # (analysis/07_compactness_refit.py) and passed in.
        self.compactness_params = compactness_params or (ENVELOPE_DR2 if release == 'dr2' else None)
        if release == 'dr3':
            self.file_location = type(self).file_location_dr3
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
        self.area = self.area_dr2 if release == 'dr2' else None   # DR3: set from the footprint on load
        self.central_frequency = 144 * u.MHz
        self.frequency_band = None

        self.set_survey_specific_columns()
        if release == 'dr3':
            # DR3 drops the mosaic bookkeeping and adds per-source resolution and island flux
            for gone in ('mosaic_id', 'masked_fraction'):
                delattr(self.cols.radio, gone)
            self.cols.radio.resolution = Column('Resolution', u.arcsec, 'restoring beam: 6" or 9"')
            self.cols.radio.healpix = Column('HEALPIX', None, 'nside-16 mosaic tile; not a per-source footprint')
            self.cols.radio.isl_total_flux = Column('Isl_Total_flux', u.mJy)
            self.cols.radio.e_isl_total_flux = Column('E_Isl_Total_flux', u.mJy)

        self.cols.measure = Column()
        self.cols.measure.total_flux = self.cols.radio.total_flux
        self.cols.measure.total_flux_err = self.cols.radio.e_total_flux

        self.cols.measure.peak_flux = self.cols.radio.peak_flux
        self.cols.measure.peak_flux_err = self.cols.radio.e_peak_flux

        if load_data:
            self.load_data(constrain=constrain, smin=smin, compact_only=compact_only, single_only=single_only,
                           min_snr=min_snr)

    def load_data(self, constrain=False, smin=None, compact_only=False, single_only=False, min_snr=None):
        """Load the catalogue; with constrain=True apply the requested cuts.

        smin: Total_flux >= smin (mJy). min_snr: Total_flux / Isl_rms >= min_snr, the local
        significance (DR3 detects >95% of compact sources above 9, Shimwell+2026).
        compact_only: DR2 R99.9 envelope (refused for DR3). single_only: S_Code == 'S'.
        """
        # LOAD LOTSS MAIN TABLE
        self.df = load_fits_as_dataframe(self.file_location)

        if self.release == 'dr3':
            # The survey footprint comes from every catalogued source, before any cuts
            self.footprint = Footprint.from_positions(self.df[self.cols.ra.label].values,
                                                      self.df[self.cols.dec.label].values,
                                                      nside=self.footprint_nside,
                                                      values=self.resolution())
            self.area = self.footprint.area

        if constrain:
            if compact_only and self.compactness_params is None:
                raise NotImplementedError(
                    f'no R99.9 envelope for release {self.release}: measure one with '
                    'analysis/07_compactness_refit.py and pass compactness_params')
            const = np.ones(len(self.df), dtype=bool)
            if smin is not None:
                const &= (self.df[self.cols.measure.total_flux.label] >= smin).values  # on Total_flux
            if min_snr is not None:
                snr = self.df[self.cols.measure.total_flux.label] / self.df[self.cols.radio.isl_rms.label]
                const &= (snr >= min_snr).values
            if compact_only:
                self.df['compact'] = self.compact_mask()
                const &= self.df['compact'].values
            if single_only:
                const &= (self.df[self.cols.radio.s_code.label] == 'S').values        # single component
            df = self.df.loc[const]
            df = df.reset_index().rename(
                {'index': 'index_df_original'},
                axis = 'columns'
            )
            self.df = df

        if self.release == 'dr2':
            self.set_dr2_regions_boxes()
        self.n_source = len(self.df)

        # LOAD LOTSS FIELDS LIST
        self.fields = Column()
        # obslist.csv is optional (only get_closest_field uses it); guard as the paper-time version did
        try:
            self.fields.df = pd.read_csv(self.field_file_location)
        except FileNotFoundError:
            self.fields.df = None

    def resolution(self, mask=None):
        """Restoring-beam FWHM in arcsec per row: the Resolution column for DR3, 6 for DR2.

        DR3 stores the 6" beam with conversion noise (6.000000000000001, 5.999999999999997), so values
        are rounded to 0.1" to make them comparable.
        """
        df = self.df if mask is None else self.df.iloc[mask]
        if self.release == 'dr3':
            return np.round(df[self.cols.radio.resolution.label].to_numpy(dtype=float), 1)
        return np.full(len(df), 6.0)

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
            for _l, s in zip(lum_w_hz, sfr):
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

    def compact_mask(self, mask=None):
        """True where ln(S_I/S_P) lies below the R99.9 envelope (Shimwell+2019, 2022).

        Uses self.compactness_params, which may be a single (a, b, c, d) tuple or a dict keyed by
        resolution in arcsec, as DR3 needs (its 6" and 9" imaging have different envelopes).
        """
        if self.compactness_params is None:
            raise NotImplementedError(f'no R99.9 envelope for release {self.release}')
        df = self.df if mask is None else self.df.iloc[mask]
        ln_si_sp = np.log(df[self.cols.measure.total_flux.label] / df[self.cols.measure.peak_flux.label]).values
        snr = (df[self.cols.measure.total_flux.label] / df[self.cols.measure.total_flux_err.label]).values
        if not isinstance(self.compactness_params, dict):
            return is_compact(ln_si_sp, snr, self.compactness_params)
        out = np.zeros(len(df), dtype=bool)
        res = self.resolution(mask=mask)
        for resolution, params in self.compactness_params.items():
            here = res == float(resolution)
            if here.any():
                out[here] = is_compact(ln_si_sp[here], snr[here], params)
        return out

    def filter_mask_to_r99_sigmoid(self, catalog_name):
        """Index (into the matched set) of the matches whose LoTSS source is compact.

        filtered_compact_idx is used together with mask_idx and filtered_idx, e.g.
        self.df.iloc[self.matches[name].filtered_idx[self.matches[name].filtered_compact_idx]]
        """
        compact = self.compact_mask(mask=self.matches[catalog_name].filtered_idx)
        self.matches[catalog_name].filtered_compact_idx = np.where(compact)[0]

    def r99_sigmoid(self, mask=None, return_indices=False):
        """Backwards-compatible wrapper for compact_mask."""
        compact = self.compact_mask(mask=mask)
        return np.where(compact)[0] if return_indices else compact

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
