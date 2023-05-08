import numpy as np
import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..model import wise

class Clu(Catalog):
    file_location = DATA_BASE_PATH + 'clu/CLU_20190708_marshalFormat.csv'
    name = 'clu'

    def __init__(self, load_data=False, constrain=True, use_a=False, use_distance=True, constrain_mass=True, constrain_sfr=True):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg), use_distance=use_distance)
        self.cols.a = Column('a', u.arcmin)
        self.cols.b = Column('b', u.arcmin)
        self.cols.b2a = Column('b2a')
        self.cols.pa = Column('pa', u.deg)
        self.cols.z = Column('z')
        self.cols.mstar = Column('mstar', u.solMass)
        self.cols.distance = Column('distmpc', u.Mpc)

        self.set_survey_specific_columns()
        self.cols.sfr = self.cols.measure.sfr_fuv

        self.central_frequency = 670 * u.nm
        self.area_dr2 = (4178 * u.deg**2) + (1457 * u.deg**2)

        if load_data:
            self.load_data(constrain=constrain, use_a=use_a, constrain_mass=constrain_mass, constrain_sfr=constrain_sfr)

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self, constrain=False, constrain_mass=True, constrain_sfr=True, use_a=True):
        df = pd.read_csv(self.file_location)
        self.df = df
        if constrain:
            if use_a:
                df = self.df.loc[
                    # (self.df['mstar'].apply(np.log10) > 7) &
                    # (self.df['cluhamag'] < 14) &
                    (self.df['mstar'] <= 3e9) &
                    (~self.df['a'].isna())
                ]
            else:
                if constrain_mass:
                    const = (self.df['mstar'].apply(np.log10) > 7) & (self.df['mstar'] <= 3e9)
                if constrain_sfr:
                    try:
                        const &= (~self.df['sfr_fuv'].isna())
                    except:
                        const = (~self.df['sfr_fuv'].isna())
                df = self.df.loc[const]
            self.df = df

    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        is_array = (type(mask) is np.ndarray)

        a = self.df[self.cols.a.label].values if mask is None else \
            self.df.iloc[mask][self.cols.a.label].values if is_array else self.df.iloc[mask][self.cols.a.label]
        a = a * self.cols.a.unit
        if to_unit is not None:
            a = a.to(to_unit)

        return a if with_unit else a.value

    def semi_minor(self, mask=None, with_unit=False, to_unit=None):
        if type(mask) is not np.ndarray:
            is_array = False

        a = self.df[self.cols.a.label].values if mask is None else \
            self.df.iloc[mask][self.cols.a.label].values if is_array else self.df.iloc[mask][self.cols.a.label]
        b2a = self.df[self.cols.b2a.label].values if mask is None else self.df.iloc[mask][self.cols.b2a.label]
        if is_array:
            b2a.fillna(1, inplace=True)
            b2a = b2a.values
        minor = (a * b2a) * self.cols.a.unit

        if to_unit is not None:
            minor = minor.to(to_unit)

        return minor if with_unit else minor.value
    #
    # def ww(self, mask=None):
    #     df = self.df if mask is None else self.df.iloc[mask]
    #     w1w2 = df['w1mpro'] - df['w2mpro']
    #     w2w3 = df['w2mpro'] - df['w3mpro']
    #
    #     # This needs to change to actual AGN box (wise.AGN_color_color_box)
    #     return False if len(np.where(w1w2 < 0.75)[0]) > 0 else True

    def non_agn_mask(self, mask=None):
        from matplotlib.path import Path
        df = self.df if mask is None else self.df.iloc[mask]
        w1w2 = df['w1mpro'] - df['w2mpro']
        w2w3 = df['w2mpro'] - df['w3mpro']

        if type(w2w3) in [float, np.float32, np.float64]:
            return ~Path(wise.AGN_color_color_box).contains_point([w2w3, w1w2])
        else:
            return ~Path(wise.AGN_color_color_box).contains_points(np.array([w2w3, w1w2]).T)


    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.source = Column()
        self.cols.source.cluid = Column('cluid', None, "unique id number")
        self.cols.source.id_other = Column('id_other', None, "pointing to old ID numbers in original version")
        self.cols.source.name = Column('name', None, "only use names from actual surveys; use name where we got RA/Dec from")

        self.cols.coords = Column()
        self.cols.coords.ra = Column('ra', u.deg, "right ascension")
        self.cols.coords.dec = Column('dec', u.deg, "declination")

        self.cols.ellipse = Column()
        self.cols.ellipse.a = Column('a', u.arcmin, "semi-major axis (This is a D25 size)")
        self.cols.ellipse.b2a = Column('b2a', None, "axis ratio")
        self.cols.ellipse.pa = Column('pa', None, "position angle")

        self.cols.classification = Column()
        self.cols.classification.type_ned = Column('type_ned', None, "from NED (G, Irr, SB, SBa, etc)")

        self.cols.galex = Column()
        self.cols.galex.name_galex = Column('name_galex', None, "nan")
        self.cols.galex.ra_galex = Column('ra_galex', None, "nan")
        self.cols.galex.dec_galex = Column('dec_galex', None, "nan")
        self.cols.galex.fuv = Column('fuv', None, "GALEX kron mag in Far UV (FUV) filter")
        self.cols.galex.fuverr = Column('fuverr', None, "nan")
        self.cols.galex.nuv = Column('nuv', None, "GALEX kron mag in Near UV (NUV) filter")
        self.cols.galex.nuverr = Column('nuverr', None, "nan")

        self.cols.sdss = Column()
        self.cols.sdss.name_sdss = Column('name_sdss', None, "nan")
        self.cols.sdss.ra_sdss = Column('ra_sdss', None, "nan")
        self.cols.sdss.dec_sdss = Column('dec_sdss', None, "nan")
        self.cols.sdss.modelmag_u = Column('modelmag_u', None, "SDSS model mag in u-band")
        self.cols.sdss.modelmagerr_u = Column('modelmagerr_u', None, "nan")
        self.cols.sdss.modelmag_g = Column('modelmag_g', None, "SDSS model mag in g-band")
        self.cols.sdss.modelmagerr_g = Column('modelmagerr_g', None, "nan")
        self.cols.sdss.modelmag_r = Column('modelmag_r', None, "SDSS model mag in r-band")
        self.cols.sdss.modelmagerr_r = Column('modelmagerr_r', None, "nan")
        self.cols.sdss.modelmag_i = Column('modelmag_i', None, "SDSS model mag in i-band")
        self.cols.sdss.modelmagerr_i = Column('modelmagerr_i', None, "nan")
        self.cols.sdss.modelmag_z = Column('modelmag_z', None, "SDSS model mag in z-band")
        self.cols.sdss.modelmagerr_z = Column('modelmagerr_z', None, "nan")

        self.cols.ps1 = Column()
        self.cols.ps1.name_ps1 = Column('name_ps1', None, "nan")
        self.cols.ps1.ra_ps1 = Column('ra_ps1', None, "nan")
        self.cols.ps1.dec_ps1 = Column('dec_ps1', None, "nan")

        self.cols.kron = Column()
        self.cols.kron.kronmag_g = Column('kronmag_g', None, "ps1 mag in g-band")
        self.cols.kron.kronmagerr_g = Column('kronmagerr_g', None, "nan")
        self.cols.kron.kronmag_r = Column('kronmag_r', None, "ps1 mag in r-band")
        self.cols.kron.kronmagerr_r = Column('kronmagerr_r', None, "nan")
        self.cols.kron.kronmag_i = Column('kronmag_i', None, "ps1 mag in i-band")
        self.cols.kron.kronmagerr_i = Column('kronmagerr_i', None, "nan")
        self.cols.kron.kronmag_z = Column('kronmag_z', None, "ps1 mag in z-band")
        self.cols.kron.kronmagerr_z = Column('kronmagerr_z', None, "nan")
        self.cols.kron.kronmag_y = Column('kronmag_y', None, "ps1 mag in y-band")
        self.cols.kron.kronmagerr_y = Column('kronmagerr_y', None, "nan")

        self.cols.b_band = Column()
        self.cols.b_band.btc = Column('btc', None, "the total B-band mag [generally from RC3 catalog as tabulated by NED.]")
        self.cols.b_band.btcerr = Column('btcerr', None, "nan")
        self.cols.b_band.b_r25 = Column('b_r25', None, "The magnitude within the isophotal ellipse at 25 mag/arcsec^2 [from the RC3 catalog as tabulated by NED.]")
        self.cols.b_band.b_r25err = Column('b_r25err', None, "nan")
        self.cols.b_band.magb = Column('magb', None, "AVS need to know how Dave converted SDSS+btc into this")
        self.cols.b_band.magberr = Column('magberr', None, "nan")

        self.cols.twomass = Column()
        self.cols.twomass.name_2mass = Column('name_2mass', None, "nan")
        self.cols.twomass.ra_2mass = Column('ra_2mass', None, "nan")
        self.cols.twomass.dec_2mass = Column('dec_2mass', None, "nan")

        self.cols.k20 = Column()
        self.cols.k20.r_k20fe = Column('r_k20fe', None, "nan")
        self.cols.k20.j_m_k20fe = Column('j_m_k20fe', None, "2mass J-band 20 mag/sq” isophotal mag")
        self.cols.k20.j_msig_k20fe = Column('j_msig_k20fe', None, "mag error")
        self.cols.k20.j_flg_k20fe = Column('j_flg_k20fe', None, "nan")
        self.cols.k20.h_m_k20fe = Column('h_m_k20fe', None, "2mass H-band 20 mag/sq” isophotal mag")
        self.cols.k20.h_msig_k20fe = Column('h_msig_k20fe', None, "mag error")
        self.cols.k20.h_flg_k20fe = Column('h_flg_k20fe', None, "nan")
        self.cols.k20.k_m_k20fe = Column('k_m_k20fe', None, "2mass K-band 20 mag/sq” isophotal mag")
        self.cols.k20.k_msig_k20fe = Column('k_msig_k20fe', None, "mag error")
        self.cols.k20.k_flg_k20fe = Column('k_flg_k20fe', None, "nan")

        self.cols.wise = Column()
        self.cols.wise.name_wise = Column('name_wise', None, "nan")
        self.cols.wise.ra_wise = Column('ra_wise', None, "nan")
        self.cols.wise.dec_wise = Column('dec_wise', None, "nan")
        self.cols.wise.w1mpro = Column('w1mpro', None, "profile fitted magnitude in WISE1 band")
        self.cols.wise.w1sigmpro = Column('w1sigmpro', None, "mag error")
        self.cols.wise.w1snr = Column('w1snr', None, "nan")
        self.cols.wise.w2mpro = Column('w2mpro', None, "profile fitted magnitude in WISE2 band")
        self.cols.wise.w2sigmpro = Column('w2sigmpro', None, "mag error")
        self.cols.wise.w2snr = Column('w2snr', None, "nan")
        self.cols.wise.w3mpro = Column('w3mpro', None, "profile fitted magnitude in WISE3 band")
        self.cols.wise.w3sigmpro = Column('w3sigmpro', None, "mag error")
        self.cols.wise.w3snr = Column('w3snr', None, "nan")
        self.cols.wise.w4mpro = Column('w4mpro', None, "profile fitted magnitude in WISE4 band")
        self.cols.wise.w4sigmpro = Column('w4sigmpro', None, "mag error")
        self.cols.wise.w4snr = Column('w4snr', None, "nan")

        self.cols.HI = Column()
        self.cols.HI.lum_b = Column('lum_b', None, "Luminosity in the B-band from the MAGB column given the distance in DM_KIN column")
        self.cols.HI.m21 = Column('m21', None, "HI mass (21cm)")
        self.cols.HI.m21err = Column('m21err', None, "HI mass (21cm) error")

        self.cols.measure = Column()
        self.cols.measure.dm_kin = Column('dm_kin', None, "distance modulus based on redshift")
        self.cols.measure.dm = Column('dm', None, "the \"best\" distance modulus (can be either redshift-based or redshift-independent)")
        self.cols.measure.dm_method = Column('dm_method', None, "associated with dm column method (e.g., redshift, Tully-Fisher, TRGB, cepheids)")
        self.cols.measure.z = Column('z', None, "Redshift")
        self.cols.measure.zerr = Column('zerr', None, "Redshift uncertainty")
        self.cols.measure.distmpc = Column('distmpc', u.Mpc, "distance in Mpc (This will be populated using the 'dm' column)")
        self.cols.measure.sfr_fuv = Column('sfr_fuv', u.Msun * u.yr**-1, "star formation rate from GALEX FUV fluxes corrected for internal dust extinction via WISE, 22um (u.Msun * u.yr**-1)")
        self.cols.measure.sfr_fuverr = Column('sfr_fuverr', u.Msun * u.yr**-1, "star formation rate uncertainty from GALEX FUV fluxes corrected for internal dust extinction via WISE, 22um (u.Msun * u.yr**-1)")
        self.cols.measure.sfr_ha = Column('sfr_ha', u.Msun * u.yr**-1, "star formation rate from CLU-Ha fluxes corrected for internal dust extinction via WISE 22um, (u.Msun * u.yr**-1)")
        self.cols.measure.mstar = Column('mstar', u.Msun, "Stellar Mass (Msun)")
        self.cols.measure.mstar_err = Column('mstarerr', u.Msun, "Stellar Mass uncertainty (Msun)")
        self.cols.measure.size_source = Column('size_source', None, "where was the radii taken from (a, b2a, pa)")
        self.cols.measure.btc_source = Column('btc_source', None, "nan")
        self.cols.measure.dm_source = Column('dm_source', None, "distance modulus source (main source to use)")
        self.cols.measure.z_source = Column('z_source', None, "redshift source")

        self.cols.source.source = Column('source', None, "copy of coord_source for backwards compatibility")
        self.cols.source.name_clu = Column('name_clu', None, "name from CLU-Ha")
        self.cols.source.name_clu_err = Column('name_clu_err', None, "nan")
        self.cols.source.ra_clu = Column('ra_clu', None, "nan")
        self.cols.source.dec_clu = Column('dec_clu', None, "nan")

        self.cols.color = Column()
        self.cols.color.maxcsig = Column('maxcsig', None, "the color sigma in the ‘On’ filters, which indicates how significant the presence of an Ha line is.  Can do a crossmatch for all galaxies in CLU compiled")

        self.cols.halpha = Column()
        self.cols.halpha.cluhamag = Column('cluhamag', None, "Mag in the ‘On’ filter")
        self.cols.halpha.cluhamagerr = Column('cluhamagerr', None, "nan")

        self.cols.other = Column()
        self.cols.other.flag = Column('flag', None, "(binary; description yet to come (e.g., 0=good, 1=bad size, etc))")
