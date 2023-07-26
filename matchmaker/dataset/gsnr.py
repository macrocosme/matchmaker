import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_fits_as_dataframe

class Gsnr(Catalog):
    file_location = DATA_BASE_PATH + 'galactic_snr/galactic_snr.fits'
    name = 'galactic_snr'

    def __init__(self, load_data=False, constrain=True, use_a=False, use_distance=True, constrain_mass=True, constrain_sfr=True):
        super().__init__(ra=Column('_RA_icrs', u.deg), dec=Column('_DE_icrs', u.deg), use_distance=use_distance)

        self.set_survey_specific_columns()

        self.cols.ra = self.cols._RA_icrs
        self.cols.dec = self.cols._DEC_icrs

        if load_data:
            self.load_data()

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self):
        df = load_fits_as_dataframe(self.file_location)
        self.df = df

    def set_survey_specific_columns(self):
        self.cols.recno = Column('recno', None, '')
        self.cols.Set = Column('Set', None, 'Set')
        self.cols.Seq = Column('Seq', None, 'Running sequence number')
        self.cols.SNR = Column('SNR', None, 'Position-based Supernova remnant name (6)')
        self.cols.Size1 = Column('Size1', u.arcmin, 'Major axis')
        self.cols.Size2 = Column('Size2', u.arcmin, 'Minor axis')
        self.cols.Type = Column('Type', None, 'Supernova remnant Type (2)')
        self.cols.l_S1GHz = Column('l_S1GHz', None, 'Limit flag on S1GHz')
        self.cols.S1GHz = Column('S1GHz', u.Jy, '1GHz Flux density')
        self.cols.q_S1GHz = Column('q_S1GHz', None, 'Quality flag on S1GHz')
        self.cols.alpha = Column('alpha', None, 'Spectral index, {alpha}')
        self.cols.u_alpha = Column('u_alpha', None, 'Quality flag on alpha')
        self.cols.e_alpha = Column('e_alpha', None, 'Lower Uncertainty on alpha')
        self.cols.E_alpha = Column('E_alpha', None, 'Upper Uncertainty on alpha')
        self.cols.r_alpha = Column('r_alpha', None, 'Reference for alpha (bibcode)')
        self.cols.l_X = Column('l_X', None, '')
        self.cols.X = Column('X', u.kpc, 'Heliocentric distance')
        self.cols.f_X = Column('f_X', None, '')
        self.cols.e_X = Column('e_X', u.kpc, 'Lower Uncertainty on X')
        self.cols.E_X = Column('E_X', u.kpc, 'Upper Uncertainty on X')
        self.cols.l_Y = Column('l_Y', None, '')
        self.cols.Y = Column('Y', u.kpc, 'Galactocentric distance')
        self.cols.f_Y = Column('f_Y', None, '')
        self.cols.e_Y = Column('e_Y', u.kpc, 'Lower Uncertainty on Y')
        self.cols.E_Y = Column('E_Y', u.kpc, 'Upper Uncertainty on Y')
        self.cols.l_Z = Column('l_Z', None, '')
        self.cols.Z = Column('Z', u.pc, 'Height above Galaxy mid-plane')
        self.cols.e_Z = Column('e_Z', u.pc, 'Lower Uncertainty on Z')
        self.cols.E_Z = Column('E_Z', u.pc, 'Upper Uncertainty on Z')
        self.cols.MCInt = Column('MCInt', None, '')
        self.cols.r_MCInt = Column('r_MCInt', None, '')
        self.cols.l_Age = Column('l_Age', None, '')
        self.cols.Age = Column('Age', u.kyr, 'Age')
        self.cols.f_Age = Column('f_Age', None, '')
        self.cols.e_Age = Column('e_Age', u.kyr, 'Lower Uncertainty on Age')
        self.cols.E_Age = Column('E_Age', u.kyr, 'Upper Uncertainty on Age')
        self.cols.n_Age = Column('n_Age', None, '')
        self.cols.r_Age = Column('r_Age', None, '')
        self.cols.SN = Column('SN', None, '')
        self.cols.r_SN = Column('r_SN', None, '')
        self.cols.l_Rad = Column('l_Rad', None, '')
        self.cols.Rad = Column('Rad', u.pc, 'Average radius (4)')
        self.cols.e_Rad = Column('e_Rad', u.pc, 'Lower Uncertainty on Rad')
        self.cols.E_Rad = Column('E_Rad', u.pc, 'Upper Uncertainty on Rad')
        self.cols.SB = Column('SB', 1e-21*u.W * u.m**-2 * u.Hz**-1 * u.sr**-1, 'Surface brightness')
        self.cols.u_SB = Column('u_SB', None, '')
        self.cols.Oval = Column('Oval', 100**-1, 'Ovality')
        self.cols.SNRcat = Column('SNRcat', None, '')
        self.cols.SimbadName = Column('SimbadName', None, '')
        self.cols._Glon = Column('_Glon', u.deg, 'Positions from SNR names (longitude part)')
        self.cols._Glat = Column('_Glat', u.deg, 'Positions from SNR names (latitude part)')
        self.cols._RA_icrs = Column('_RA_icrs', u.deg, 'Right ascension (ICRS) (computed by VizieR, not part of the original data. The format may include more digits than the original data because of internal accuracy requirements in VizieR and across other CDS services)')
        self.cols._DE_icrs = Column('_DE_icrs', u.deg, 'Declination (ICRS) (computed by VizieR, not part of the original data. The format may include more digits than the original data because of internal accuracy requirements in VizieR and across other CDS services)')

