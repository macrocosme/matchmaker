import numpy as np
import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..model import wise

class Glade(Catalog):
    file_location = DATA_BASE_PATH + 'data/glade/glade+_trimmed.csv'
    name = 'glade'

    def __init__(self, load_data=False, constrain=True, use_distance=True):
        super().__init__(ra=Column('RA', u.deg), dec=Column('RA', u.deg), use_distance=use_distance)

        self.set_survey_specific_columns()

        self.cols.z = self.cols.all.z_helio
        self.cols.mstar = self.cols.all.mstar
        self.cols.distance = self.cols.all.d_L

        # self.central_frequency = 670 * u.nm
        self.area_dr2 = (4178 * u.deg**2) + (1457 * u.deg**2)

        if load_data:
            self.load_data(constrain=constrain)

    @property
    def n_source(self):
        return len(self.df)

    def load_data(self, constrain=False):
        cols = ['GLADE_no', 'PGC_no', 'GWGC_name', 'HyperLEDA_name', 'twoMASS_name', 'WISExSCOS_name',
                'SDSS_DR16Q_name', 'Object_type_flag', 'RA', 'Dec', 'B', 'B_err', 'B', 'B_Abs', 'J', 'J_err', 'H',
                'H_err', 'K', 'K_err', 'W1', 'W1_err', 'W2', 'W2_err', 'W1_flag', 'B_J', 'B_J_err', 'z_helio', 'z_cmb',
                'z_flag', 'v_err', 'z_err', 'd_L', 'd_L_err', 'dist_flag', 'mstar', 'mstar_err', 'merger_rate',
                'merger_rate_error']
        df = pd.read_csv(self.file_location, )
        self.df = df
        if constrain:
            df = self.df.loc[
                ((self.df['mstar'] * self.cols.all.mstar.unit).values < 1e9)
            ]
            self.df = df

    def set_survey_specific_columns(self):
        self.cols.all = Column()
        
        self.cols.all.GLADE_no = Column('GLADE_no', None, 'GLADE+ catalog number')
        self.cols.all.PGC_no = Column('PGC_no', None, 'Principal Galaxies Catalogue number')
        self.cols.all.GWGC_name = Column('GWGC_name', None, 'Name in the GWGC catalog')
        self.cols.all.HyperLEDA_name = Column('HyperLEDA_name', None, 'Name in the HyperLEDA catalog')
        self.cols.all.twoMASS_name = Column('2MASS_name', None, 'Name in the 2MASS XSC catalog')
        self.cols.all.WISExSCOS_name = Column('WISExSCOS_name', None, 'Name in the WISExSuperCOSMOS catalog (wiseX)')
        self.cols.all.SDSS_DR16Q_name = Column('SDSS_DR16Q_name', None, 'Name in the SDSS-DR16Q catalog')
        self.cols.all.Object_type_flag = Column('Object_type_flag', None, 'Q: the source is from the SDSS-DR16Q catalog; G:the source is from another catalog and has not been identified as a quasar')
        self.cols.all.RA = Column('RA', u.deg, 'Right ascension in degrees')
        self.cols.all.Dec = Column('Dec', u.deg, 'Declination in degrees')
        self.cols.all.B = Column('B', u.mag, 'Apparent B magnitude')
        self.cols.all.B_err = Column('B_err', u.mag, 'Absolute error of apparent B magnitude')
        self.cols.all.B_flag = Column('B_flag', u.mag, '0: the B magnitude is measured; 1: the B magnitude is calculated from the B_J magnitude')
        self.cols.all.B_Abs = Column('B_Abs', u.mag, 'Absolute B magnitude')
        self.cols.all.J = Column('J', u.mag, 'Apparent J magnitude')
        self.cols.all.J_err = Column('J_err', u.mag, 'Absolute error of apparent J magnitude')
        self.cols.all.H = Column('H', u.mag, 'Apparent H magnitude')
        self.cols.all.H_err = Column('H_err', u.mag, 'Absolute error of apparent H magnitude')
        self.cols.all.K = Column('K', u.mag, 'Apparent K_s magnitude')
        self.cols.all.K_err = Column('K_err', u.mag, 'Absolute error of apparent K_s magnitude')
        self.cols.all.W1 = Column('W1', u.mag, 'Apparent W1 magnitude')
        self.cols.all.W1_err = Column('W1_err', u.mag, 'Absolute error of apparent W1 magnitude')
        self.cols.all.W2 = Column('W2', u.mag, 'Apparent W2 magnitude')
        self.cols.all.W2_err = Column('W2_err', u.mag, 'Absolute error of apparent W2 magnitude')
        self.cols.all.W1_flag = Column('W1_flag', u.mag, '0: the W1 magnitude is measured; 1: the W1 magnitude is calculated from the K_s magnitude')
        self.cols.all.B_J = Column('B_J', u.mag, 'Apparent B_J magnitude')
        self.cols.all.B_J_err = Column('B_J_err', u.mag, 'Absolute error of apparent B_J magnitude')
        self.cols.all.z_helio = Column('z_helio', None, 'Redshift in the heliocentric frame')
        self.cols.all.z_cmb = Column('z_cmb', None, 'Redshift converted to the Cosmic Microwave Background (CMB) frame')
        self.cols.all.z_flag = Column('z_flag', None, '0: the CMB frame redshift and luminosity distance values given '
                                                      'in columns 25 and 28 are not corrected for the peculiar velocity; 1: they are corrected values')
        self.cols.all.v_err = Column('v_err', None, 'Error of redshift from the peculiar velocity estimation')
        self.cols.all.z_err = Column('z_err', None, 'Measurement error of heliocentric redshift')
        self.cols.all.d_L = Column('d_L', u.Mpc, 'Luminosity distance in Mpc units')
        self.cols.all.d_L_err = Column('d_L_err', u.Mpc, 'Error of luminosity distance in Mpc units')
        self.cols.all.dist_flag = Column('dist_flag', None, '0: the galaxy has no measured redshift or distance value; 1: it has a measured photometric redshift from which we have calculated its luminosity distance; 2: it has a measured luminosity distance value from which we have calculated its redshift; 3: it has a measured spectroscopic redshift from which we have calculated its luminosity distance')
        self.cols.all.mstar = Column('mstar', 1e10 * u.Msun, 'Stellar mass in 10^10 M_Sun units')
        self.cols.all.mstar_err = Column('mstar_err', 1e10 * u.Msun, 'Absolute error of stellar mass in 10^10 M_Sun units')
        self.cols.all.merger_rate = Column('merger_rate', u.Gyr**-1, 'Base-10 logarithm of estimated BNS merger rate in the galaxy in Gyr^-1 units')
        self.cols.all.merger_rate_error = Column('merger_rate_error', u.Gyr**-1, 'Absolute error of estimated BNS merger rate in the galaxy')
