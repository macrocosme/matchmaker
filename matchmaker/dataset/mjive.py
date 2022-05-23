import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from ..utils import load_apj_table
from astropy.coordinates import SkyCoord

class Mjive(Catalog):
    file_location = DATA_BASE_PATH + 'data/mJIVE/MJIVE20-31MAR2014.txt'
    name = 'mjive'

    def __init__(self, load_data=False, constrain=True):
        super().__init__(ra=Column('RA', u.deg), dec=Column('DEC', u.deg))
        self.frequency_band = None
        self.central_frequency = 5 * u.MHz

        self.set_survey_specific_columns()
        self.cols.measure = Column()

        self.cols.major = self.cols.survey.dmaj
        self.cols.minor = self.cols.survey.dmin
        self.cols.pa = self.cols.survey.dpa

        self.cols.measure.total_flux = self.cols.survey.vint
        self.cols.measure.peak_flux = self.cols.survey.vpk

        if load_data:
            self.load_data()

    def load_data(self):
        df = load_apj_table(self.file_location)

        # Remove lines without VLBI information
        df.dropna(subset = ["VRAh"], inplace=True)

        # need to convert ra and dec input to deg.
        df['RA'] = df.apply(lambda row: SkyCoord("{} {} {} {}{} {} {}".format(row['VRAh'], row['VRAm'], row['VRAs'],
                                                                              row['VDE-'], row['VDEd'], row['VDEm'], row['VDEs']),
                                                 unit=(u.hourangle, u.deg)).ra.deg,
                            axis=1)

        df['DEC'] = df.apply(lambda row: SkyCoord("{} {} {} {}{} {} {}".format(row['VRAh'], row['VRAm'], row['VRAs'],
                                                                               row['VDE-'], row['VDEd'], row['VDEm'], row['VDEs']),
                                                  unit=(u.hourangle, u.deg)).dec.deg, axis=1)

        self.df = df

    # Should simply fetch the same code from the Lotss class as it is identical.
    def semi_major(self, mask=None, with_unit=False, to_unit=None):
        a = self.prop_to_unit('major', self.cols.major.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                a = (a * self.cols.major.unit).to(to_unit).value
            else:
                a = a.to(to_unit)

        return a

    # idem
    def semi_minor(self, mask=None, with_unit=False, to_unit=None):
        b = self.prop_to_unit('minor', self.cols.minor.unit, mask=mask, with_unit=with_unit) / 2

        if to_unit:
            if not with_unit:
                b = (b * self.cols.minor.unit).to(to_unit).value
            else:
                b = b.to(to_unit)

        return b

    def set_survey_specific_columns(self):
        self.cols.survey = Column()
        self.cols.survey.mid = Column('mID', None, 'The mJIVE-20 designation')
        self.cols.survey.frah = Column('FRAh', u.h, 'FIRST Hour of Right Ascension (J2000)')
        self.cols.survey.fram = Column('FRAm', u.min, 'FIRST Minute of Right Ascension (J2000)')
        self.cols.survey.fras = Column('FRAs', u.s, 'FIRST Second of Right Ascension (J2000)')
        self.cols.survey.fde = Column('FDE', None, 'FIRST Sign of the Declination (J2000)')
        self.cols.survey.fded = Column('FDEd', u.deg, 'FIRST Degree of Declination (J2000)')
        self.cols.survey.fdem = Column('FDEm', u.arcmin, 'FIRST Arcminute of Declination (J2000)')
        self.cols.survey.fdes = Column('FDEs', u.arcsec, 'FIRST Arcsecond of Declination (J2000)')
        self.cols.survey.fpk = Column('FPk', None, 'FIRST peak flux density; mJy/beam')
        self.cols.survey.fint = Column('FInt', u.mJy, 'FIRST integrated flux density')
        self.cols.survey.bmaj = Column('BMaj', u.mas, 'VLBI synthesized beam major axis')
        self.cols.survey.bmin = Column('BMin', u.mas, 'VLBI synthesized beam minor axis')
        self.cols.survey.bpa = Column('BPA', u.deg, 'VLBI synthesized beam position angle')
        self.cols.survey.l_vpk = Column('l_VPk', None, 'Limit flag on VPk')
        self.cols.survey.vpk = Column('VPk', None, 'VLBI peak flux density; mJy/beam')
        self.cols.survey.e_vpk = Column('e_VPk', None, 'Error in VPk (1)')
        self.cols.survey.vrah = Column('VRAh', u.hour, 'VLBI Hour of Right Ascension (J2000) (1)')
        self.cols.survey.vram = Column('VRAm', u.min, 'VLBI Minute of Right Ascension (J2000) (1)')
        self.cols.survey.vras = Column('VRAs', u.s, 'VLBI Second of Right Ascension (J2000) (1)')
        self.cols.survey.e_vras = Column('e_VRAs', u.s, 'Error in VRAs (1)')
        self.cols.survey.vde = Column('VDE', None, 'VLBI Sign of the Declination (J2000) (1)')
        self.cols.survey.vded = Column('VDEd', u.deg, 'VLBI Degree of Declination (J2000) (1)')
        self.cols.survey.vdem = Column('VDEm', u.arcmin, 'VLBI Arcminute of Declination (J2000) (1)')
        self.cols.survey.vdes = Column('VDEs', u.arcsec, 'VLBI Arcsecond of Declination (J2000) (1)')
        self.cols.survey.e_vdes = Column('e_VDEs', u.arcsec, 'Error in VDEs (1)')
        self.cols.survey.pkratio = Column('PkRatio', None, 'Ratio of VLBI to FIRST peak flux densities (1)')
        self.cols.survey.vint = Column('VInt', u.mJy, 'VLBI integrated flux density (1)')
        self.cols.survey.e_vint = Column('e_VInt', u.mJy, 'Error in VInt (1)')
        self.cols.survey.inratio = Column('InRatio', None, 'Ratio of VLBI integrated flux density to FIRST peak flux density (1)')
        self.cols.survey.dmaj = Column('DMaj', u.mas, 'VLBI single-gaussian deconvolved major axis (1)')
        self.cols.survey.e_dmaj = Column('e_DMaj', u.mas, 'Error in DMaj (1)')
        self.cols.survey.dmin = Column('DMin', u.mas, 'VLBI single-gaussian deconvolved minor axis (1)')
        self.cols.survey.e_dmin = Column('e_DMin', u.mas, 'Error in DMin (1)')
        self.cols.survey.dpa = Column('DPA', u.deg, 'VLBI single-gaussian deconvolved position angle (1)')
        self.cols.survey.e_dpa = Column('e_DPA', u.deg, 'Error in DPA (1)')
        self.cols.survey.cplx = Column('CPLX', None, 'Complex flag (1)')
