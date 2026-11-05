"""RACS-low, the Rapid ASKAP Continuum Survey at 887.5 MHz.

The expected file is the RACS-low DR1 source catalogue (Hale et al. 2021, PASA 38, e058):
2,123,638 sources over declinations -80 to +30 degrees, with PyBDSF column names. The centre
frequency is hardcoded here because the file does not carry it, so `load_data` checks that
fingerprint and refuses anything else. RACS-mid (Duchesne et al. 2023) is the likely mistake:
it reaches declination +49 with about 3.1 million sources and a centre of 1367.5 MHz, which
would corrupt every spectral index computed against LoTSS.
"""
import os

import astropy.units as u

from ..utils import load_fits_as_dataframe
from . import DATA_BASE_PATH, Catalog, Column

# Fingerprint of the RACS-low DR1 source catalogue.
EXPECTED_ROWS = 2123638
EXPECTED_DEC_MAX = 30.0
EXPECTED_DEC_MIN = -85.0


class Racs(Catalog):
    file_location = DATA_BASE_PATH + 'racs/racs.fits'
    name = 'racs'
    boxes = None

    def __init__(self, load_data=True, constrain=False, check_release=True):
        # Required and practical columns
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))

        self.set_survey_specific_columns()
        # McConnell et al. (2020), PASA 37, e048
        self.central_frequency = 887.5 * u.MHz
        self.bandwidth = 288 * u.MHz
        self.check_release = check_release

        if load_data:
            self.load_data(constrain=constrain)

    def load_data(self, constrain=False):
        if not os.path.exists(self.file_location):
            raise FileNotFoundError(
                f'RACS catalogue not found at {self.file_location}. The expected file is '
                'the RACS-low DR1 source catalogue (Hale+2021), also kept in this repository at '
                'archive/reproduction_package/data/racs.fits. Point MATCHMAKER_DATA_PATH at a '
                'directory holding racs/racs.fits, or set Racs.file_location directly.')
        self.df = load_fits_as_dataframe(self.file_location)
        self.n_source = len(self.df)
        if self.check_release:
            self.verify_release()

    def verify_release(self):
        """Refuse a file that is not RACS-low DR1, since the frequency is hardcoded for it."""
        problems = []
        if abs(self.n_source - EXPECTED_ROWS) > 0.02 * EXPECTED_ROWS:
            problems.append(f'{self.n_source:,} sources, expected about {EXPECTED_ROWS:,}')
        dec = self.df[self.cols.dec.label]
        if dec.max() > EXPECTED_DEC_MAX + 1 or dec.min() < EXPECTED_DEC_MIN - 1:
            problems.append(f'declinations {dec.min():.1f} to {dec.max():.1f}, expected '
                            f'{EXPECTED_DEC_MIN:.0f} to {EXPECTED_DEC_MAX:.0f}')
        missing = [c for c in ('total_flux_source', 'ra', 'e_ra') if c not in self.df.columns]
        if missing:
            problems.append(f'missing columns {missing}')
        if problems:
            raise ValueError(
                'This does not look like the RACS-low DR1 source catalogue: '
                + '; '.join(problems)
                + f'. central_frequency is hardcoded to {self.central_frequency}; a different '
                  'release needs a different value (RACS-mid 1367.5 MHz, RACS-high 1655.5 MHz). '
                  'Pass check_release=False only if you have set central_frequency yourself.')

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.ra = Column('ra', u.deg, 'RA J2000 in degrees')
        self.cols.dec = Column('dec', u.deg, 'DEC J2000 in degrees')
        self.cols.ra_err = Column('e_ra', u.arcsec, 'RA uncertainty in arcsec')
        self.cols.dec_err = Column('e_dec', u.arcsec, 'DEC uncertainty in arcsec')
        self.cols.measure = Column()
        # RACS-low DR1 source catalogue (Hale+2021) column names
        self.cols.measure.total_flux = Column('total_flux_source', u.mJy)
        self.cols.measure.total_flux_err = Column('e_total_flux_source', u.mJy)
        self.cols.measure.peak_flux = Column('peak_flux', u.mJy * u.beam**-1)
        self.cols.measure.peak_flux_err = Column('e_peak_flux', u.mJy * u.beam**-1)
