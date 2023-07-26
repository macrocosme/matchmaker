import numpy as np
import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord

from . import (Catalog, Column, DATA_BASE_PATH)


class Rfc(Catalog):
    file_location = DATA_BASE_PATH + 'rfc/rfc_2023a_cat.txt'
    name = 'rfc'
    boxes = None

    def __init__(self, load_data=True, minFlux=1.0):
        # Required and practical columns
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))

        self.set_survey_specific_columns()

        self.total_flux = self.cols.flux.unresolved_flux_S_band

        if load_data:
            self.load_data(minFlux=minFlux)

    def load_data(self, minFlux=1.0):
        # Adapted from Benito Marcote's Seffers sources.py (https://github.com/bmarcote/seffers)
        with open(self.file_location, 'rt') as f:
            source_names = []
            ivs_names = []
            categories = []
            coordStrings = []
            ras = []
            decs = []
            ra_errs = []
            dec_errs = []
            resolved_fluxes_S = []
            unresolved_fluxes_S = []
            resolved_fluxes_C = []
            unresolved_fluxes_C = []
            resolved_fluxes_X = []
            unresolved_fluxes_X = []
            resolved_fluxes_U = []
            unresolved_fluxes_U = []
            resolved_fluxes_K = []
            unresolved_fluxes_K = []

            for line in f:
                if not (line.startswith('#') or line.startswith('U')): # or line.startswith('N'):
                    cols = line.split()
                    cols[13:23] = [float(f) if '<' not in f else 0.0 for f in cols[13:23]]
                    if cols[16] > minFlux:
                        # Get individual values
                        source_name = cols[2]
                        ivs_name = cols[1]
                        category = cols[0]
                        coord = SkyCoord("{}h{}m{}s {}d{}m{}s".format(*cols[3:9]))
                        ra_err = cols[9]
                        dec_err = cols[10]
                        resolved_flux_S, unresolved_flux_S = cols[13], cols[14]
                        resolved_flux_C, unresolved_flux_C = cols[15], cols[16]
                        resolved_flux_X, unresolved_flux_X = cols[17], cols[18]
                        resolved_flux_U, unresolved_flux_U = cols[19], cols[20]
                        resolved_flux_K, unresolved_flux_K = cols[21], cols[22]

                        # Add to lists
                        source_names.append(source_name)
                        ivs_names.append(ivs_name)
                        categories.append(category)
                        coordStrings.append(coord)
                        ras.append(coord.ra.deg)
                        decs.append(coord.dec.deg)
                        ra_errs.append(ra_err)
                        dec_errs.append(dec_err)
                        resolved_fluxes_S.append(resolved_flux_S)
                        unresolved_fluxes_S.append(unresolved_flux_S)
                        resolved_fluxes_C.append(resolved_flux_C)
                        unresolved_fluxes_C.append(unresolved_flux_C)
                        resolved_fluxes_X.append(resolved_flux_X)
                        unresolved_fluxes_X.append(unresolved_flux_X)
                        resolved_fluxes_U.append(resolved_flux_U)
                        unresolved_fluxes_U.append(unresolved_flux_U)
                        resolved_fluxes_K.append(resolved_flux_K)
                        unresolved_fluxes_K.append(unresolved_flux_K)

            self.df = pd.DataFrame(
                {'source_name': source_names,
                 'ivs_name': ivs_names,
                 'categorie': categories,
                 'coordString': coordStrings,
                 'ra': ras,
                 'dec': decs,
                 'ra_err': ra_errs,
                 'dec_err': dec_errs,
                 'resolved_flux_S_band': resolved_fluxes_S,
                 'unresolved_flux_S_band': unresolved_fluxes_S,
                 'resolved_flux_C_band': resolved_fluxes_C,
                 'unresolved_flux_C_band': unresolved_fluxes_C,
                 'resolved_flux_X_band': resolved_fluxes_X,
                 'unresolved_flux_X_band': unresolved_fluxes_X,
                 'resolved_flux_U_band': resolved_fluxes_U,
                 'unresolved_flux_U_band': unresolved_fluxes_U,
                 'resolved_flux_K_band': resolved_fluxes_K,
                 'unresolved_flux_K_band': unresolved_fluxes_K}
            )

    def set_survey_specific_columns(self):
        self.cols.source_name = Column('Source_name', None, 'IAU name (J2000.0)')
        self.cols.ivs_name = Column('ivs_name', None, 'IVS name (B1950)')

        self.cols.category = Column('category', None, 'Category: C (calibrator), N (non-calibrator), U (unreliable coordinates)')

        self.cols.ra = Column('ra', u.deg, 'Right ascension (J2000.0)')
        self.cols.dec = Column('dec', u.deg, 'Declination (J2000.0)')

        self.cols.coords = Column('coords', None, 'RA Dec (hms dms)')
        self.cols.ra_err = Column('ra_err', u.mas, 'Right ascension inflated error (mas)')
        self.cols.dec_err = Column('dec_err', u.mas, 'Declination inflated error (mas)')

        self.cols.flux = Column()
        # S band
        self.cols.flux.resolved_flux_S_band = Column('resolved_flux_S_band', u.Jy, 'S-band total flux density integrated over entire map,  Jy')
        self.cols.flux.unresolved_flux_S_band = Column('unresolved_flux_S_band', u.Jy, 'S-band unresolved flux density at long VLBA baselines+,  Jy')
        # C band
        self.cols.flux.resolved_flux_C_band = Column('resolved_flux_C_band', u.Jy, 'C-band total flux density integrated over entire map,  Jy')
        self.cols.flux.unresolved_flux_C_band = Column('unresolved_flux_C_band', u.Jy, 'C-band unresolved flux density at long VLBA baselines+,  Jy')
        # X band
        self.cols.flux.resolved_flux_X_band = Column('resolved_flux_X_band', u.Jy, 'X-band total flux density integrated over entire map,  Jy')
        self.cols.flux.unresolved_flux_X_band = Column('unresolved_flux_X_band', u.Jy, 'X-band unresolved flux density at long VLBA baselines+,  Jy')
        # U band
        self.cols.flux.resolved_flux_U_band = Column('resolved_flux_U_band', u.Jy, 'U-band total flux density integrated over entire map,  Jy')
        self.cols.flux.unresolved_flux_U_band = Column('unresolved_flux_U_band', u.Jy, 'U-band unresolved flux density at long VLBA baselines+,  Jy')
        # K band
        self.cols.flux.resolved_flux_K_band = Column('resolved_flux_K_band', u.Jy, 'K-band total flux density integrated over entire map,  Jy')
        self.cols.flux.unresolved_flux_K_band = Column('unresolved_flux_K_band', u.Jy, 'K-band unresolved flux density at long VLBA baselines+,  Jy')
