import numpy as np
import pandas as pd
import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from astroquery.gaia import Gaia

class Andrews(Catalog):
    """A collection of BH and NS in Gaia DR3"""
    name = 'andrews'

    def __init__(self, load_data=False, use_distance=False):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg), use_distance=use_distance)

        self.set_survey_specific_columns()

        self.cols.ra_err = self.cols.survey.ra_err
        self.cols.dec_err = self.cols.survey.dec_err
        if load_data:
            self.load_data()



    def load_data(self):
        df = Gaia.launch_job("""SELECT source_id, ref_epoch, ra, dec, ra_error, dec_error
                                FROM gaiadr3.gaia_source 
                                WHERE source_id in (5681911574178198400, 3649963989549165440, 747174436620510976, 
                                                    1581117310088807552, 1525829295599805184, 4271998639836225920, 
                                                    1695294922548180224, 1854241667792418304, 1058875159778407808, 
                                                    1947292821452944896, 2397135910639986304, 1144019690966028928, 
                                                    6593763230249162112, 5590962927271507712, 809741149368202752, 
                                                    5847919241396757888, 5580526947012630912, 1350295047363872512, 
                                                    4744087975990080896, 6001459821083925120, 1749013354127453696, 
                                                    4314242838679237120, 5593444799901901696, 6328149636482597888)
                                ORDER BY source_id""").get_results().to_pandas()
        self.df = df

    def set_survey_specific_columns(self):
        self.cols.survey = Column()
        self.cols.survey.source_id = Column('source_id', None, None)
        self.cols.survey.ref_epoch = Column('ref_epoch', None, None)
        self.cols.survey.ra = Column('ra', u.deg, None)
        self.cols.survey.dec = Column('dec', u.deg, None)
        self.cols.survey.ra_err = Column('ra_error', u.mas, None)
        self.cols.survey.dec_err = Column('dec_error', u.mas, None)


