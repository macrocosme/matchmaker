import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
import pandas as pd
import sqlite3
import requests

class Grb(Catalog):
    file_location = DATA_BASE_PATH + 'data/grb_web/GRBweb2.sqlite' # limited to 15mag
    name = 'grb'

    def __init__(self, load_data=False, query=False):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))

        self.set_survey_specific_columns()

        self.cols.dec = self.cols.grb_web.decl
        self.cols.ra_err = self.cols.grb_web.pos_error
        self.cols.dec_err = self.cols.grb_web.pos_error

        self.cols.name = self.cols.grb_web.GRB_name
        # self.cols.obj_type = None

        if load_data:
            self.load_data(query=query)

    def load_data(self, query=False):
        # Code from https://user-web.icecube.wisc.edu/~grbweb_public/Loading_data.html

        # Download the SQLite file from the GRBweb webpage
        if query:
            r = requests.get("https://icecube.wisc.edu/~grbweb_public/GRBweb2.sqlite")
            f = open(self.file_location, 'wb').write(r.content)

        # Load the database with the sqlite3 module
        db = sqlite3.connect('GRBweb2.sqlite' if query else self.file_location)

        # Print the names of all the tables
        df_table_names = pd.read_sql_query("SELECT * from sqlite_sequence", db)
        # print("Table names:\n", df_table_names, "\n\n")

        # Let's use the summary table.
        # see the df_table_names example above to seek for individual sources.
        self.df = pd.read_sql_query("SELECT * from Summary", db)

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.grb_web = Column()
        self.cols.grb_web.id = Column('id')
        self.cols.grb_web.GRB_name = Column('GRB_name')
        self.cols.grb_web.GRB_name_Fermi = Column('GRB_name_Fermi')
        self.cols.grb_web.T0 = Column('T0', None, 'UTC')
        self.cols.grb_web.T0_source = Column('T0_source')
        self.cols.grb_web.ra = Column('ra', u.deg, 'J2000')
        self.cols.grb_web.ra_source = Column('ra_source', u.deg, 'J2000')
        self.cols.grb_web.decl = Column('decl', u.deg, 'J2000')
        self.cols.grb_web.decl_source = Column('decl_source', u.deg, 'J2000')
        self.cols.grb_web.pos_error = Column('pos_error', u.deg, '1-sigma')
        self.cols.grb_web.pos_error_source = Column('pos_error_source',)
        self.cols.grb_web.T90 = Column('T90', u.second)
        self.cols.grb_web.T90_source = Column('T90_source', u.second)
        self.cols.grb_web.T90_error = Column('T90_error', u.second)
        self.cols.grb_web.T90_error_source = Column('T90_error_source')
        self.cols.grb_web.T90_start = Column('T90_start', None, 'UTC')
        self.cols.grb_web.T90_start_source = Column('T90_start_source')
        self.cols.grb_web.fluence = Column('fluence', u.erg * u.cm**-2)
        self.cols.grb_web.fluence_source = Column('fluence_source', u.erg * u.cm**-2)
        self.cols.grb_web.fluence_error = Column('fluence_error', u.erg * u.cm**-2)
        self.cols.grb_web.fluence_error_source = Column('fluence_error_source', u.erg * u.cm**-2)
        self.cols.grb_web.redshift = Column('redshift')
        self.cols.grb_web.redshift_source = Column('redshift_source')
        self.cols.grb_web.T100 = Column('T100', u.second)
        self.cols.grb_web.GBM_located = Column('GBM_located')
        self.cols.grb_web.mjd = Column('mjd', u.day)
        self.cols.grb_web.mjd_source = Column('mjd_source', u.day)







