import astropy.units as u
from . import (Catalog, Column, DATA_BASE_PATH)
from pandas import read_csv
from astropy.coordinates import SkyCoord
from .tns import get_results_as_dataframe

class TnsSN(Catalog):
    # file_location = DATA_BASE_PATH + 'data/tns-sn/tns_search-16mag.csv' # limited to 15mag
    file_location = DATA_BASE_PATH + 'data/tns-sn/tns_sn.csv'
    name = 'tns-sn'

    def __init__(self, load_data=False, query=False, save_csv=False):
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))

        self.set_survey_specific_columns()
        self.cols.name = self.cols.tns.name
        self.cols.obj_type = self.cols.tns.obj_type

        if load_data:
            self.load_data(query=query, save_csv=save_csv)

    def load_data(self, query=False, save_csv=False):
        if not query:
            df = read_csv(self.file_location)
            df['ra'] = df.apply(lambda row: SkyCoord('{}h{}m{}s'.format(row['RA'].split(':')[0], row['RA'].split(':')[1], row['RA'].split(':')[2]),
                                               '{}d{}m{}s'.format(row['DEC'].split(':')[0], row['DEC'].split(':')[1], row['DEC'].split(':')[2]), frame='icrs').ra.deg, axis=1)
            df['dec'] = df.apply(lambda row: SkyCoord('{}h{}m{}s'.format(row['RA'].split(':')[0], row['RA'].split(':')[1], row['RA'].split(':')[2]),
                                                        '{}d{}m{}s'.format(row['DEC'].split(':')[0], row['DEC'].split(':')[1], row['DEC'].split(':')[2]), frame='icrs').dec.deg, axis=1)
        else:
            df = get_results_as_dataframe()
            if save_csv:
                df.to_csv(self.file_location, index=False)

        self.df = df

    def set_survey_specific_columns(self):
        # Radio columns
        self.cols.tns = Column()
        self.cols.tns.id = Column('ID')
        self.cols.tns.name = Column('Name')
        self.cols.tns.raj = Column('RA')
        self.cols.tns.decj = Column('DEC')
        self.cols.tns.obj_type = Column('Obj. Type')
        self.cols.tns.redshift = Column('Redshift')
        self.cols.tns.host_name = Column('Host Name')
        self.cols.tns.host_redshift = Column('Host Redshift')
        self.cols.tns.reporting_group_s = Column('Reporting Group/s')
        self.cols.tns.discovery_data_source_s = Column('Discovery Data Source/s')
        self.cols.tns.classifying_group_s = Column('Classifying Group/s')
        self.cols.tns.associated_group_s = Column('Associated Group/s')
        self.cols.tns.disc_internal_name = Column('Disc. Internal Name')
        self.cols.tns.disc_instrument_s = Column('Disc. Instrument/s')
        self.cols.tns.class_instrument_s = Column('Class. Instrument/s')
        self.cols.tns.tns_at = Column('TNS AT')
        self.cols.tns.public = Column('Public')
        self.cols.tns.end_prop_period = Column('End Prop. Period')
        self.cols.tns.discovery_mag_flux = Column('Discovery Mag/Flux')
        self.cols.tns.discovery_filter = Column('Discovery Filter')
        self.cols.tns.discovery_date_ut = Column('Discovery Date (UT)')
        self.cols.tns.sender = Column('Sender')
        self.cols.tns.remarks = Column('Remarks')
        self.cols.tns.ext_catalog_s = Column('Ext. catalog/s')
