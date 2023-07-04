import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord

from . import (Catalog, Column, DATA_BASE_PATH)


class Lbcs(Catalog):
    file_location = DATA_BASE_PATH + 'lbcs/lbcs_stats.cons'
    name = 'lbcs'
    boxes = None

    def __init__(self, load_data=True):
        # Required and practical columns
        super().__init__(ra=Column('ra', u.deg), dec=Column('dec', u.deg))

        self.set_survey_specific_columns()

        if load_data:
            self.load_data()

    def load_data(self, minFlux=1.0):
        with open(self.file_location, 'rt') as f:
            #0) Observation number
            #1) RA of the pointing (in most cases determined from the WENSS catalogue,
                #   from which the sources are mainly selected)
            #2) Dec
            #3) and 4) Observation date and time
            #5) Goodness of the calibrator on the nine international stations:
                #   where "P" means a clear detection of fringes, "S" a marginal detection,
                #   "X" no detection, and "-" problems with the antenna.
            #6) Quality of the observation, judged by the percentage of source detections
                #   on the shorter baselines. Above 20 is probably fine, below 10 is almost
                #   certainly a bad observation.

            observation_numbers = []
            coordStrings = []
            ras = []
            decs = []
            observation_dates = []
            observation_times = []
            goodness_on_nine_stations = []
            quality_of_observations = []

            for line in f:
                if not line.startswith('#'):
                    cols = line.split()
                    # Get individual values
                    if 'P' in cols[5] or 'S' in cols[5]:
                        h, m, s = cols[1].split(':')
                        d, _m, _s = cols[2].split(':')
                        coord = SkyCoord("{}h{}m{}s {}d{}m{}s".format(h,m,s, d, _m, _s))

                        # Add to lists
                        observation_numbers.append(cols[0])
                        coordStrings.append(coord)
                        ras.append(coord.ra.deg)
                        decs.append(coord.dec.deg)
                        observation_dates.append(cols[3])
                        observation_times.append(cols[4])
                        goodness_on_nine_stations.append(cols[5])
                        quality_of_observations.append(cols[6])
                        # There are more columns but not sure what they are

            self.df = pd.DataFrame(
                {'observation_number': observation_numbers,
                 'coordString': coordStrings,
                 'ra': ras,
                 'dec': decs,
                 'observation_date': observation_dates,
                 'observation_time': observation_times,
                 'goodness_on_nine_station': goodness_on_nine_stations,
                 'quality_of_observation': quality_of_observations}
            )

    def set_survey_specific_columns(self):
        self.cols.observation_number = Column('observation_number', None, None)

        self.cols.ra = Column('ra', u.deg, 'Right ascension (J2000.0)')
        self.cols.dec = Column('dec', u.deg, 'Declination (J2000.0)')

        self.cols.coords = Column('coords', None, 'RA Dec (hms dms)')

        self.observation_date = Column('observation_date', None, 'Date of observation')
        self.observation_time = Column('observation_time', None, 'Time of observation')

        self.goodness_on_nine_station = Column('goodness_on_nine_station', None, """Goodness of the calibrator on the nine international stations:  
        where "P" means a clear detection of fringes, "S" a marginal detection, "X" no detection, and "-" problems with the antenna.""")
        self.quality_of_observation = Column('quality_of_observation', None, """Quality of the observation, judged by the percentage of source detections
        on the shorter baselines. Above 20 is probably fine, below 10 is almost certainly a bad observation.""")
