"""dataset.__init__.py
"""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

# todo: set BASEPATH somewhere else
DATA_BASE_PATH = '/Users/vohl/Documents/code/EOSC/'

def box2d(obj, box_xxyy):
    df = obj.df.loc[
        (obj.df[obj.cols.ra.label] >= box_xxyy[0]) &
        (obj.df[obj.cols.ra.label] <= box_xxyy[1]) &
        (obj.df[obj.cols.dec.label] >= box_xxyy[2]) &
        (obj.df[obj.cols.dec.label] <= box_xxyy[3])
    ]
    return df

def filter_dataframe_to_box(obj, boxes):
    df_boxed = pd.concat([box2d(obj, b) for b in boxes])
    df_boxed = df_boxed.reset_index().rename(
        {'index': 'index_df_{}'.format(obj.name)},
        axis = 'columns'
    )
    return df_boxed

class Column:
    def __init__(self, label=None, unit=None, description=None):
        self.label = label
        self.unit = unit
        self.description = description

    def list_columns(self):
        # This could be done better. Currently, it only prints the first level.
        for k, _ in self.__dict__.items():
            print(k)

class Columns:
    def __init__(self, ra=None, dec=None):
        self.ra = ra
        self.dec = dec

    def prop(self, prop):
        return getattr(self, prop)

    def list_columns(self):
        # This could be done better. Currently, it only prints the first level.
        for k, _ in self.__dict__.items():
            print(k)

class File:
    def __init__(self, location=None, local:str=None, remote=None):
        if location is not None:
            self.location = location
        else:
            self.location = Column('location')

        if local is not None:
            self.location.local = local

        if remote is not None:
            self.location.remote = remote

class Image:
    def __init__(self, source_catalog=None):
        self.source_catalog = source_catalog
        self.filters = {}

class Chance:
    def __init__(self, total, expected):
        self.total = total
        self.expected = expected

class Match:
    def __init__(self, catalog_name, from_to=0, mask=None, idx=None, sep2d=None, sep3d=None):
        self.catalog_name = catalog_name
        self.from_to = from_to
        self._chance = None
        self.images = {}

        self.further_filters = {}

        if mask is not None:
            self.mask = mask
        if idx is not None:
            self.idx = idx
        if sep2d is not None:
            self.sep2d = sep2d
        if sep3d is not None:
            self.sep3d = sep3d

    @property
    def chance(self):
        return self._chance

    def set_chance(self, total=None, expected=None):
        self._chance = Chance(total, expected)

    @property
    def filtered_idx(self):
        return self.idx[self.mask]

    @property
    def anti_filtered_idx(self):
        return self.idx[~self.mask]

    @property
    def mask_idx(self):
        return np.where(self.mask)[0]

    @property
    def anti_mask_idx(self):
        return np.where(~self.mask)[0]

    @property
    def filtered_sep2d(self):
        return self.sep2d[self.mask]

    # @property
    # def filtered_sep3d(self):
    #     return self.sep3d # already includes mask filtering (see lotss.py)
    #                       # writing this comment highlights that sep3d could be handled differently
    #                       # (what if it's not from lotss?)

class Catalog:
    df = None
    name = None

    def __init__(self, ra=None, dec=None, use_distance=False):
        self.cols = Columns(ra=ra, dec=dec)
        self.use_distance = use_distance
        self.matches = {}

    def add_match(self, catalog_name, from_to=0, mask=None, idx=None, sep2d=None, sep3d=None):
        self.matches[catalog_name] = Match(catalog_name=catalog_name,
                                           from_to=from_to,
                                           mask=mask,
                                           idx=idx,
                                           sep2d=sep2d,
                                           sep3d=sep3d)

    def load_data(self):
        pass

    def filter_box(self, boxes):
        # Removes any object outside RA/DEC boxes
        # Filtering happens in place
        df = filter_dataframe_to_box(self, boxes)
        self.df = df

    def as_SkyCoord(self, mask=None, distance=None, distance_unit=None):
        if mask is None:
            df = self.df
            is_array = (len(df)>1 is np.ndarray)
        else:
            try:
                df = self.df.iloc[mask]
            except NotImplementedError:
                df = self.df.loc[mask]
            is_array = (type(mask) in [pd.Series, np.ndarray])

        ra_label = self.cols.ra.label
        dec_label = self.cols.dec.label

        ra = df[ra_label].values if is_array else df[ra_label]
        dec = df[dec_label].values if is_array else df[dec_label]

        if not self.use_distance and distance is None:
            sc = SkyCoord(ra=ra * self.cols.ra.unit,
                          dec=dec * self.cols.dec.unit)
        elif self.use_distance and distance is None:
            distance_label = self.cols.distance.label
            distance = df[distance_label].values if is_array else df[distance_label]
            sc = SkyCoord(ra=ra * self.cols.ra.unit,
                          dec=dec * self.cols.dec.unit,
                          distance=distance * self.cols.distance.unit)
        else:
            # Case where we want to specify what to use as distance based on previous cross-match
            # TODO: assert all is in order
            sc = SkyCoord(ra=df[ra_label].values * self.cols.ra.unit,
                          dec=df[dec_label].values * self.cols.dec.unit,
                          distance=distance * distance_unit)
            sc = None

        return sc

    def prop(self, prop, default=None):
        return getattr(self, prop, default)

    def prop_to_unit(self, prop, unit, mask=None, with_unit=False):
        if mask is None:
            try:
                converted_prop = (self.df[self.cols.prop(prop).label].values * self.cols.prop(prop).unit).to(unit)
            except:
                converted_prop = (self.df[self.cols.prop(prop).label] * self.cols.prop(prop).unit).to(unit)
        else:
            try:
                converted_prop = (self.df.iloc[mask][self.cols.prop(prop).label].values * self.cols.prop(prop).unit).to(unit)
            except AttributeError:
                converted_prop = (self.df.iloc[mask][self.cols.prop(prop).label] * self.cols.prop(prop).unit).to(unit)

        if with_unit:
            return converted_prop
        else:
            return converted_prop.value
