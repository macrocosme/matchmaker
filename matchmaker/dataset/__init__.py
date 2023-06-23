"""dataset.__init__.py
"""
import astropy.units
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

#   
#  Set DATA_BASE_PATH in accordance to the location of your datasets
DATA_BASE_PATH = '../../data/'

class Column:
    """Column object"""
    def __init__(self, label:str=None, unit:astropy.units.Unit=None, description:str=None):
        """Initialize Column class

        Parameters
        ----------
        label:str  (Optional)
            Column label in the Pandas.DataFrame
        unit:str  (Optional)
            Unit of values in the column
        description:str  (Optional)
            Description of column's content
        """
        self.label = label
        self.unit = unit
        self.description = description

    def list_properties(self):
        """List all properties"""

        # This could be done better. Currently, it only prints the first level.
        for k, _ in self.__dict__.items():
            print(k)

class Columns:
    """Columns object: gathering of multiple Column objects"""
    def __init__(self, ra=None, dec=None):
        """Initialization

        Parameters
        ----------
        ra:float  (Optional)
            Right ascension (deg)
        dec:float  (Optional)
            Declination (deg)
        """
        self.ra = ra
        self.dec = dec

    def prop(self, prop:str):
        """Get a given property

        Parameter
        ---------
        prop:str
            Property variable name
        """
        return getattr(self, prop)

    def list_columns(self):
        """List all columns"""

        # This could be done better. Currently, it only prints the first level.
        for k, _ in self.__dict__.items():
            print(k)

class File:
    """File object"""
    def __init__(self, location:Column=None, local:str=None, remote=None):
        """Initialization"""
        if location is not None:
            self.location = location
        else:
            self.location = Column('location')

        if local is not None:
            self.location.local = local

        if remote is not None:
            self.location.remote = remote

class Image:
    """Image object"""
    def __init__(self, source_catalog=None):
        """Initialization

        Parameters
        ----------
        source_catalog:str  (optional)
        """
        self.source_catalog = source_catalog
        self.filters = {}

class Chance:
    """Chance association object"""
    def __init__(self, total, expected):
        """Initialization

        Parameters
        ----------
        total:int
            Total number of matches
        expected:int
            Expected number of chance association
        """
        self.total = total
        self.expected = expected

class Match:
    """Match object"""
    def __init__(self, catalog_name, from_to=0, mask=None, idx=None, sep2d=None, sep3d=None):
        """Initialization

        Parameters
        ----------
        catalog_name:str
            Name of the matched catalog
        from_to:int
            Indicates if the self is passed as obj1 or obj2 in matchmaker.matchmaker.crossmatch
            TODO: check which one was which
            TODO: remove this -- I don't think it is used anywhere, really
        mask:np.array
        idx:np.array
        sep2d:np.array of astropy.coordinates.angles.Angle
        sep3d:np.array of astropy.coordinates.distances.Distance
        """
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
        """chance property"""
        return self._chance

    def set_chance(self, total=None, expected=None):
        """Set chance property

        Parameters
        ----------
        total:int (optional)
        expected:int (optional)
        """
        self._chance = Chance(total, expected)

    @property
    def filtered_idx(self):
        """Filtered indices property

        Filtered indices are match.idx[mask]
        """
        return self.idx[self.mask]

    @property
    def anti_filtered_idx(self):
        """Anti-filtered indices property

        Anti-filtered indices are match.idx[~mask]
        """
        return self.idx[~self.mask]

    @property
    def mask_idx(self):
        """Indices where mask is True

        Returns
        -------
        array of indices: numpy.array
        """
        return np.where(self.mask)[0]

    @property
    def anti_mask_idx(self):
        """Indices where mask is False

        Returns
        -------
        array of indices: numpy.array
        """
        return np.where(~self.mask)[0]

    @property
    def filtered_sep2d(self):
        """Apply self.mask to sep2d"""
        return self.sep2d[self.mask]

    # @property
    # def filtered_sep3d(self):
    #     return self.sep3d # already includes mask filtering (see lotss.py)
    #                       # writing this comment highlights that sep3d could be handled differently
    #                       # (what if it's not from lotss?)

class Catalog:
    """Catalog object"""
    df = None
    name = None

    def __init__(self, ra=None, dec=None, use_distance=False):
        """Initialization

        ra:float   (optional)
            Right ascension in degree
        dec:float   (optional)
            Declination in degree
        use_distance:bool   (optional)
            Use distance (self.as_SkyCoord() will include distance if available)
        """
        self.cols = Columns(ra=ra, dec=dec)
        self.use_distance = use_distance
        self.matches = {}

    def add_match(self, catalog_name, from_to=0, mask=None, idx=None, sep2d=None, sep3d=None):
        """Add match from Catalog catalog_name

        Parameters
        ----------
        catalog_name:str
            Name of the matched catalog
        from_to:int   (optional)
            Indicates if the self is passed as obj1 or obj2 in matchmaker.matchmaker.crossmatch
            TODO: check which one was which
            TODO: remove this -- I don't think it is used anywhere, really
        mask:np.array  (optional)
        idx:np.array   (optional)
        sep2d:np.array of astropy.coordinates.angles.Angle   (optional)
        sep3d:np.array of astropy.coordinates.distances.Distance   (optional)
        """
        self.matches[catalog_name] = Match(catalog_name=catalog_name,
                                           from_to=from_to,
                                           mask=mask,
                                           idx=idx,
                                           sep2d=sep2d,
                                           sep3d=sep3d)

    def get_mask_to_third_catalog(self, cat1, cat2):
        """Get mask for a third catalog (cat2) given self and cat1

        Parameters
        ----------
        cat1:Catalog
        cat2:Catalog

        Returns
        -------
        cat2_mask:numpy.array
        """
        def to_array(cat1, cat2):
            return np.array([
                np.where(i == self.matches[cat1].mask_filtered[f'{cat2}_{self.name}'])[0][0] if \
                np.where(i == self.matches[cat1].mask_filtered[f'{cat2}_{self.name}'])[0].size > 0 else \
                False \
                for i in self.matches[cat1].mask_filtered[self.name]
            ])

        return to_array(cat1, cat2)[np.where(to_array(cat1, cat2))[0]]

    def load_data(self):
        """Load data (to be overriden by survey dataset class)"""
        pass

    def filter_box(self, boxes):
        """Removes any object outside RA/DEC boxes

        Filtering happens in place

        Parameter
        ---------
        boxes:list of boxes
        """
        df = filter_dataframe_to_box(self, boxes)
        self.df = df

    def as_SkyCoord(self, mask=None, distance=None, distance_unit=None):
        """Get source coordinates as astropy.coordinates.SkyCoord

        Parameters
        ----------
        mask:numpy.array   (optional)
            Index mask for rows to be returned in the dataframe
        distance:bool   (optional)
            Include distance in SkyCoord or not
        distance_unit   (optional)
            Distance unit to be returned if distance is True
        """
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
        """Get property by variable name"""
        return getattr(self, prop, default)

    def prop_to_unit(self, prop:str, unit:astropy.units.Unit, mask=None, with_unit=False):
        """Get property in given unit

        Properties
        ----------
        prop: str
        unit: astropy.units.Unit
        mask: numpy.array   (optional)
        with_unit: bool   (optional, default: False)
            Return property with unit or not

        Returns
        -------
        converted_prop
        """
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
            assert isinstance(converted_prop.value, object)
            return converted_prop.value

def box2d(obj:Catalog, box_xxyy:list):
    """Filter obj.df to box_xxyy limits

    Parameters
    ----------
    obj:Catalog
        Catalog object
    box_xxyy:list
        List of box corners. The list is of size 4,
        with indices 0 (min) and 1 (max) for ra and 2 (min) and 3 (max) for dec

    Returns
    -------
    filtered_df:Pandas.DataFrame
        The obj dataframe filtered to ra and dec inside the box_xxyy
    """
    df = obj.df.loc[
        (obj.df[obj.cols.ra.label] >= box_xxyy[0]) &
        (obj.df[obj.cols.ra.label] <= box_xxyy[1]) &
        (obj.df[obj.cols.dec.label] >= box_xxyy[2]) &
        (obj.df[obj.cols.dec.label] <= box_xxyy[3])
    ]
    return df

def filter_dataframe_to_box(obj:Catalog, boxes:list):
    """Filter obj.df to multiple boxes

    Parameters
    ----------
    obj:Catalog
        Catalog object
    boxes:list
        List of boxes lists. A given box list is of size 4,
        with indices 0 (min) and 1 (max) for ra and 2 (min) and 3 (max) for dec.

    Returns
    -------
    filtered_df:Pandas.DataFrame
        The obj dataframe filtered to ra and dec inside the boxes
    """
    df_boxed = pd.concat([box2d(obj, b) for b in boxes])
    df_boxed = df_boxed.reset_index().rename(
        {'index': 'index_df_{}'.format(obj.name)},
        axis = 'columns'
    )
    return df_boxed
