"""dataset.__init__.py
"""
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

#   
#  Set DATA_BASE_PATH in accordance to the location of your datasets
DATA_BASE_PATH = '../../data/'

def box2d(obj, box_xxyy):
    """
    Filter the DataFrame `obj.df` based on the given bounding box coordinates.

    Args:
        obj (object): The object containing the DataFrame to be filtered.
        box_xxyy (list): The bounding box coordinates in the format [xmin, xmax, ymin, ymax].

    Returns:
        pandas.DataFrame: The filtered DataFrame containing rows within the specified bounding box.
    """
    df = obj.df.loc[
        (obj.df[obj.cols.ra.label] >= box_xxyy[0]) &
        (obj.df[obj.cols.ra.label] <= box_xxyy[1]) &
        (obj.df[obj.cols.dec.label] >= box_xxyy[2]) &
        (obj.df[obj.cols.dec.label] <= box_xxyy[3])
    ]
    return df

def filter_dataframe_to_box(obj, boxes):
    """
    Filters a dataframe to include only the rows that fall within the specified boxes.

    Args:
        obj: The dataframe to filter.
        boxes: A list of boxes (bounding boxes) to filter the dataframe with.

    Returns:
        A new dataframe containing only the rows that fall within the specified boxes.
    """
    df_boxed = pd.concat([box2d(obj, b) for b in boxes])
    df_boxed = df_boxed.reset_index().rename(
        {'index': 'index_df_{}'.format(obj.name)},
        axis='columns'
    )
    return df_boxed

class Column:
    """
    Represents a column in a dataset.

    Attributes:
        label (str): The label of the column.
        unit (str): The unit of measurement for the column.
        description (str): A description of the column.
    """

    def __init__(self, label=None, unit=None, description=None):
            """
            Initialize a new instance of the Dataset class.

            Args:
                label (str): The label of the dataset.
                unit (str): The unit of measurement for the dataset.
                description (str): A description of the dataset.
            """
            self.label = label
            self.unit = unit
            self.description = description

    def list_columns(self):
        """
        Prints the names of the attributes of the Column object.

        Note:
            This method currently only prints the first level of attributes.

        Returns:
            None
        """
        for k, _ in self.__dict__.items():
            print(k)

class Columns:
    """
    A class representing columns in a dataset.

    Attributes:
        ra (str): The right ascension column.
        dec (str): The declination column.
    """

    def __init__(self, ra=None, dec=None):
        """
        Initialize a Dataset object.

        Args:
            ra (float): The right ascension value.
            dec (float): The declination value.
        """
        self.ra = ra
        self.dec = dec

    def prop(self, prop):
            """
            Returns the value of the specified property.

            Args:
                prop (str): The name of the property.

            Returns:
                Any: The value of the property.

            """
            return getattr(self, prop)

    def list_columns(self):
        """
        Prints the names of the columns in the dataset.

        Note:
            This method currently only prints the first level of columns.
        """
        for k, _ in self.__dict__.items():
            print(k)

class File:
    """
    Represents a file in the dataset.

    Args:
        location (str, optional): The location of the file. Defaults to None.
        local (str, optional): The local path of the file. Defaults to None.
        remote (str, optional): The remote path of the file. Defaults to None.
    """
    def __init__(self, location=None, local:str=None, remote=None):
        """
        Initialize a Dataset object.

        Args:
            location (str, optional): The location of the dataset. Defaults to None.
            local (str, optional): The local path of the dataset. Defaults to None.
            remote (str, optional): The remote path of the dataset. Defaults to None.
        """
        if location is not None:
            self.location = location
        else:
            self.location = Column('location')

        if local is not None:
            self.location.local = local

        if remote is not None:
            self.location.remote = remote

class Image:
    """
    Represents an image in the dataset.

    Attributes:
        source_catalog (str): The source catalog of the image.
        filters (dict): A dictionary of filters applied to the image.
    """

    def __init__(self, source_catalog=None):
            """
            Initialize a new instance of the Dataset class.

            Args:
                source_catalog (str): The source catalog for the dataset.

            Attributes:
                source_catalog (str): The source catalog for the dataset.
                filters (dict): A dictionary of filters for the dataset.

            """
            self.source_catalog = source_catalog
            self.filters = {}

class Chance:
    """
    Represents a chance object with total and expected values.

    Attributes:
        total (int): The total value.
        expected (float): The expected value.
    """

    def __init__(self, total, expected):
        """
        Initializes a new instance of the Dataset class.

        Args:
            total (int): The total number of items in the dataset.
            expected (int): The expected number of items in the dataset.

        Returns:
            None
        """
        self.total = total
        self.expected = expected

class Match:
    """
    Represents a match in the dataset.

    Args:
        catalog_name (str): The name of the catalog.
        from_to (int, optional): The range of matches. Defaults to 0.
        mask (ndarray, optional): The mask for filtering matches. Defaults to None.
        idx (ndarray, optional): The indices of the matches. Defaults to None.
        sep2d (ndarray, optional): The 2D separations of the matches. Defaults to None.
        sep3d (ndarray, optional): The 3D separations of the matches. Defaults to None.
    """

    def __init__(self, catalog_name, from_to=0, mask=None, idx=None, sep2d=None, sep3d=None):
        """
        Initialize a Dataset object.

        Parameters:
        - catalog_name (str): The name of the catalog.
        - from_to (int): The range of items to load from the catalog.
        - mask (optional): The mask to apply to the dataset.
        - idx (optional): The index of the dataset.
        - sep2d (optional): The 2D separation.
        - sep3d (optional): The 3D separation.
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
            """
            Returns the chance value associated with the dataset.

            Returns:
                float: The chance value of the dataset.
            """
            return self._chance

    def set_chance(self, total=None, expected=None):
            """
            Sets the chance for the dataset.

            Args:
                total (int, optional): The total number of items in the dataset. Defaults to None.
                expected (float, optional): The expected probability of an item being selected. Defaults to None.
            """
            self._chance = Chance(total, expected)

    @property
    def filtered_idx(self):
            """
            Returns the filtered indices based on the current mask.

            Returns:
                numpy.ndarray: The filtered indices.
            """
            return self.idx[self.mask]

    @property
    def anti_filtered_idx(self):
            """
            Returns the indices of the elements in the dataset that are not filtered out.

            Returns:
                numpy.ndarray: The indices of the elements that are not filtered out.
            """
            return self.idx[~self.mask]

    @property
    def mask_idx(self):
            """
            Returns the indices of the masked elements in the dataset.

            Returns:
                numpy.ndarray: An array containing the indices of the masked elements.
            """
            return np.where(self.mask)[0]

    @property
    def anti_mask_idx(self):
            """
            Returns the indices of elements in the dataset that are not masked.

            Returns:
                numpy.ndarray: An array containing the indices of unmasked elements.
            """
            return np.where(~self.mask)[0]

    @property
    def filtered_sep2d(self):
        """
        Returns a filtered version of the sep2d array based on the mask.

        The method filters the sep2d array using the mask attribute and returns the filtered version.

        Returns:
            numpy.ndarray: The filtered sep2d array.

        """
        return self.sep2d[self.mask]

    # @property
    # def filtered_sep3d(self):
    #     return self.sep3d # already includes mask filtering (see lotss.py)
    #                       # writing this comment highlights that sep3d could be handled differently
    #                       # (what if it's not from lotss?)

class Catalog:
    """
    Represents a catalog of astronomical objects.

    Attributes:
        df: The DataFrame containing the catalog data.
        name: The name of the catalog.
        cols: An instance of the Columns class representing the column names and units.
        use_distance: A boolean indicating whether distance information is used in the catalog.
        matches: A dictionary containing information about matches with other catalogs.
    """

    df = None
    name = None

    def __init__(self, ra=None, dec=None, use_distance=False):
        """
        Initializes a Catalog object.

        Args:
            ra: The right ascension column name.
            dec: The declination column name.
            use_distance: A boolean indicating whether distance information is used in the catalog.
        """
        self.cols = Columns(ra=ra, dec=dec)
        self.use_distance = use_distance
        self.matches = {}

    def add_match(self, catalog_name, from_to=0, mask=None, idx=None, sep2d=None, sep3d=None):
        """
        Adds information about a match with another catalog.

        Args:
            catalog_name: The name of the other catalog.
            from_to: The direction of the match (0 for self to other catalog, 1 for other catalog to self).
            mask: The mask indicating the matched objects.
            idx: The indices of the matched objects.
            sep2d: The 2D separation between the matched objects.
            sep3d: The 3D separation between the matched objects.
        """
        self.matches[catalog_name] = Match(catalog_name=catalog_name,
                                           from_to=from_to,
                                           mask=mask,
                                           idx=idx,
                                           sep2d=sep2d,
                                           sep3d=sep3d)

    def get_mask_to_third_catalog(self, cat1, cat2):
        """
        Returns the mask to select objects in the third catalog that are matched with objects in the first catalog.

        Args:
            cat1: The name of the first catalog.
            cat2: The name of the second catalog.

        Returns:
            The mask to select objects in the third catalog.
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
        """
        Loads the catalog data.
        """
        pass

    def filter_box(self, boxes):
        """
        Removes any object outside RA/DEC boxes.

        Args:
            boxes: The RA/DEC boxes to filter the catalog data.
        """
        df = filter_dataframe_to_box(self, boxes)
        self.df = df

    def as_SkyCoord(self, mask=None, distance=None, distance_unit=None):
        """
        Converts the catalog data to a SkyCoord object.

        Args:
            mask: The mask indicating the objects to include in the SkyCoord object.
            distance: The distance values to include in the SkyCoord object.
            distance_unit: The unit of the distance values.

        Returns:
            A SkyCoord object representing the catalog data.
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
        """
        Returns the value of a property of the catalog.

        Args:
            prop: The name of the property.
            default: The default value to return if the property is not found.

        Returns:
            The value of the property.
        """
        return getattr(self, prop, default)

    def prop_to_unit(self, prop, unit, mask=None, with_unit=False):
        """
        Converts a property of the catalog to a specified unit.

        Args:
            prop: The name of the property.
            unit: The unit to convert the property to.
            mask: The mask indicating the objects to include in the conversion.
            with_unit: A boolean indicating whether to include the unit in the result.

        Returns:
            The converted property value.
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
            return converted_prop.value
