import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

class Ellipse:
    def __init__(self, x, y, a, b, pa):
        """Initialise an ellipse

        Parameters
        ----------
        x:float
            x coordinate
        y:float
            y coordinate
        a:float
            semi-major axis
        b:float
            semi-minor axis
        pa:float
            Position angle
        """
        self.x = x
        self.y = y
        self.a = a  # semi-major axis
        self.b = b  # semi-minor axis
        self.pa = pa

    def is_point_inside(self, x, y):
        """Check if point (x,y) is inside ellipse

        Parameters
        ----------
            x:float
                x coordinate
            y:float
                y coordinate

        Returns
        -------
            True if (x,y) in ellipse else False
        """
        cos_sin = (np.cos(self.pa) * (x-self.x) +
                   np.sin(self.pa) * (y-self.y))
        sin_cos = (np.sin(self.pa) * (x-self.x) +
                   np.cos(self.pa) * (y-self.y))

        return cos_sin**2/self.a**2 + sin_cos**2/self.b**2 <= 1

    # TODO: code a way to check for overlapping ellipses
    #  see:
    #   - https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect
    #   - https://math.stackexchange.com/questions/1447730/drawing-ellipse-from-eigenvalue-eigenvector
    # @property
    # def E(self):
    #     X = np.array(
    #         [1, 1],
    #         [self.x, self.y],
    #         [1, 1]
    #     )
    #     Q = np.array([
    #         [1/self.a**2, 0],
    #         [0, 1/self.b**2]
    #     ])
    #     rot = np.array([
    #         [np.cos(self.pa), np.cos(self.pa)],
    #         [-np.sin(self.pa), np.cos(self.pa)]
    #     ])
    #     X_prime = rot * X
    #
    #     return X_prime * Q * X_prime.T
    #
    # @property
    # def sigma(self):
    #     return np.linalg.inv(self.E)

def crossmatch(obj1, obj2,
               sep=6, unit_sep=u.arcsec, ellipse_filter=False,
               mask1=None, mask2=None):
    """ Cross-matching Catalog 1 (source) to Catalog 2 (target)

    In addition to the two Catalog objects to be cross-matched,
    the function accepts a separation threshold (default: 6 arcsec),
    the separation unit (default: arcmin), whether to apply ellipse filter or not, and
    provide an index mask for obj1.

    The ellipse filter currently requires information about semi-major and semi-minor axes.
    TODO: make the ellipse filter usage less confusing... perhaps add options to choose the source by param or move this to a different function

    Results are stored in the match dictionary of both obj1 and obj2. E.g. obj1.matches[obj2.name] or obj2.matches[obj1.name]

    Parameters
    ----------
    obj1:Catalog
        Source catalog of sources
    obj2:Catalog
        Target catalog to be matched against source catalog
    sep:float
        Separation threshold for match association (default: 0.5/60)
    unit_sep:astropy.units
        Unit of separation value in astropy unit format (default: astropy.units.arcmin)
    ellipse_filter:bool
        Apply ellipse filter: yes or no (default: False)
    mask1:numpy.array
        Index mask for obj1 collection (default: None)
    mask2:numpy.array
        Index mask for obj2 collection (default: None)

    """
    source_coords = obj1.as_SkyCoord(mask1)
    _mask = np.zeros(len(source_coords), dtype=bool)

    unit = obj2.cols.ra.unit

    _ra = obj2.df[obj2.cols.ra.label].values if mask2 is None else obj2.df.iloc[mask2][obj2.cols.ra.label].values
    _dec = obj2.df[obj2.cols.dec.label].values if mask2 is None else obj2.df.iloc[mask2][obj2.cols.dec.label].values
    try:
        _ra_err = obj2.prop_to_unit('ra_err', unit, mask2)
        _dec_err = obj2.prop_to_unit('dec_err', unit, mask2)
    except (KeyError, AttributeError) as e:
        # TODO: would be best to avoid loop if no err provided
        _ra_err = 0.
        _dec_err = 0.

    # Astropy based first cross-match on coordinates and position error
    # TODO: check for parallelization
    for ra_err, dec_err in ((2,1), (0,1), (1,1),
                            (2,0),        (1,0),
                            (2,2), (0,2), (1,2), (0,0)) : # finishing with 0,0 to save d2d of central position.
        '''
        ra_err, dec_err description:

        x = ra
        y = dec

        Logic used with index:
        0:0
        1:1
        2:-1

        Coordinate space:
        (-1,1)  (0,1)  (1,1)
        (-1,0)  (0,0)  (1,0)
        (-1,-1) (0,-1) (1,-1)

        Constraint space:
        (2,1) (0,1) (1,1)
        (2,0) (0,0) (1,0)
        (2,2) (0,2) (1,2)
        
        'Logical or' motivation (|=):
        True | True == True
        True | False == True
        False | True == True
        False | False == False
        '''
        # print (_ra_err, _dec_err)
        if ra_err == 1:
            ra = (_ra + _ra_err)
        elif ra_err == 2:
            ra = (_ra - _ra_err)
        else:
            ra = _ra

        if dec_err == 1:
            dec = (_dec + _dec_err)
        elif dec_err == 2:
            dec = (_dec - _dec_err)
        else:
            dec = _dec

        # Verify we don't go beyond what's acceptable
        ra[np.where(ra > 360)[0]] = ra[np.where(ra > 360)[0]] - 360
        ra[np.where(ra < 0)[0]] = 360 + ra[np.where(ra < 0)[0]]
        ra = ra * unit
        
        dec[np.where(dec > 90)[0]] = -(dec[np.where(dec > 90)[0]] - 90)
        dec[np.where(dec < -90)[0]] = -(dec[np.where(dec < -90)[0]] + 90)
        dec = dec * unit
        try:
            target_coords = SkyCoord(ra=ra, dec=dec)
        except ValueError:
            print (dec.min(), dec.max())
        idx, sep2d, _ = source_coords.match_to_catalog_sky(target_coords)

        _mask |= (sep2d <= sep * unit_sep)

    if ellipse_filter:
        sc_obj1 = obj1.as_SkyCoord(mask=mask1)

        # Post quick astropy cross-match (circular matching) to now check if source is within a defined ellipse
        for i, obj1_id in enumerate(np.where(_mask)[0]):
            pa = obj1.df.iloc[obj1_id][obj1.cols.pa.label]
            gal = Ellipse(
                x = sc_obj1[obj1_id].ra.deg,
                y = sc_obj1[obj1_id].dec.deg,
                a = obj1.semi_major(mask=obj1_id, with_unit=False, to_unit=u.deg),
                b = obj1.semi_minor(mask=obj1_id, with_unit=False, to_unit=u.deg),
                pa = pa if not pd.isna(pa) else 0
            )
            obj2_mask = idx[obj1_id]
            sc_obj2 = obj2.as_SkyCoord(mask=obj2_mask)

            x, y = sc_obj2.ra.deg, \
                   sc_obj2.dec.deg,
            try:
                e_x, e_y = obj2.prop_to_unit('ra_err', unit, mask=obj2_mask), \
                           obj2.prop_to_unit('dec_err', unit, mask=obj2_mask)
            except AttributeError:
                e_x, e_y = obj2.semi_major(mask=obj2_mask, with_unit=False, to_unit=u.deg), \
                           obj2.semi_minor(mask=obj2_mask, with_unit=False, to_unit=u.deg)

            # Update mask  -- this line could be parallelized (vectorized?)
            b = gal.is_point_inside(x, y) | \
                    gal.is_point_inside(x-e_x, y-e_y) | gal.is_point_inside(x-e_x, y+e_y) | \
                    gal.is_point_inside(x+e_x, y-e_y) | gal.is_point_inside(x+e_x, y+e_y) | \
                    gal.is_point_inside(x, y-e_y) | gal.is_point_inside(x, y+e_y) | \
                    gal.is_point_inside(x-e_x, y) | gal.is_point_inside(x+e_x, y)

            _mask[obj1_id] = b

    obj1.add_match(catalog_name=obj2.name, from_to=0, mask=_mask)
    obj2.add_match(catalog_name=obj1.name, from_to=1, idx=idx, mask=_mask, sep2d=sep2d)

