import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord


class Ellipse:
    """Ellipse on the sky.

    Parameters
    ----------
    x, y:float
        Center RA, DEC in degrees
    a, b:float
        Semi-major and semi-minor axes in degrees (true angular size)
    pa:float
        Position angle of the major axis in degrees, measured East of North
        (standard astronomical convention)
    """
    def __init__(self, x, y, a, b, pa):
        self.x = x
        self.y = y
        self.a = a  # semi-major axis
        self.b = b  # semi-minor axis
        self.pa = pa

    def is_point_inside(self, x, y):
        """Check whether (ra, dec) in degrees falls inside the ellipse.

        Works in the local tangent plane: RA offsets are scaled by cos(dec)
        so distances are true angular distances.
        """
        pa_rad = np.radians(self.pa)
        # Tangent-plane offsets in degrees (RA wrap-safe)
        dx = (((x - self.x) + 180.) % 360. - 180.) * np.cos(np.radians(self.y))
        dy = y - self.y
        # Rotate into the ellipse frame: major axis at PA East of North
        along_major = dx * np.sin(pa_rad) + dy * np.cos(pa_rad)
        along_minor = dx * np.cos(pa_rad) - dy * np.sin(pa_rad)

        return along_major**2 / self.a**2 + along_minor**2 / self.b**2 <= 1

    # TODO: code a way to check for overlapping ellipses
    #  see:
    #   - https://math.stackexchange.com/questions/1114879/detect-if-two-ellipses-intersect
    #   - https://math.stackexchange.com/questions/1447730/drawing-ellipse-from-eigenvalue-eigenvector


class _LegacyEllipse:
    """Pre-2026 Ellipse, kept verbatim for the crossmatch legacy path.

    Known issues (why it was replaced): pa is used as radians although
    catalogs store degrees; the rotation matrix is not a true rotation
    (missing minus sign); RA offsets are not scaled by cos(dec); pa is
    measured from the RA axis instead of East of North.
    """
    def __init__(self, x, y, a, b, pa):
        self.x = x
        self.y = y
        self.a = a
        self.b = b
        self.pa = pa

    def is_point_inside(self, x, y):
        cos_sin = (np.cos(self.pa) * (x-self.x) +
                   np.sin(self.pa) * (y-self.y))
        sin_cos = (np.sin(self.pa) * (x-self.x) +
                   np.cos(self.pa) * (y-self.y))

        return cos_sin**2/self.a**2 + sin_cos**2/self.b**2 <= 1


def _positional_error_arcsec(obj, mask, n):
    """Radial 1-sigma positional uncertainty per row in arcsec.

    Zeros when the catalog declares no ra_err/dec_err columns; NaN errors
    are treated as zero.
    """
    try:
        ra_err = obj.prop_to_unit('ra_err', u.arcsec, mask)
        dec_err = obj.prop_to_unit('dec_err', u.arcsec, mask)
    except (KeyError, AttributeError):
        return np.zeros(n)
    return np.hypot(np.nan_to_num(ra_err), np.nan_to_num(dec_err))


def _ellipse_inflation_deg(obj2, j, sigma_target_arcsec, n_sigma):
    """How much to grow the source ellipse axes (deg) for target j.

    Uses the target's positional uncertainty when available; otherwise falls
    back to the target's own angular extent (mirrors the legacy behavior of
    testing the target's semi-axes against the ellipse).
    """
    if sigma_target_arcsec[j] > 0:
        return (n_sigma * sigma_target_arcsec[j] * u.arcsec).to(u.deg).value
    try:
        return max(obj2.semi_major(mask=j, with_unit=False, to_unit=u.deg),
                   obj2.semi_minor(mask=j, with_unit=False, to_unit=u.deg))
    except AttributeError:
        return 0.


def crossmatch(obj1, obj2,
               sep=0.5/60, unit_sep=u.arcmin, ellipse_filter=False,
               mask1=None, mask2=None, n_sigma=1.0, legacy=False, target_sep=None):
    """ Cross-matching Catalog 1 (source) to Catalog 2 (target)

    A single nearest-neighbour sky match; a source is matched when

        sep2d <= sep * unit_sep + n_sigma * sigma       (or target_sep[idx] * unit_sep + ...)

    where sigma is the combined radial positional uncertainty of the pair
    (source and matched target ra_err/dec_err added in quadrature; zero for
    catalogs that declare no error columns). This grows the acceptance
    radius smoothly with positional uncertainty instead of the historical
    9-point shift grid, whose union-of-circles acceptance region had holes
    (a closer target could be rejected while a farther one was accepted).

    With ellipse_filter=True, matched pairs are additionally required to
    have the target inside the source's ellipse (semi_major/semi_minor/pa),
    with axes inflated by the target's positional uncertainty (or the
    target's own extent when it has no error columns).

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
    n_sigma:float
        How many combined positional-error sigmas to add to the acceptance
        radius (default: 1.0; set 0 for a pure fixed-radius match)
    legacy:bool
        Reproduce the historical 9-point shift-grid matching and legacy
        ellipse filter verbatim, quirks included (default: False). Intended
        for diffing results against previous runs.
    target_sep:array, optional
        Per-target acceptance radius in unit_sep, one value per target row
        being matched (after mask2). Replaces sep: each source is tested
        against the radius of its nearest target, e.g. LoTSS DR3's
        per-source resolution (6" or 9"). Not available with legacy=True.

    Results are stored on both catalogs via Catalog.add_match():
    obj1.matches[obj2.name].mask (bool over source rows) and
    obj2.matches[obj1.name].idx/mask/sep2d.
    """
    if legacy:
        if target_sep is not None:
            raise ValueError('target_sep is not supported with legacy=True')
        return _crossmatch_legacy(obj1, obj2, sep=sep, unit_sep=unit_sep,
                                  ellipse_filter=ellipse_filter,
                                  mask1=mask1, mask2=mask2)

    source_coords = obj1.as_SkyCoord(mask1)

    unit = obj2.cols.ra.unit
    _ra = obj2.df[obj2.cols.ra.label].values if mask2 is None else obj2.df.iloc[mask2][obj2.cols.ra.label].values
    _dec = obj2.df[obj2.cols.dec.label].values if mask2 is None else obj2.df.iloc[mask2][obj2.cols.dec.label].values
    target_coords = SkyCoord(ra=_ra * unit, dec=_dec * unit)

    idx, sep2d, _ = source_coords.match_to_catalog_sky(target_coords)

    sigma_source = _positional_error_arcsec(obj1, mask1, len(source_coords))
    sigma_target = _positional_error_arcsec(obj2, mask2, len(target_coords))
    sigma = np.hypot(sigma_source, sigma_target[idx]) * u.arcsec

    if target_sep is None:
        radius = sep * unit_sep
    else:
        target_sep = np.asarray(target_sep, dtype=float)
        if target_sep.shape != (len(target_coords),):
            raise ValueError(f'target_sep has shape {target_sep.shape}; expected ({len(target_coords)},), '
                             'one radius per target row being matched')
        radius = target_sep[idx] * unit_sep

    _mask = sep2d <= (radius + n_sigma * sigma)

    if ellipse_filter:
        pa_unit = obj1.cols.pa.unit if obj1.cols.pa.unit is not None else u.deg
        for obj1_id in np.where(_mask)[0]:
            j = idx[obj1_id]

            pa = obj1.df.iloc[obj1_id][obj1.cols.pa.label]
            pa = 0. if pd.isna(pa) else (pa * pa_unit).to(u.deg).value

            inflation = _ellipse_inflation_deg(obj2, j, sigma_target, n_sigma)
            gal = Ellipse(
                x = source_coords[obj1_id].ra.deg,
                y = source_coords[obj1_id].dec.deg,
                a = obj1.semi_major(mask=obj1_id, with_unit=False, to_unit=u.deg) + inflation,
                b = obj1.semi_minor(mask=obj1_id, with_unit=False, to_unit=u.deg) + inflation,
                pa = pa
            )
            _mask[obj1_id] = gal.is_point_inside(target_coords[j].ra.deg,
                                                 target_coords[j].dec.deg)

    obj1.add_match(catalog_name=obj2.name, from_to=0, mask=_mask)
    obj2.add_match(catalog_name=obj1.name, from_to=1, idx=idx, mask=_mask, sep2d=sep2d)


def _crossmatch_legacy(obj1, obj2,
                       sep=0.5/60, unit_sep=u.arcmin, ellipse_filter=False,
                       mask1=None, mask2=None):
    """Historical shift-grid cross-match, kept verbatim for comparison runs.

    Known issues (why it was replaced): the union-of-9-circles acceptance
    region has holes; idx/sep2d come from the unshifted pass and can be
    inconsistent with mask; RA error shifts are not scaled by cos(dec);
    declination overflow past the poles maps to the wrong hemisphere; the
    ellipse filter uses _LegacyEllipse (see its docstring).
    """
    source_coords = obj1.as_SkyCoord(mask1)
    _mask = np.zeros(len(source_coords), dtype=bool)

    unit = obj2.cols.ra.unit

    _ra = obj2.df[obj2.cols.ra.label].values if mask2 is None else obj2.df.iloc[mask2][obj2.cols.ra.label].values
    _dec = obj2.df[obj2.cols.dec.label].values if mask2 is None else obj2.df.iloc[mask2][obj2.cols.dec.label].values
    try:
        _ra_err = obj2.prop_to_unit('ra_err', unit, mask2)
        _dec_err = obj2.prop_to_unit('dec_err', unit, mask2)
    except (KeyError, AttributeError):
        _ra_err = 0.
        _dec_err = 0.

    # Astropy based first cross-match on coordinates and position error
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
        # Copies: the wraparound below writes in place, and _ra/_dec view the catalog's data
        if ra_err == 1:
            ra = (_ra + _ra_err)
        elif ra_err == 2:
            ra = (_ra - _ra_err)
        else:
            ra = _ra.copy()

        if dec_err == 1:
            dec = (_dec + _dec_err)
        elif dec_err == 2:
            dec = (_dec - _dec_err)
        else:
            dec = _dec.copy()

        # Verify we don't go beyond what's acceptable
        ra[np.where(ra > 360)[0]] = ra[np.where(ra > 360)[0]] - 360
        ra[np.where(ra < 0)[0]] = 360 + ra[np.where(ra < 0)[0]]
        ra = ra * unit

        dec[np.where(dec > 90)[0]] = -(dec[np.where(dec > 90)[0]] - 90)
        dec[np.where(dec < -90)[0]] = -(dec[np.where(dec < -90)[0]] + 90)
        dec = dec * unit
        try:
            target_coords = SkyCoord(ra=ra, dec=dec)
        except ValueError as e:
            # Re-raise with coordinate ranges: proceeding without target_coords would fail anyway
            raise ValueError(
                f'Could not build target SkyCoord '
                f'(ra range: {ra.min()}..{ra.max()}, dec range: {dec.min()}..{dec.max()})'
            ) from e
        idx, sep2d, _ = source_coords.match_to_catalog_sky(target_coords)

        _mask |= (sep2d <= sep * unit_sep)

    if ellipse_filter:
        sc_obj1 = obj1.as_SkyCoord(mask=mask1)

        # Post quick astropy cross-match (circular matching) to now check if source is within a defined ellipse
        for _i, obj1_id in enumerate(np.where(_mask)[0]):
            pa = obj1.df.iloc[obj1_id][obj1.cols.pa.label]
            gal = _LegacyEllipse(
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
