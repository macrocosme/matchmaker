"""Survey footprints built from source positions.

A footprint is the set of fine HEALPix cells that contain at least one catalogue source. This
replaces RA/Dec bounding boxes (which overcount any footprint that is not a rectangle) and the
coarse mosaic tiling some catalogues ship (LoTSS DR3's HEALPIX column is nside 16, 13 deg^2).

At LoTSS DR3 density (~0.2 sources per arcmin^2) an nside-256 cell (0.052 deg^2) holds ~30
sources, so a truly empty cell inside the survey is vanishingly rare (P ~ e^-30). Cells on the
survey edge count in full, so the area is overestimated by about one cell width times the
perimeter. For sparser catalogues, lower nside so a cell holds several sources.
"""
import astropy.units as u
import numpy as np


def _healpix(nside):
    # Imported here so DR2-only use of the package does not require astropy-healpix
    from astropy_healpix import HEALPix
    return HEALPix(nside=nside, order='nested')


class Footprint:
    """Sky coverage as a set of occupied HEALPix cells, optionally with one value per cell.

    Parameters
    ----------
    cells : array of int
        Occupied cell indices (nested ordering).
    nside : int
    values : array, optional
        One value per cell (same order as cells), e.g. the survey resolution there.
    """

    def __init__(self, cells, nside=256, values=None):
        order = np.argsort(cells)
        self.cells = np.asarray(cells, dtype=np.int64)[order]
        self.nside = nside
        self.values = None if values is None else np.asarray(values)[order]
        self._grid = _healpix(nside)

    @classmethod
    def from_positions(cls, ra, dec, nside=256, values=None):
        """Footprint of the cells holding the given positions (degrees).

        values, if given (one per position), are reduced to the most common value in each cell.
        """
        grid = _healpix(nside)
        cell = np.asarray(grid.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg), dtype=np.int64)
        cells, inverse = np.unique(cell, return_inverse=True)
        if values is None:
            return cls(cells, nside)
        values = np.asarray(values)
        kinds = np.unique(values)
        counts = np.zeros((len(cells), len(kinds)), dtype=np.int64)
        np.add.at(counts, (inverse, np.searchsorted(kinds, values)), 1)
        return cls(cells, nside, kinds[np.argmax(counts, axis=1)])

    def _lookup(self, ra, dec):
        cell = np.asarray(self._grid.lonlat_to_healpix(np.asarray(ra) * u.deg, np.asarray(dec) * u.deg),
                          dtype=np.int64)
        pos = np.clip(np.searchsorted(self.cells, cell), 0, len(self.cells) - 1)
        return pos, self.cells[pos] == cell

    def contains(self, ra, dec):
        """True for positions (degrees) inside the footprint."""
        return self._lookup(ra, dec)[1]

    def value_at(self, ra, dec, default=np.nan):
        """Per-cell value at each position, or default outside the footprint."""
        if self.values is None:
            raise ValueError('this footprint was built without values')
        pos, inside = self._lookup(ra, dec)
        out = np.full(len(pos), default, dtype=np.result_type(self.values, np.asarray(default)))
        out[inside] = self.values[pos[inside]]
        return out

    @property
    def cell_area(self):
        return self._grid.pixel_area.to(u.deg**2)

    @property
    def area(self):
        """Total area, as an astropy Quantity in deg^2."""
        return len(self.cells) * self.cell_area

    def area_where(self, value):
        """Area of the cells whose value equals `value`."""
        if self.values is None:
            raise ValueError('this footprint was built without values')
        same = np.isclose(self.values, value) if np.issubdtype(self.values.dtype, np.number) else self.values == value
        return int(np.sum(same)) * self.cell_area
