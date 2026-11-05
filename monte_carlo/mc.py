"""Chance-association Monte Carlo: how many L-SFR outliers appear from unrelated pairs?

The original version estimated this by shifting every LoTSS position by a random offset, re-running the
cross-match and the +3 sigma outlier selection, and counting what survived. This rewrite keeps
that experiment and generalises it to DR3 (per-source matching radius, footprint), with two
changes:

- **The galaxies are shifted, not the radio sources.** Only the relative offset matters, and a
  galaxy catalogue is small: the radio positions can then be put in one search tree that is
  reused across realisations. On DR3 that is ~37,000 shifted positions against a fixed tree of
  5.6 million sources, which takes a second per realisation instead of rebuilding the tree.
- **Right-ascension offsets are divided by cos(dec)**, so the displacement is a true angle. The
  original version shifted RA in degrees directly, which at dec = 60 is half the intended
  offset.

The original version also deep-copied and pickled the whole dataset dictionary for every
realisation; this one keeps arrays and returns them.
"""
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord

from matchmaker.measures import luminosity_distance
from matchmaker.model.lofar import power_law
from matchmaker.model.lofar import sigma as LSFR_SIGMA

LOG_LSFR_SIGMA = np.log10(LSFR_SIGMA)


def random_offsets(n, shift_arcmin, dec_deg, rng):
    """Uniform offsets in a +/- shift_arcmin box, as true angles (RA scaled by 1/cos(dec))."""
    off = rng.uniform(-shift_arcmin, shift_arcmin, size=(2, n)) / 60.0
    return off[0] / np.cos(np.radians(dec_deg)), off[1]


def shift_positions(ra, dec, d_ra, d_dec):
    """Apply offsets in degrees, wrapping RA and reflecting at the poles."""
    new_ra = (np.asarray(ra) + d_ra) % 360.0
    new_dec = np.asarray(dec) + d_dec
    over = new_dec > 90
    new_dec[over] = 180 - new_dec[over]
    under = new_dec < -90
    new_dec[under] = -180 - new_dec[under]
    new_ra[over | under] = (new_ra[over | under] + 180) % 360.0
    return new_ra, new_dec


def run(galaxy_ra, galaxy_dec, galaxy_z, galaxy_sfr, radio_coords, radio_flux,
        radius_arcsec, radio_sigma_arcsec=None, n_sigma_position=1.0, n_sigma_lsfr=3.0,
        alpha=-0.8, n_realisations=1000, shift_arcmin=10.0, footprint=None, region=None,
        rng=None, progress=False, relation=None, scatter_dex=None, galaxy_mstar=None):
    """Counts of chance matches and chance L-SFR outliers, one entry per realisation.

    Parameters mirror the real selection, so pass exactly what the search used:
    radius_arcsec is a scalar (DR2's 6") or one value per radio source (DR3's resolution);
    radio_sigma_arcsec is the radio positional uncertainty per source, added as
    n_sigma_position * sigma, as matchmaker.crossmatch does. Luminosities use the galaxy's
    redshift with the unrelated source's flux, which is the point of the exercise.

    footprint (matchmaker.footprint.Footprint) drops shifted galaxies that land outside the
    survey; region, if given with a footprint, is a per-galaxy label whose counts are returned
    separately (e.g. the 6" and 9" parts of DR3).

    relation and scatter_dex replace the Gurkan+2018 luminosity - star formation relation and
    its 0.161 dex scatter, so the chance expectation can be computed against a relation measured
    from the data rather than the assumed one. Pass the same pair the selection used, or the
    excess is a comparison between two different definitions of outlier. With galaxy_mstar the
    relation is called as relation(sfr, mstar), for a relation that carries a stellar-mass term;
    the mass belongs to the galaxy, as the redshift does, while the flux belongs to the
    unrelated radio source, which is the whole point of the experiment.

    Returns a dict of arrays: n_galaxies, n_matches, n_outliers (plus n_matches_<label> and
    n_outliers_<label> per region).
    """
    rng = np.random.default_rng() if rng is None else rng
    ra = np.asarray(galaxy_ra, dtype=float)
    dec = np.asarray(galaxy_dec, dtype=float)
    z = np.asarray(galaxy_z, dtype=float)
    sfr = np.asarray(galaxy_sfr, dtype=float)
    mstar = None if galaxy_mstar is None else np.asarray(galaxy_mstar, dtype=float)
    flux = np.asarray(radio_flux, dtype=float)
    radius = np.asarray(radius_arcsec, dtype=float)
    sigma = None if radio_sigma_arcsec is None else np.asarray(radio_sigma_arcsec, dtype=float)
    predict = power_law if relation is None else relation
    log_sigma = LOG_LSFR_SIGMA if scatter_dex is None else float(scatter_dex)
    labels = [] if region is None else sorted(set(np.asarray(region)))
    out = {k: [] for k in ['n_galaxies', 'n_matches', 'n_outliers']
           + [f'{k}_{v}' for v in labels for k in ('n_matches', 'n_outliers')]}

    steps = range(n_realisations)
    if progress:
        from tqdm import tqdm
        steps = tqdm(steps)
    for _ in steps:
        d_ra, d_dec = random_offsets(len(ra), shift_arcmin, dec, rng)
        new_ra, new_dec = shift_positions(ra, dec, d_ra, d_dec)
        keep = np.ones(len(ra), dtype=bool) if footprint is None else footprint.contains(new_ra, new_dec)
        idx, sep, _ = SkyCoord(new_ra[keep] * u.deg, new_dec[keep] * u.deg).match_to_catalog_sky(radio_coords)

        accept = radius[idx] if radius.ndim else radius
        if sigma is not None:
            accept = accept + n_sigma_position * sigma[idx]
        matched = sep.arcsec <= accept

        z_m, sfr_m = z[keep][matched], sfr[keep][matched]
        with np.errstate(divide='ignore', invalid='ignore'):
            log_l = np.log10(luminosity_distance(z_m, flux[idx[matched]], alpha=alpha))
            expected = (predict(sfr_m) if mstar is None
                        else predict(sfr_m, mstar[keep][matched]))
            sigma_lsfr = (log_l - np.log10(expected)) / log_sigma
        outlier = sigma_lsfr >= n_sigma_lsfr

        out['n_galaxies'].append(int(keep.sum()))
        out['n_matches'].append(int(matched.sum()))
        out['n_outliers'].append(int(np.nansum(outlier)))
        if labels:
            reg = np.asarray(region)[keep][matched]
            for v in labels:
                here = reg == v
                out[f'n_matches_{v}'].append(int(here.sum()))
                out[f'n_outliers_{v}'].append(int(np.nansum(outlier & here)))
    return {k: np.array(v) for k, v in out.items()}


def from_catalogs(clu, lotss, radius_arcsec=None, **kwargs):
    """Run the Monte Carlo from loaded Catalog objects (either release).

    The radius defaults to each LoTSS source's resolution for DR3 and 6" for DR2, matching what
    the selection used. Pass the same cuts: the catalogues must already be constrained.
    """
    if radius_arcsec is None:
        radius_arcsec = lotss.resolution() if getattr(lotss, 'release', 'dr2') == 'dr3' else 6.0
    radio = SkyCoord(lotss.df[lotss.cols.ra.label].values * u.deg,
                     lotss.df[lotss.cols.dec.label].values * u.deg)
    sigma = np.hypot(lotss.df[lotss.cols.ra_err.label].values, lotss.df[lotss.cols.dec_err.label].values)
    if kwargs.get('galaxy_mstar') is None and getattr(clu.cols, 'mstar', None) is not None:
        label = clu.cols.mstar.label
        if label in clu.df:
            kwargs['galaxy_mstar'] = clu.df[label].values
    return run(clu.df[clu.cols.ra.label].values, clu.df[clu.cols.dec.label].values,
               clu.df[clu.cols.z.label].values, clu.df[clu.cols.sfr.label].values,
               radio, lotss.df[lotss.cols.measure.total_flux.label].values,
               radius_arcsec, radio_sigma_arcsec=sigma, **kwargs)


def p_value(counts, observed):
    """Fraction of realisations with at least `observed` chance outliers (0 means < 1/N)."""
    counts = np.asarray(counts)
    return float(np.mean(counts >= observed))
