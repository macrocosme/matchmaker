"""Radio spectral indices, S ~ nu^alpha.

Ported from the paper-time module (archive/reproduction_package/matchmaker/matchmaker/spectral_index.py),
keeping its fitting convention (power law in linear flux, absolute errors) and its two indices:

  alpha       fitted to every point
  alpha_high  fitted without the lowest-frequency point; alpha - alpha_high > 0 means the spectrum
              flattens toward low frequency (a turnover / absorption signature)

Not ported: add_flux_from_other_survey and get_arrays, which were tied to a dict of Catalogs
and carried bugs (a broken scalar branch; total flux copied into the peak-flux fields).
Undefined uncertainties are returned as NaN, never 0.
"""
import numpy as np

from .fitting import fit_power_law


def fit_spectral_index(freqs, fluxes, flux_errs):
    """alpha and its 1-sigma error from points with positive finite flux and error.

    freqs in any unit (only ratios matter); returns (nan, nan) with fewer than 2 usable points
    or if the fit fails. With exactly 2 points the fit has zero degrees of freedom; the error
    still propagates the flux errors because absolute_sigma=True.
    """
    f, s, e = (np.asarray(x, dtype=float) for x in (freqs, fluxes, flux_errs))
    ok = np.isfinite(f) & np.isfinite(s) & np.isfinite(e) & (f > 0) & (s > 0) & (e > 0)
    f, s, e = f[ok], s[ok], e[ok]
    if len(f) < 2:
        return np.nan, np.nan
    f = f / np.median(f)   # rescale so the amplitude starts near the data; alpha is unaffected
    try:
        _, k, _, dk = fit_power_law(f, s, e)
    except (RuntimeError, ValueError):
        return np.nan, np.nan
    return float(k), (float(dk) if np.isfinite(dk) else np.nan)


def spectral_indices(freqs, fluxes, flux_errs, min_high_points=2):
    """alpha over all points and alpha_high without the lowest frequency.

    Returns dict(alpha, alpha_err, alpha_high, alpha_high_err, n_points). alpha_high needs at
    least min_high_points usable points above the lowest frequency.
    """
    f, s, e = (np.asarray(x, dtype=float) for x in (freqs, fluxes, flux_errs))
    usable = np.isfinite(f) & np.isfinite(s) & np.isfinite(e) & (s > 0) & (e > 0)
    order = np.argsort(f)
    f, s, e, usable = f[order], s[order], e[order], usable[order]
    alpha, alpha_err = fit_spectral_index(f, s, e)
    high = usable.copy()
    if usable.any():
        high[np.argmax(usable)] = False          # drop the lowest usable frequency
    if high.sum() >= min_high_points:
        a_h, a_h_err = fit_spectral_index(f[high], s[high], e[high])
    else:
        a_h, a_h_err = np.nan, np.nan
    return dict(alpha=alpha, alpha_err=alpha_err, alpha_high=a_h, alpha_high_err=a_h_err,
                n_points=int(usable.sum()))
