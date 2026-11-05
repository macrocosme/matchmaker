"""Point-source (compactness) envelope in ln(S_I / S_P) versus signal-to-noise.

A source whose integrated flux exceeds its peak brightness by more than noise allows is
resolved. Shimwell et al. (2019, 2022) draw an envelope R_99.9(S/N) below which 99.9% of
unresolved sources fall, of the form

    R_99.9(x) = a + b / (1 + (x / c)^d),        x = S_I / sigma_S_I

and call a source compact when ln(S_I / S_P) <= R_99.9. The coefficients are release-specific
(they absorb the smearing and the noise behaviour of that imaging), so this module measures them
from a catalogue rather than hard-coding them. DR3 publishes none.

Measurement, per S/N bin:
  - unresolved sources scatter symmetrically about a centre c0, which is positive because time-
    and bandwidth-smearing lower the peak brightness. Resolved sources only add a positive tail,
    so the centre is the *mode*, not the median (in DR3 the median rises with S/N as bright
    extended sources come to dominate, and would give a nonsensically wide envelope);
  - the low tail is therefore free of resolved sources: mirroring the 0.1st percentile about the
    mode gives the value below which 99.9% of unresolved sources lie.
"""
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import norm

# Shimwell+2022 (LoTSS DR2). Kept for reference and for reproducing the published DR2 selection.
ENVELOPE_DR2 = (0.42, 1.08, 96.57, 2.49)


def envelope(snr, a, b, c, d):
    """R_99.9 at a given signal-to-noise."""
    return a + b / (1 + (np.asarray(snr, dtype=float) / c) ** d)


def is_compact(ln_si_sp, snr, params=ENVELOPE_DR2):
    """True where ln(S_I/S_P) lies below the envelope."""
    return np.asarray(ln_si_sp, dtype=float) <= envelope(snr, *params)


def _mode(x, lo=-1.0, hi=1.5, bins=250):
    """Peak of the ln(S_I/S_P) histogram: the centre of the unresolved population."""
    x = x[(x > lo) & (x < hi)]
    if len(x) < 10:
        return np.nan
    counts, edges = np.histogram(x, bins=bins)
    i = int(np.argmax(counts))
    return 0.5 * (edges[i] + edges[i + 1])


def measure_envelope(ln_si_sp, snr, target_percentile=99.9, tail_percentiles=(1.0, 10.0),
                     min_per_bin=500, n_bins=25, snr_min=5.0, method='tail'):
    """R_99.9 per S/N bin, from the low tail of ln(S_I/S_P), which resolved sources cannot reach.

    method='tail' (default): fit a Gaussian to two low percentiles and read off the matching upper
    percentile. The pair (1, 10) is used rather than the extreme tail: on simulated catalogues it
    is the most accurate of the pairs tried, and on DR3's 9" imaging it gives the smoothest curve.
    method='mode': mirror the low percentile about the histogram mode. Model-free and slightly
    better on clean simulations, but the mode stops tracking unresolved sources wherever extended
    ones dominate (DR3 9" above S/N ~ 50), so it is kept as a diagnostic rather than the default.

    Both assume the low percentiles are free of resolved sources, which fails if more than roughly
    85% of a bin is resolved; check fraction_compact in the diagnostics if that is a concern.

    Returns (snr_centre, r999, n_per_bin, mode), one entry per usable bin.
    """
    r = np.asarray(ln_si_sp, dtype=float)
    x = np.asarray(snr, dtype=float)
    ok = np.isfinite(r) & np.isfinite(x) & (x >= snr_min)
    r, x = r[ok], x[ok]
    z_lo, z_hi = norm.ppf(np.asarray(tail_percentiles) / 100)
    z_target = norm.ppf(target_percentile / 100)
    edges = np.geomspace(snr_min, max(np.percentile(x, 99.99), snr_min * 10), n_bins + 1)
    centre, r999, counts, modes = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        b = (x >= lo) & (x < hi)
        if b.sum() < min_per_bin:
            continue
        c0 = _mode(r[b])
        p_lo, p_hi = np.percentile(r[b], tail_percentiles)
        if method == 'tail':
            sigma = (p_hi - p_lo) / (z_hi - z_lo)
            value = (p_hi - z_hi * sigma) + z_target * sigma
        elif method == 'mode':
            if not np.isfinite(c0):
                continue
            value = c0 + (c0 - np.percentile(r[b], 100 - target_percentile))
        else:
            raise ValueError(f"method must be 'tail' or 'mode', not {method!r}")
        centre.append(np.sqrt(lo * hi))
        modes.append(c0)
        r999.append(value)
        counts.append(int(b.sum()))
    return (np.array(centre), np.array(r999), np.array(counts), np.array(modes))


def fit_envelope(ln_si_sp, snr, p0=ENVELOPE_DR2, **kwargs):
    """Fit a + b/(1+(x/c)^d) to the measured per-bin envelope, weighted by bin population.

    Returns (params, diagnostics) where diagnostics holds the binned measurements and the
    fraction of sources the fitted envelope calls compact.
    """
    centre, r999, counts, modes = measure_envelope(ln_si_sp, snr, **kwargs)
    if len(centre) < 5:
        raise ValueError(f'only {len(centre)} usable S/N bins; need at least 5')
    params, _ = curve_fit(envelope, centre, r999, p0=p0, sigma=1 / np.sqrt(counts),
                          absolute_sigma=False, maxfev=20000,
                          bounds=([-1.0, 0.0, 0.1, 0.1], [2.0, 20.0, 1e5, 10.0]))
    compact = is_compact(ln_si_sp, snr, params)
    return tuple(float(p) for p in params), {
        'snr_bin_centre': centre.tolist(), 'r999_measured': r999.tolist(),
        'n_per_bin': counts.tolist(), 'mode_per_bin': modes.tolist(),
        'r999_fitted': envelope(centre, *params).tolist(),
        'fraction_compact': float(np.nanmean(compact)),
    }
