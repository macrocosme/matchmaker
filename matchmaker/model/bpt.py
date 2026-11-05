"""Emission-line (BPT) diagnostics.

Line ratios with propagated uncertainties (`log_ratio`), and the standard demarcation lines.
"""
from math import log

import numpy as np

LN10 = log(10.0)


def log_ratio(a, a_err, b, b_err):
    """log10(a/b) and its 1-sigma uncertainty (1/ln10) sqrt((a_err/a)^2 + (b_err/b)^2)."""
    a, a_err, b, b_err = (np.asarray(x, dtype=float) for x in (a, a_err, b, b_err))
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.log10(a / b)
        err = np.sqrt((a_err / a) ** 2 + (b_err / b) ** 2) / LN10
    return ratio, err


def kauffmann03_nii(log_nii_ha):
    """Kauffmann+2003 empirical star-forming limit on the [NII]/Ha diagram (valid for x < 0.05)."""
    return 0.61 / (np.asarray(log_nii_ha, dtype=float) - 0.05) + 1.30


def kewley01_nii(log_nii_ha):
    """Kewley+2001 maximum-starburst line on the [NII]/Ha diagram (valid for x < 0.47)."""
    return 0.61 / (np.asarray(log_nii_ha, dtype=float) - 0.47) + 1.19


def kewley01_sii(log_sii_ha):
    """Kewley+2001 maximum-starburst line on the [SII]/Ha diagram (valid for x < 0.32)."""
    return 0.72 / (np.asarray(log_sii_ha, dtype=float) - 0.32) + 1.30


def kewley06_sii_seyfert_liner(log_sii_ha):
    """Kewley+2006 Seyfert/LINER division on the [SII]/Ha diagram."""
    return 1.89 * np.asarray(log_sii_ha, dtype=float) + 0.76


#: Where the Seyfert/LINER line meets the maximum-starburst line on the [SII] diagram.
#: 1.89x + 0.76 = 0.72/(x - 0.32) + 1.30 has roots at -0.3146 and 0.9203; the first is the
#: one on the branch of the starburst line that bounds the star-forming region.
SII_SEYFERT_LINER_START = -0.3145942791


def kewley06_sii_valid(log_sii_ha):
    """True where the Seyfert/LINER division applies: above the maximum-starburst line.

    The division separates two kinds of ionisation among galaxies that are *already* above
    Kewley+2001. Drawing it below that line puts a boundary through the star-forming region,
    where it means nothing.
    """
    return np.asarray(log_sii_ha, dtype=float) >= SII_SEYFERT_LINER_START


def kewley01_oi(log_oi_ha):
    """Kewley+2001 maximum-starburst line on the [OI]/Ha diagram (valid for x < -0.59)."""
    return 0.73 / (np.asarray(log_oi_ha, dtype=float) + 0.59) + 1.33


def classify_nii(log_nii_ha, log_oiii_hb):
    """'star-forming', 'composite' or 'AGN' on the [NII] diagram; '' where either ratio is NaN."""
    x = np.asarray(log_nii_ha, dtype=float)
    y = np.asarray(log_oiii_hb, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        sf = (x < 0.05) & (y < kauffmann03_nii(x))
        comp = ~sf & (x < 0.47) & (y < kewley01_nii(x))
    out = np.where(sf, 'star-forming', np.where(comp, 'composite', 'AGN'))
    return np.where(np.isfinite(x) & np.isfinite(y), out, '')
