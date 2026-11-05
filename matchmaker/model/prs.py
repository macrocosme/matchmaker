"""Reference persistent radio sources (PRS): the data the PRS-similarity search is centred on.

Single source of truth for these objects' measured properties. Every number carries its
reference; values not yet extracted from the literature are None, and callers must skip them
rather than guess. Plot styling for the same objects lives in matchmaker/dataset/sources.py,
which takes its shared values from here.

Conventions: luminosities stored at 144 MHz (log10 W/Hz), with the flux, frequency and alpha
used to get there; emission lines corrected for Galactic extinction, in 1e-16 erg/s/cm^2.
"""
import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM

COSMO = FlatLambdaCDM(H0=70, Om0=0.3)   # same as matchmaker.measures.DEFAULT_COSMO (Gurkan+2018 convention)


def log_l144(flux_uJy, freq_MHz, z, alpha):
    """log10 L_nu at 144 MHz in W/Hz, scaling flux from freq_MHz with S ~ nu^alpha and k-correcting."""
    s144 = flux_uJy * (144.0 / freq_MHz) ** alpha
    d = COSMO.luminosity_distance(z).to(u.m)
    s = (s144 * u.uJy).to(u.W / u.m**2 / u.Hz)
    return float(np.log10((4 * np.pi * d**2 * s / (1 + z) ** (1 + alpha)).to(u.W / u.Hz).value))


# Galactic extinction toward FRB 20121102A used by Tendulkar+2017: A_V = 2.42 (Schlegel+1998).
_AV_121102 = 2.42


def _unabsorb(flux, a_ratio, av=_AV_121102):
    return flux * 10 ** (0.4 * a_ratio * av)


# Tendulkar+2017 (ApJL 834, L7), Table 1: observed fluxes in 1e-16 erg/s/cm^2, and A_lambda/A_V.
# Upper limits are 3-sigma. Corrected here for Galactic extinction, although the BPT ratios move
# by <0.03 dex because each ratio pairs nearby wavelengths.
_T17 = {
    #            flux,  err,   A_l/A_V, upper_limit
    'h_beta':   (0.118, 0.011, 0.941, False),
    'oiii_5007': (0.575, 0.011, 0.911, False),
    'oi_6300':  (0.009, None,  0.670, True),
    'h_alpha':  (0.652, 0.009, 0.622, False),
    'nii_6584': (0.030, None,  0.619, True),
    'sii_6717': (0.040, 0.006, 0.596, False),
    'sii_6731': (0.024, 0.006, 0.593, False),
}
LINES_121102 = {
    k: dict(flux=_unabsorb(f, a), err=None if e is None else _unabsorb(e, a), upper_limit=ul)
    for k, (f, e, a, ul) in _T17.items()
}

FRB20121102A = dict(
    name='FRB 20121102A',
    z=0.19273, z_ref='Tendulkar+2017',
    # PRS: 203.5 +/- 33.6 uJy in uGMRT band 3 (200-300 MHz, Resmi+2021), flat alpha = -0.07 +/- 0.03.
    # Taking the band centre 250 MHz. Scaling sources.py's 2e29 erg/s/Hz at 1.4 GHz instead gives 22.37.
    log_l144=log_l144(203.5, 250.0, 0.19273, -0.07),
    log_l144_ref='Resmi+2021 250 MHz, alpha=-0.07',
    sfr=0.4, sfr_err=0.17, sfr_ref='Tendulkar+2017 (H-alpha); err from sources.py',
    mstar=1.3e8, mstar_ref='Bassa+2017 (Tendulkar+2017 give 4-7e7)',
    alpha=-0.07, alpha_err=0.03, alpha_ref='Resmi+2021',
    lines=LINES_121102, lines_ref='Tendulkar+2017 Table 1, Galactic-extinction corrected',
    wise=None,   # host is ~5.5 mag below the AllWISE W1 limit: cannot be placed in the WISE plane
)

FRB20190520B = dict(
    name='FRB 20190520B',
    z=0.241, z_ref='Niu+2022',
    # 202 +/- 8 uJy at 3 GHz (Niu+2022), alpha = -0.41 (Zhao & Wang 2021), as in sources.py.
    log_l144=log_l144(202.0, 3000.0, 0.241, -0.41),
    log_l144_ref='Niu+2022 3 GHz, alpha=-0.41',
    sfr=0.41, sfr_err=None, sfr_ref='Niu+2022',
    mstar=6e8, mstar_ref='Niu+2022',
    alpha=-0.41, alpha_err=0.04, alpha_ref='Zhao & Wang 2021',
    lines=None,  # TODO: host emission-line fluxes from Niu+2022 (Methods / Extended Data). Not yet extracted.
    lines_ref=None,
    wise=None,
)

ANCHORS = [FRB20121102A, FRB20190520B]
