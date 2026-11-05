"""The 144 MHz radio luminosity - star formation rate relation, and its scatter.

`power_law` is the relation Vohl+2023 selected on: L144 = 10^22.06 psi^1.07 (Gurkan+2018),
with `sigma` the scatter used to define the +3 sigma outlier cut.

`broken_power_law` adds the stellar-mass term and the change of slope at low star formation
rate. It is kept here because the mass term matters for dwarfs, where the mass-independent
relation over-predicts, but the constants below came from a work-in-progress citation in the
original code and should be checked against the published paper before any result rests on
them; pass them explicitly if you have better values.
"""
import numpy as np

sigma = 1.45
sigma_err = 0.04


def power_law(psi, beta=1.07, l1=10 ** 22.06):
    """L144 in W/Hz from the star formation rate in Msun/yr. Gurkan+2018, mass-independent."""
    return l1 * np.asarray(psi, dtype=float) ** beta


def broken_power_law(psi, mstar, log_lc=22.02, gamma=0.44, beta_low=0.52, beta_high=1.01,
                     psi_break=0.01):
    """L144 in W/Hz with a stellar-mass term and a break in slope at `psi_break`.

        L = 10^log_lc * psi^beta * (mstar / 1e10)^gamma

    with beta = beta_low below the break and beta_high above it, the high branch carrying
    psi_break^(beta_low - beta_high) so the two meet at the break. Both arguments may be arrays.
    """
    psi = np.asarray(psi, dtype=float)
    mstar = np.asarray(mstar, dtype=float)
    mass_term = (mstar / 1e10) ** gamma
    with np.errstate(divide='ignore', invalid='ignore'):
        low = 10 ** log_lc * psi ** beta_low * mass_term
        high = 10 ** log_lc * psi ** beta_high * mass_term * psi_break ** (beta_low - beta_high)
    return np.where(psi <= psi_break, low, high)


def measured_mass_power_law(psi, mstar, log_l1, beta, gamma, m_ref=1e9):
    """log L = log_l1 + beta log(psi) + gamma log(mstar / m_ref), as a luminosity in W/Hz.

    Fitted by `analysis/followup/fit_mass_dependence.py`. Separating the mass term brings the
    star formation slope back to about 1, so a selection that ignores the mass term
    over-predicts the faintest hosts by an order of magnitude.
    """
    psi = np.asarray(psi, dtype=float)
    mstar = np.asarray(mstar, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return 10 ** log_l1 * psi ** beta * (mstar / m_ref) ** gamma


def measured_power_law(psi, log_l1, beta):
    """The same form as `power_law`, with a normalisation and slope measured rather than assumed.

    `analysis/followup/fit_radio_sfr.py` fits these from the stack; the DR3 dwarfs give a
    steeper slope than 1.07, which changes what counts as an outlier.
    """
    return 10 ** log_l1 * np.asarray(psi, dtype=float) ** beta
