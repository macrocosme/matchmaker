import astropy.units as u

sigma = 1.45
sigma_err = 0.04

power_law = lambda psi, beta=1.07, l1=10**22.06: l1*psi**beta

def broken_power_law(psi, mstar):
    """
    Calculate the radio luminosity using a broken power law model.

    Parameters:
    - psi (float): The input parameter.
    - mstar (float): The mass of the star.

    Returns:
    - lum (float): The calculated radio luminosity.
    """
    # Mass corrected with values from Gürkan et al. (in prep., LoTSS/H-Atlas radio lum-sfr relation)
    gamma = 0.44
    lc = 22.02
    beta_low = 0.52
    beta_high = 1.01
    psi_break = 0.01

    if psi <= psi_break:
        lum = lc * psi**beta_low * (mstar / 10**10 * u.solMass)**gamma
    else:
        lum = lc * psi**beta_high * (mstar / 10**10 * u.solMass)**gamma * psi_break**(beta_low - beta_high)

    return lum
