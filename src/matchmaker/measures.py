import math
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM #, WMAP9 as cosmo

luminosity_radio = lambda flux, lum_dist, z, alpha=-0.7: (flux * 4 * math.pi * lum_dist**2) / (1+z)**(1+alpha)

def powerlaw_scale(freq1, freq2, input1, alpha):
    '''flux at freq2 from powerlaw index alpha and flux1 at freq1

    Parameters
    ----------
    freq1:float
    freq2:float
    input1:float
        Can be flux, luminosity, etc.
    alpha:float

    Returns
    -------
    scaled_flux2:float
    '''
    return input1 * (freq2/freq1)**alpha

def w_hz_to_erg_s_hz(x):
    """
    Converts a value from watts per hertz to ergs per second per hertz.

    Parameters:
    x (float): The value to be converted.

    Returns:
    float: The converted value in ergs per second per hertz.
    """
    return (x * (u.W * u.Hz**-1)).to(u.erg * u.s**-1 * u.Hz**-1).value

def erg_s_hz_to_w_hz(x):
    """
    Converts energy flux density from erg/s/Hz to W/Hz.

    Parameters:
    - x: float or array-like
        Energy flux density in erg/s/Hz.

    Returns:
    - float or array-like
        Energy flux density in W/Hz.
    """
    return (x * (u.erg * u.s**-1 * u.Hz**-1)).to(u.W * u.Hz**-1).value

# TODO: these has fixed units -- should be provided through params instead

def luminosity_distance_series(redshifts:pd.Series,
                               fluxes:pd.Series,
                               output_units='W_Hz',
                               alpha=-0.7,
                               cosmo=FlatLambdaCDM(H0=70, Om0=0.3)
                               ):
    """
    Calculate the luminosity distance series for a given set of redshifts and fluxes.

    Parameters:
    - redshifts (pd.Series): Series of redshift values.
    - fluxes (pd.Series): Series of flux values.
    - output_units (str, optional): Output units for the luminosity distance. Default is 'W_Hz'.
    - alpha (float, optional): Spectral index. Default is -0.7.
    - cosmo (astropy.cosmology.FlatLambdaCDM, optional): Cosmology object. Default is FlatLambdaCDM(H0=70, Om0=0.3).

    Returns:
    - np.array: Array of calculated luminosity distances.

    """
    # Possibly deprecated. TBC.
    print (redshifts)
    print (fluxes)
    if output_units == 'W_Hz':
        return np.array([
            luminosity_radio(
                (f * u.mJy).to(u.W * u.m**-2 * u.Hz**-1),
                cosmo.luminosity_distance(z).to('m'),
                z,
                alpha=alpha
            ).value for f, z in zip(fluxes, redshifts)
        ])
    elif output_units == 'erg_s_Hz':
        return np.array(
            [
                luminosity_radio(
                    (f * u.mJy).to(u.W * u.m**-2 * u.Hz**-1),
                    cosmo.luminosity_distance(z).to('m'),
                    z,
                    alpha=alpha
                ).to(u.erg * u.s**-1 * u.Hz**-1).value for f, z in zip(fluxes, redshifts)
            ]
        )

def distance_atomic(redshift: float,
                    output_unit='Mpc',
                    cosmo=FlatLambdaCDM(H0=70, Om0=0.3),
                    with_unit=False):
    """
    Calculate the atomic distance based on the given redshift.

    Parameters:
    redshift (float): The redshift value.
    output_unit (str, optional): The desired output unit. Defaults to 'Mpc'.
    cosmo (astropy.cosmology.FlatLambdaCDM, optional): The cosmology model. Defaults to FlatLambdaCDM(H0=70, Om0=0.3).
    with_unit (bool, optional): Whether to include the unit in the output. Defaults to False.

    Returns:
    float or astropy.units.Quantity: The calculated atomic distance.

    """
    return cosmo.luminosity_distance(redshift).to(output_unit) if with_unit \
        else cosmo.luminosity_distance(redshift).to(output_unit).value

def luminosity_distance_atomic(redshift:float,
                               flux:float,
                               output_units='W_Hz',
                               alpha=-0.7,
                               cosmo=FlatLambdaCDM(H0=70, Om0=0.3),
                               with_unit=False):
    """
    Calculate the luminosity distance for atomic sources.

    Parameters:
    - redshift (float): The redshift of the source.
    - flux (float): The flux of the source.
    - output_units (str, optional): The desired output units. Default is 'W_Hz'.
    - alpha (float, optional): The spectral index. Default is -0.7.
    - cosmo (astropy.cosmology.FlatLambdaCDM, optional): The cosmology model. Default is FlatLambdaCDM(H0=70, Om0=0.3).
    - with_unit (bool, optional): Whether to return the result with units. Default is False.

    Returns:
    - The luminosity distance (float or astropy.units.Quantity) based on the specified parameters.
    """
    if output_units == 'W_Hz':
        ld = luminosity_radio(
            (flux * u.mJy).to(u.W * u.m**-2 * u.Hz**-1),
            cosmo.luminosity_distance(redshift).to('m'),
            redshift,
            alpha=alpha
        )
    elif output_units == 'erg_s_Hz':
        ld = luminosity_radio(
            (flux * u.mJy).to(u.W * u.m**-2 * u.Hz**-1),
            cosmo.luminosity_distance(redshift).to('m'),
            redshift,
            alpha=alpha
        ).to(u.erg * u.s**-1 * u.Hz**-1)
    else:
        ld = luminosity_radio(
            (flux * u.mJy).to(u.W * u.m**-2 * u.Hz**-1),
            cosmo.luminosity_distance(redshift).to('m'),
            redshift,
            alpha=alpha
        ).to(output_units)

    if with_unit:
        return ld
    else:
        return ld.value

def luminosity_distance(redshift,
                        flux,
                        alpha=-0.7,
                        output_units='W_Hz',
                        cosmo=FlatLambdaCDM(H0=70, Om0=0.3),
                        with_unit=False):
    """
    Calculate the luminosity distance for a given redshift and flux.

    Parameters:
    - redshift (float, np.array, pd.Series): The redshift value(s) for which to calculate the luminosity distance.
    - flux (float, np.array, pd.Series): The flux value(s) corresponding to the redshift(s).
    - alpha (float, optional): The spectral index. Default is -0.7.
    - output_units (str, optional): The desired output units for the luminosity distance. Default is 'W_Hz'.
    - cosmo (astropy.cosmology.FlatLambdaCDM, optional): The cosmology model to use for the calculation. Default is FlatLambdaCDM(H0=70, Om0=0.3).
    - with_unit (bool, optional): Whether to return the luminosity distance with units. Default is False.

    Returns:
    - luminosity_distance (float, np.array, pd.Series): The calculated luminosity distance(s).

    If the redshift parameter is not of type np.array or pd.Series, the function will call the luminosity_distance_atomic function with the provided parameters.
    If the redshift parameter is of type pd.Series, the function will call the luminosity_distance_atomic function with the values of the redshift and flux Series.
    If the redshift parameter is of type np.array, the function will call the luminosity_distance_atomic function with the redshift and flux converted to lists.
    """
    if type(redshift) not in [np.array, pd.Series]:
        return luminosity_distance_atomic(
            redshift=redshift,
            flux=flux,
            output_units=output_units,
            alpha=alpha,
            cosmo=cosmo,
            with_unit=with_unit
        )
    else:
        if type(redshift) is pd.Series:
            return luminosity_distance_atomic(
                redshift=redshift.values,
                flux=flux.values,
                output_units=output_units,
                alpha=alpha,
                cosmo=cosmo,
                with_unit=with_unit
            )
        elif type(redshift) is np.array:
            return luminosity_distance_atomic(
                redshift=list(redshift),
                flux=list(flux),
                output_units=output_units,
                alpha=alpha,
                cosmo=cosmo,
                with_unit=with_unit
            )
