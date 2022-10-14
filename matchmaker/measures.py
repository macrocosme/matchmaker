import math
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

luminosity_radio = lambda flux, lum_dist, z, alpha=-0.7: (flux * 4 * math.pi * lum_dist**2) / (1+z)**(1+alpha)


def w_hz_to_erg_s_hz(x):
    return (x * (u.W * u.Hz**-1)).to(u.erg * u.s**-1 * u.Hz**-1).value

def erg_s_hz_to_w_hz(x):
    return (x * (u.erg * u.s**-1 * u.Hz**-1)).to(u.W * u.Hz**-1).value

# TODO: these has fixed units -- should be provided through params instead

def luminosity_distance_series(redshifts:pd.Series,
                               fluxes:pd.Series,
                               output_units='W_Hz',
                               alpha=-0.7,
                               cosmo=FlatLambdaCDM(H0=70, Om0=0.3)):
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

def luminosity_distance_atomic(redshift:float,
                               flux:float,
                               output_units='W_Hz',
                               alpha=-0.7,
                               cosmo=FlatLambdaCDM(H0=70, Om0=0.3),
                               with_unit=False):

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
