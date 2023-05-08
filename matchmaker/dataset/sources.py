import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from . import DATA_BASE_PATH

pink = (230/255, 29/255, 95/255, 1)
pink_translucid = (230/255, 29/255, 95/255, .2)
blue = (47/255, 161/255, 214/255, 0.2)
blue_full = (47/255, 161/255, 214/255, 1)

frb121102 = {
    # http://simbad.u-strasbg.fr/simbad/sim-coo?protocol=html&NbIdent=1&Radius=1&Radius.unit=arcmin&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&Coord=83.0375%2B33.086944
        'coords': SkyCoord(ra=83.04 * u.deg,
                           dec=33.087056 * u.deg),
        'freq': 1.4 * u.GHz,
        # 'lum' : 1.6e29 * (u.erg * u.s**-1 * u.Hz**-1),  # Eftekhari+2020 (10 GHz)
        'lum' : 2e29 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1),  # Law Connor Aggarwal 2021 (1.4 GHz)
        'flux': 180 * u.uJy, # Chatterjee+ 2017
        'z_host': 0.32,
        'Msun' : 1.3e8,  # Bassa+ 2017
        'mstar_err': 0.4e8, # Bassa+ 2017
        'sfr' : 0.4,  # Tendulkar 2017 quotes  0.4 solarmass per year
        'sfr_err': 0.17,
        # 'alpha' : -0.27, # Marcote et al (2017) -- EVN   (Spectral index)
        # 'alpha_err': 0.24, # Marcote et al (2017) -- EVN
        'alpha': -0.07, # Resmi+ (2021)
        'alpha_err': 0.03, # Resmi+ (2021)
        'color' : pink,
        's' : 100,
        'zorder': 90,
        'name': 'FRB 20121102A',
        'label': 'FRB 20121102A',
    }

frb190520b = {
    # http://simbad.u-strasbg.fr/simbad/sim-id?Ident=FRB+20190520B&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id
    'coords': SkyCoord(ra=240.5178 * u.deg,
                       dec=-11.2881444 * u.deg),
    'freq': 3 * u.GHz,
    'lum' : 3e29 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1),  # From Niu+ 2022
    'flux': 202 * u.uJy,
    'flux_err': 8 * u.uJy,
    'z_host': 0.241,
    'z_host_err': 0.001,
    'Msun': 6e8,  # From Niu+ 2022
    'mstar_err': 0,
    'sfr': 0.41, # From Niu+ 2022
    'sfr_err': 0., # From Niu+ 2022
    'alpha': -0.41, # From Zhao & Wang 2021
    'alpha_err': 0.04, # From Zhao & Wang 2021
    'color' : blue_full,
    's' : 100,
    'zorder': 90,
    'name': 'FRB 20190520B',
    'label': 'FRB 20190520B'
}

SGR1935_2154 = {
    'lum' : 2.61e27 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1), #STAR2 Chime (Wang, Zhang, Dai 2020)
    'halo_mass' : 1.5e12,
    'sfr' :  1.06, # Average of RW2010
    'sfr_l' : 0.68, # Robitaille & Whitney 2010
    'sfr_h' : 1.45,
    'color' : 'darkgrey',
    's' : 50,
    'zorder': 90,
    'label': 'SGR 1935+2154'
}
FRBs = [frb121102, frb190520b]#, SGR1935_2154]

l_ptf10hgi = 1.1e28 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1)
m_lmc = 3e9
