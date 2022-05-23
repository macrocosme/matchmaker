import astropy.units as u
import pandas as pd
from astropy.coordinates import SkyCoord
from . import DATA_BASE_PATH

def frb121102_wise_image_files():
    wise_bands = ['w1', 'w2', 'w3', 'w4']
    frb121102_wise_files = {
        band: DATA_BASE_PATH + 'data/FRB121102B_WISE_Files/L3a/0827p333_ab41/0827p333_ab41-{}-int-3_ra83.0375_dec33.086944_asec600.000.fits'.format(band) for band in wise_bands
    }
    return frb121102_wise_files

df_frb121102_wise = pd.read_csv(DATA_BASE_PATH + 'data/FRB_WISE_Files/FRB121102_WISE.csv')
frb121102 = {
    # http://simbad.u-strasbg.fr/simbad/sim-coo?protocol=html&NbIdent=1&Radius=1&Radius.unit=arcmin&CooFrame=FK5&CooEpoch=2000&CooEqui=2000&Coord=83.0375%2B33.086944
        'coords': SkyCoord(ra=83.04 * u.deg,
                           dec=33.087056 * u.deg),
        'freq': 1.4 * u.GHz,
        'w1': df_frb121102_wise['w1mag'].values[0],
        'w2': df_frb121102_wise['w2mag'].values[0],
        'w3': df_frb121102_wise['w3mag'].values[0],
        'w4': df_frb121102_wise['w4mag'].values[0],
        'w1_res': 6.1 * u.arcsec,
        'w2_res': 6.4 * u.arcsec,
        'w3_res': 6.5 * u.arcsec,
        'w4_res': 12.0 * u.arcsec,
        # 'lum' : 1.6e29 * (u.erg * u.s**-1 * u.Hz**-1),  # Eftekhari+2020 (10 GHz)
        'lum' : 2e29 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1),  # Law Connor Aggarwal 2021 (1.4 GHz)
        'flux': 180 * u.uJy, # Chatterjee+ 2017
        'z_host': 0.32,
        'Msun' : 1e8,  # Bassa+ 2017
        'sfr' : 0.8,  # Tendulkar 2017 quotes  0.4 solarmass per year
        'alpha' : -0.27, # Marcote et al (2017) -- EVN   (Spectral index)
        'alpha_err': 0.24, # Marcote et al (2017) -- EVN
        'color' : 'purple',
        's' : 50,
        'name': 'FRB 20121102A',
        'label': 'FRB 20121102A (10 GHz)',
    }

df_frb190520b_wise = pd.read_csv(DATA_BASE_PATH + 'data/FRB_WISE_Files/FRB20190520B_WISE.csv')
frb190520b = {
    # http://simbad.u-strasbg.fr/simbad/sim-id?Ident=FRB+20190520B&NbIdent=1&Radius=2&Radius.unit=arcmin&submit=submit+id
    'coords': SkyCoord(ra=240.5178 * u.deg,
                       dec=-11.2881444 * u.deg),
    'freq': 3 * u.GHz,
    'w1': df_frb190520b_wise['W1mag'].values[0],
    'w2': df_frb190520b_wise['W2mag'].values[0],
    'w3': df_frb190520b_wise['W3mag'].values[0],
    'w4': df_frb190520b_wise['W4mag'].values[0],
    'lum' : 3e29 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1),  # From Niu+ 2022
    'flux': 202 * u.uJy,
    'flux_err': 8 * u.uJy,
    'z_host': 0.241,
    'z_host_err': 0.001,
    'Msun': 6e8,  # From Niu+ 2022
    'sfr': 0.41, # From Niu+ 2022
    'alpha': -0.41, # From Zhao & Wang 2021
    'alpha_err': 0.04, # From Zhao & Wang 2021
    'color' : 'red',
    's' : 50,
    'name': 'FRB 20190520B',
    'label': 'FRB 20190520B (3 GHz)'
}

SGR1935_2154 = {
    'lum' : 2.61e27 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1), #STAR2 Chime (Wang, Zhang, Dai 2020)
    'halo_mass' : 1.5e12,
    'sfr' :  1.06, # Average of RW2010
    'sfr_l' : 0.68, # Robitaille & Whitney 2010
    'sfr_h' : 1.45,
    'color' : 'darkgrey',
    's' : 100,
    'label': 'SGR 1935+2154'
}
FRBs = [frb121102, frb190520b]#, SGR1935_2154]

l_ptf10hgi = 1.1e28 * (u.erg * u.s**-1 * u.Hz**-1).to(u.W * u.Hz**-1)
m_lmc = 3e9
