from astroquery.sdss import SDSS
from astroquery.esasky import ESASky
from astroquery.vizier import Vizier
import astropy.units as u

def search_sdss(sky_coords, mask, get_spectra=True, get_images=False):
    sdss_matches = {}
    for ii in mask:
        # print (ii, coord)
        xid = SDSS.query_region(sky_coords[ii], spectro=True)
        if xid is not None:
            sdss_matches[ii] = {}
            sdss_matches[ii]['table'] = xid
            if get_spectra:
                sdss_matches[ii]['spectra'] = SDSS.get_spectra(matches=xid)
            if get_images:
                sdss_matches[ii]['images'] = SDSS.get_images(matches=xid, band=['u', 'g', 'r', 'i', 'z'])
    return sdss_matches

def search_ESASky(obj, mask):
    matches = {}

    catalogs=['CHANDRA-SC2',
              'XMM-EPIC',
              'XMM-EPIC-STACK',
              'LAMOST',
              'HSC',
              'Herschel-HPPSC-070',
              'Herschel-HPPSC-100',
              'Herschel-HPPSC-160',
              'Herschel-SPSC-250',
              'Herschel-SPSC-350',
              'Herschel-SPSC-500',
              # 'Fermi_4FGL-DR2'
             ]

    for i in mask:
        radius = obj.df.iloc[i][obj.cols.major.label] * obj.cols.major.unit
        returned_table = ESASky.query_region_catalogs(position=obj.as_SkyCoord()[i],
                                                      radius=radius,
                                                      catalogs=catalogs)
        if len(returned_table):
            matches[i] = returned_table
            print (i, radius, matches[i])

    return matches

def search_vizier(obj, idx):
    matches = {}
    for i in idx:
        # radius = obj.df.iloc[i][obj.cols.major.label] * obj.cols.major.unit
        # print (radius)
        radius = 1 * u.arcsec
        matches[i] = result = Vizier.query_region(obj.as_SkyCoord(mask=i),
                                                  radius=radius)
    return matches

