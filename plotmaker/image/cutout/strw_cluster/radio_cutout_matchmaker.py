import os
import sys
import numpy as np
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.coordinates import SkyCoord as SkyCoord

# example usage
# python radio_cutout_matchmaker.py 9.59220833333 32.63822222222 P010+34 <dataframe id> 0.038412;
#
#
if not len(sys.argv)==7:
    print ("ARGUMENTS: ra, dec, field, name, size, outpath")
    exit(1)
    
ra = float(sys.argv[1])
dec = float(sys.argv[2])
field = sys.argv[3]         # LOTSS field name  
name = sys.argv[4]
size = float(sys.argv[5]) # in degree
outpath = sys.argv[6]

filename = '/disks/paradata/shimwell/LoTSS-DR2/archive_DR2_final/'+field+'/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
# size = float(0.5/60.0)
try:
    h = fits.getheader(filename)
    dim = fits.getdata(filename)[-1, -1, :, :]
except IOError:
    filename = '/disks/paradata/shimwell/Beyond-DR2/archive_images/'+field+'/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
    h = fits.getheader(filename)
    dim = fits.getdata(filename)[-1, -1, :, :]

position = SkyCoord(ra=ra*u.deg, dec = dec*u.deg, frame='icrs')

h['NAXIS']=2; h['WCSAXES']=2; del h['NAXIS3'], h['NAXIS4'], h['CRPIX3'], h['CDELT3'], h['CRPIX4'], h['CDELT4'], h['CRVAL3'], h['CRVAL4'], h['CTYPE3'], h['CTYPE4'], h['CUNIT4']
wcs = WCS(h)

pix = wcs.wcs_world2pix(ra, dec, 0)

x_center, y_center = int(pix[0]), int(pix[1])

npix = int(size / np.absolute(h['CDELT2'] / 2))
tp = h['CRPIX1']
h['CRPIX1'] = tp - x_center + npix

tp = h['CRPIX2']
h['CRPIX2'] = tp - y_center + npix

h['NAXIS1'], h['NAXIS2'] = 2*npix, 2*npix

data_out = dim[y_center-npix:y_center+npix, 
               x_center-npix:x_center+npix]


def check_folder_exists_or_create(folder, return_folder=True):
    if not os.path.exists(folder):
        os.makedirs(folder)
    if return_folder:
        return folder



fits.writeto(data = data_out, 
             header = h, 
             filename = "%s/%s_cutout.fits" % (check_folder_exists_or_create(outpath), name),
             overwrite = True) 



