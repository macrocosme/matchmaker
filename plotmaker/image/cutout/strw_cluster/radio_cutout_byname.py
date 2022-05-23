from astropy import units as u
from astropy.wcs import WCS
import sys
import astropy.io.fits as pf
from astropy.coordinates import SkyCoord as SkyCoord
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes import SphericalCircle
from astroquery.gaia import Gaia
from astropy.coordinates import get_icrs_coordinates
from astropy.time import Time
import os
#
#
#
if not len(sys.argv)==3:
   print ("ARGUMENTS: obj-name, field")
   exit(1)
#
field = sys.argv[2]         # LOTSS field name
name=sys.argv[1]            # Common name for annotation on plot


f = open("fields.txt","r")
lines = f.readlines()
f.close()
field_names = []
isot_times = []
for line in lines:
   w = line.strip("\n").split(",")
   field_names.append(w[0])
   isot_times.append(w[4].replace(" ","T"))

#
try:
   coord = get_icrs_coordinates(name)
except:
   print ("Could not find object!")
   exit(1)
print ("Found object at coordinates %f, %f"%(coord.ra.value,coord.dec.value))
#	
ra = coord.ra.value
dec = coord.dec.value
size = float(0.5/60.0)
try:
   fin_name = '/disks/paradata/shimwell/LoTSS-DR2/archive_DR2_final/'+field+'/image_full_ampphase_di_m.NS_shift.app.facetRestored.fits'
   h = pf.getheader(fin_name)
   dim = pf.getdata(fin_name)[-1,-1,:,:]
except:
   fin_name = '/disks/paradata/shimwell/Beyond-DR2/archive_images/'+field+'/image_full_ampphase_di_m.NS_shift.app.facetRestored.fits'
   h = pf.getheader(fin_name)
   dim = pf.getdata(fin_name)[-1,-1,:,:]
#
#
h['NAXIS']=2; h['WCSAXES']=2; del h['NAXIS3'], h['NAXIS4'], h['CRPIX3'], h['CDELT3'], h['CRPIX4'], h['CDELT4'], h['CRVAL3'], h['CRVAL4'], h['CTYPE3'], h['CTYPE4'], h['CUNIT4']
wcs = WCS(h)

pix = wcs.wcs_world2pix(ra,dec,0)
x_cent = int(pix[0]); y_cent = int(pix[1])
npix = int(size/np.absolute(h['CDELT2']/2))
tp = h['CRPIX1']; h['CRPIX1'] = tp-x_cent+npix
tp = h['CRPIX2']; h['CRPIX2'] = tp-y_cent+npix
h['NAXIS1'] = 2*npix; h['NAXIS2'] = 2*npix
data_out = dim[y_cent-npix:y_cent+npix, x_cent-npix:x_cent+npix]
stat_data = dim[y_cent-100:y_cent+100, x_cent-100:x_cent+100]*1e3
pf.writeto(data = data_out, header = h, filename = "cutout.fits", overwrite = True) 
med = np.nanmedian(stat_data)
mad = np.nanmedian(np.absolute(stat_data-med))

plt.figure(figsize=(2.5,2.5))
hdu = pf.open("cutout.fits")[0]
h = hdu.header
wcs = WCS(h)
ax = plt.subplot(projection=wcs)


ax.imshow(hdu.data*1e3, vmin=med-7*mad, vmax=med+7*mad, origin='lower',cmap="viridis")
#ax.imshow(hdu.data*1e3, vmin=-vlim, vmax=vlim, origin='lower',cmap="viridis")
c = SphericalCircle((ra*u.deg, dec*u.deg),5./3600*u.deg,
                     edgecolor='black', facecolor='none',
                     transform=ax.get_transform('fk5'))
ax.add_patch(c)
plt.gca().set_axis_off()
plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
#plt.annotate(text=name,xy=(0.3,0.9),xycoords="figure fraction",fontsize=16)
#II = field_names.index(field)
#t = Time(isot_times[II],scale="utc",format="isot")
#print (isot_times[II],t.jd)
#plt.annotate(text=r"JD$_{\rm mid}$ = %.2f"%(t.jd+4./24),xy=(0.1,0.1),xycoords="figure fraction",fontsize=16)
plt.margins(0, 0)
plt.gca().xaxis.set_major_locator(plt.NullLocator())
plt.gca().yaxis.set_major_locator(plt.NullLocator())
plt.savefig("%s-%s.png"%(name.replace(" ","-"),field), transparent=True, bbox_inches="tight", pad_inches=0)
plt.close()
