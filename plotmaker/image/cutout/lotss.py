from ..utils import set_image_radii
from ....matchmaker.utils import check_folder_exists_or_create, get_matched
from ....matchmaker.dataset import Image
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord as SkyCoord
import numpy as np
import pandas as pd

# Example use to get cutout from coord:
# python radio_cutout_matchmaker.py 9.59220833333 32.63822222222 P010+34 652 0.038412;

def make_lotss_cutouts_download_script(obj1, obj_lotss, mask1=None, mask_lotss=None,
                                       outpath='cutout', print_commands=True, saveto=None, saveto_folder='data/lotss'):
    matched_obj1, matched_lotss = get_matched(obj1, obj_lotss, mask1, mask_lotss)

    # TODO: the split with repetition here is cumbersome and prone to future errors... find a better way.
    if saveto is not None:
        save = '{}/{}'.format(check_folder_exists_or_create(saveto_folder), saveto)
        with open(save, 'w+') as f:
            for (i, row_obj1), (j, row_lotss) in zip(matched_obj1.iterrows(), matched_lotss.iterrows()):
                line = 'python radio_cutout_matchmaker.py {} {} {} {} {:f} {};'.format(
                    row_obj1[obj1.cols.ra.label],
                    row_obj1[obj1.cols.dec.label],
                    row_lotss['Mosaic_ID'],
                    i,
                    set_image_radii(obj_lotss, unit=u.deg),
                    outpath
                )
                f.write(line + '\n')
                if print_commands:
                    print(i, j, line)
                # Temporary trick as we need to go through different servers independently to gather these images:
                obj_lotss.matches[obj1.name].images[i] = Image(source_catalog='lotss dr2')
                obj_lotss.matches[obj1.name].images[i].file_location = 'data/lotss/{}/{}_cutout.fits'.format(outpath, i)

        print ('Saved to {}'.format(save))
    else:
        if print_commands:
            for (i, row_obj1), (j, row_lotss) in zip(matched_obj1.iterrows(), matched_lotss.iterrows()):
                print('python radio_cutout_matchmaker.py {} {} {} {} {:f} {};'.format(
                    row_obj1[obj1.cols.ra.label],
                    row_obj1[obj1.cols.dec.label],
                    row_lotss['Mosaic_ID'],
                    i,
                    set_image_radii(obj_lotss, unit=u.deg),
                    outpath
                ))
                # Temporary trick as we need to go through different servers independently to gather these images:
                obj_lotss.matches[obj1.name].images[i] = Image(source_catalog='lotss dr2')
                obj_lotss.matches[obj1.name].images[i].file_location = 'data/lotss/{}/{}_cutout.fits'.format(outpath, i)


def make_lotss_cutouts_obj1_download_script(obj1, obj_lotss, mask1=None,
                                            outpath='cutout', print_commands=True,
                                            core_collapse_only=False,
                                            stretch=4,
                                            saveto=None, saveto_folder='data/lotss'):
    matched_obj1 = obj1.df.iloc[mask1] if mask1 is not None else obj1.df
    if core_collapse_only:
        # This is a temporary use case that should be done some other, generic, way
        matched_obj1 = matched_obj1.loc[
            # matched_obj1[obj1.cols.tns.obj_type.label] != 'SN Ia'
            (matched_obj1[obj1.cols.tns.obj_type.label] != 'SN Ia') &
            (matched_obj1['Redshift'] <= 0.005) &
            (~matched_obj1['Redshift'].isna())
        ]


    # print (np.unique(matched_obj1[obj1.cols.tns.obj_type.label]))

    # TODO: the split with repetition here is cumbersome and prone to future errors... find a better way.
    if saveto is not None:
        save = '{}/{}'.format(check_folder_exists_or_create(saveto_folder), saveto)
        with open(save, 'w+') as f:
            for (i, row_obj1) in matched_obj1.iterrows():

                field, obs_date = obj_lotss.get_closest_field(row_obj1[obj1.cols.ra.label],
                                                    row_obj1[obj1.cols.dec.label],
                                                    row_obj1[obj1.cols.tns.discovery_date_ut.label])
                if field is not None:
                    print (i, row_obj1['Name'], row_obj1['ra'], row_obj1['dec'], obs_date)
                    line = 'python radio_cutout_byname.py "{}" {};'.format(
                        row_obj1[obj1.cols.tns.name.label],
                        field
                    )

                    line = 'python radio_cutout_matchmaker.py {} {} {} {} {:f} {};'.format(
                        row_obj1[obj1.cols.ra.label],
                        row_obj1[obj1.cols.dec.label],
                        field,
                        i,
                        set_image_radii(obj_lotss, unit=u.deg, stretch=stretch),
                        outpath
                    )
                    f.write(line + '\n')
                    if print_commands:
                        print(i, line)
                    # Temporary trick as we need to go through different servers independently to gather these images:
                    obj_lotss.matches[obj1.name].images[i] = Image(source_catalog='lotss dr2')
                    obj_lotss.matches[obj1.name].images[i].file_location = 'data/lotss/{}/{}_cutout.fits'.format(outpath, i)

        print ('Saved to {}'.format(save))

def make_cutout_from_field_fits(ra, dec, filename, size, name, outpath):
    # ra = float(sys.argv[1])
    # dec = float(sys.argv[2])
    # field = sys.argv[3]         # LOTSS field name
    # name = sys.argv[4]
    # size = float(sys.argv[5]) # in degree
    # outpath = sys.argv[6]

    # filename = '/disks/paradata/shimwell/LoTSS-DR2/archive_DR2_final/'+field+'/image_full_ampphase_di_m.NS_shift.int.facetRestored.fits'
    # size = float(0.5/60.0)
    h = fits.getheader(filename)

    print (h['CDELT2'])

    # position = SkyCoord(ra=ra*u.deg, dec = dec*u.deg, frame='icrs')

    h['NAXIS']=2; h['WCSAXES']=2
    del h['NAXIS3'], h['NAXIS4'], h['CRPIX3'], h['CDELT3'], h['CRPIX4'], h['CDELT4'], h['CRVAL3'], h['CRVAL4'], h['CTYPE3'], h['CTYPE4'], h['CUNIT4']
    wcs = WCS(h)

    dim = fits.getdata(filename)[-1, -1, :, :]

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

    fits.writeto(data = data_out,
                 header = h,
                 filename = "%s/%s_cutout.fits" % (outpath, name),
                 overwrite = True)
