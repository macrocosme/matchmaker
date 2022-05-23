import os
import numpy as np
import wget
from tqdm import tqdm

import astropy.units as u
from ..utils import set_image_radii
from ....matchmaker.utils import (
    check_folder_exists_or_create,
    remove_file_if_exists,
    remove_file_from_folder,
    get_matched
)
from ....matchmaker.dataset import Image, Column

# PanSTARSS
def get_name(ra, dec):
    rah = int(np.floor(ra/15))
    ram = int(np.floor(60*(ra/15-rah)))
    decd = int(np.floor(dec))
    decm = int(np.floor(60*(dec-decd)))
    return "BDC-J%02d%02d+%02d%02d"%(rah,ram,decd,decm)

arcsec_to_pixel = lambda arcsec: np.round(arcsec / 0.25).astype(int) # 0.25 arcsec / pixel in PanSTARRS

def download_images(obj1, obj2, mask1=None, mask2=None,
                    stretch=16,
                    base_url = 'http://ps1images.stsci.edu/cgi-bin',
                    output_path='data/panstarrs/',
                    only_missing=False,
                    check_files_are_in=False):

    matched_obj1 = obj1.df.iloc[mask1] if mask1 is not None else obj1.df

    if not only_missing:
        obj1.matches[obj2.name].missing_images = []

    for i, obj1_row in tqdm(matched_obj1.iterrows(), total=len(matched_obj1)):
        if only_missing:
            if i in obj1.matches[obj2.name].missing_images:
                download = True
            else:
                download = False
        else:
            download = True

        if download:
            try:
                ra = obj1_row[obj1.cols.ra.label]
                dec = obj1_row[obj1.cols.dec.label]

                size = arcsec_to_pixel(set_image_radii(obj2, unit=u.arcsec, stretch=stretch))
                check_folder_exists_or_create(output_path, return_folder=False)

                fname = wget.download("{}/ps1filenames.py?ra={:f}&dec={:f}".format(base_url, ra, dec),
                                      out=check_folder_exists_or_create('tmp/'))

                image = obj1.matches[obj2.name].images[i] = Image(source_catalog='panstarrs')
                with open(fname, "r") as f:
                    for line in f.readlines()[1:]:
                        w = line.strip("\n").split(" ")
                        filter = w[4]
                        red = w[7]
                        image.filters[filter] = Column()
                        image.filters[filter].remote_download = "{}/fitscut.cgi?" \
                                                                "ra={:f}&" \
                                                                "dec={:f}&" \
                                                                "red={}&" \
                                                                "wcs=true&" \
                                                                "format=fits&" \
                                                                "size={}".format(base_url, ra, dec, red, size)
                        try:
                            remove_file_if_exists(image.filters[filter].file_location)
                        except AttributeError:
                            if check_files_are_in:
                                image.filters[filter].file_location = wget.download(
                                    remove_file_if_exists(image.filters[filter].remote_download, return_filename=True),
                                    out=output_path
                                )
                            else:
                                pass

                        if not check_files_are_in:
                            image.filters[filter].file_location = wget.download(
                                remove_file_if_exists(image.filters[filter].remote_download, return_filename=True),
                                out=output_path
                            )

                    remove_file_from_folder(fname)
                    if i in obj1.matches[obj2.name].missing_images:
                        del obj1.matches[obj2.name].missing_images[np.where(obj1.matches[obj2.name].missing_images == i)[0][0]]

            except: # TODO: except too general: should be only when download failed
                obj1.matches[obj2.name].missing_images.append(i)

