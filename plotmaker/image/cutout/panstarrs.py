import astropy.units as u
import numpy as np
import wget
from tqdm import tqdm

from ....matchmaker.dataset import Column, Image
from ....matchmaker.utils import check_folder_exists_or_create, remove_file_from_folder, remove_file_if_exists
from ..utils import set_image_radii


# PanSTARSS
def get_name(ra, dec):
    rah = int(np.floor(ra/15))
    ram = int(np.floor(60*(ra/15-rah)))
    decd = int(np.floor(dec))
    decm = int(np.floor(60*(dec-decd)))
    return f"BDC-J{rah:02d}{ram:02d}+{decd:02d}{decm:02d}"

arcsec_to_pixel = lambda arcsec: np.round(arcsec / 0.25).astype(int) # 0.25 arcsec / pixel in PanSTARRS

def download_images(obj1, obj2, mask1=None, mask2=None,
                    stretch=16,
                    base_url = 'http://ps1images.stsci.edu/cgi-bin',
                    output_path='data/panstarrs/',
                    only_missing=False):

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
            ra = obj1_row[obj1.cols.ra.label]
            dec = obj1_row[obj1.cols.dec.label]

            size = arcsec_to_pixel(set_image_radii(obj2, unit=u.arcsec, stretch=stretch))
            check_folder_exists_or_create(output_path, return_folder=False)

            fname = wget.download(f"{base_url}/ps1filenames.py?ra={ra:f}&dec={dec:f}",
                                  out=check_folder_exists_or_create('tmp/'))

            image = obj1.matches[obj2.name].images[i] = Image(source_catalog='panstarrs')

            with open(fname) as f:
                for line in f.readlines()[1:]:
                    w = line.strip("\n").split(" ")
                    filter = w[4]
                    red = w[7]
                    image.filters[filter] = Column()
                    image.filters[filter].remote_download = f"{base_url}/fitscut.cgi?" \
                                                            f"ra={ra:f}&" \
                                                            f"dec={dec:f}&" \
                                                            f"red={red}&" \
                                                            "wcs=true&" \
                                                            "format=fits&" \
                                                            f"size={size}"

                    image.filters[filter].file_location = wget.download(
                        remove_file_if_exists(image.filters[filter].remote_download, return_filename=True),
                        out=f"{output_path}{i}.fits"
                    )

                remove_file_from_folder(fname)
                if i in obj1.matches[obj2.name].missing_images:
                    del obj1.matches[obj2.name].missing_images[np.where(obj1.matches[obj2.name].missing_images == i)[0][0]]


