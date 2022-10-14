import astropy.units as u

def construct_image_url(dataset, idx, size=0.4, maxint=50, fits=1):
    """Construct the image url.
    Parameters
    ----------
    dataset:
        matchmaker.matchmaker.dataset.Dataset
    idx:
        object index in dataset
    size:
        Size in arcmin
    maxint:
        If the image is returned as a GIF image, this parameter determines the contrast enhancement.
    fits:
        If 1, returns the image as a fits file. If 0, returns image as GIF.
    :return:
    """
    ra = dataset.as_SkyCoord(idx).ra.to_string(u.hour, sep=":", precision=5)
    dec = dataset.as_SkyCoord(idx).dec.to_string(u.degree, alwayssign=True, sep=":", precision=5)

    url = "https://third.ucllnl.org/cgi-bin/firstimage?" \
          "RA={}&" \
          "Dec={}&Equinox=J2000&" \
          "ImageSize={}&" \
          "MaxInt={}&" \
          "FITS={}".format(ra, dec, size, maxint, fits)

    return url
