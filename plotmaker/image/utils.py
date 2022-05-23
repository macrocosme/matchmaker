import astropy.units as u

def set_image_radii(obj2, unit=u.arcsec, stretch=4,):
    # Let's face it, this is about lotss... Select either beam size or object size if bigger than beam
    beam = obj2.beam['bmaj'].to(unit).value
    return stretch * beam

def set_image_radii_from_obj_size(obj1, obj2, i, j, unit=u.arcsec, stretch=4,):
    beam = obj2.beam['bmaj'].to(unit).value
    semi_major = obj2.semi_major(j, to_unit=unit)
    obj2_size = semi_major if semi_major > beam else beam

    return obj1.semi_major(i, to_unit=unit) * stretch \
        if obj1.semi_major(i, to_unit=unit) > obj2_size \
        else obj2.semi_major(j, to_unit=unit) * stretch
