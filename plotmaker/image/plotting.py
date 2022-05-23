import numpy as np
import pandas as pd
from astropy.io import fits
import astropy.units as u
import aplpy
from tqdm import tqdm

from .cutout.panstarrs import get_name
from .utils import set_image_radii
from ...matchmaker.utils import check_folder_exists_or_create, get_matched

pink_rgba = (230/255, 29/255, 95/255, 1)
pink_rgb = (230/255, 29/255, 95/255)
pink_translucid = (230/255, 29/255, 95/255, .2)
blue = (47/255, 161/255, 214/255, 0.2)
blue_full = (47/255, 161/255, 214/255, 1)
blue_rgb = (47/255, 161/255, 214/255)

def collect_statistics(filename):
    d = fits.open(filename)
    med = np.nanmedian(d[0].data)
    mad = np.nanmedian(np.absolute(d[0].data-med))
    min = np.nanmin(d[0].data)
    return d, med, mad, min

def get_beam_info(filename):
    hh = fits.getheader(filename)
    bmaj, bmin, bpa = hh['BMAJ'], hh['BMIN'], hh['BPA']
    return bmaj, bmin, bpa

def plot_panstarrs_filter(obj1, obj2,
                          mask1=None, mask2=None,
                          d2d_obj1=None, d2d_obj2=None,
                          stretch=4,
                          filter='r', extra_path=None):
    # PanSTARRS r filter fits most of CLU's 4 filters (H-alpha 656, 663, 672, 681)

    matched_obj1, matched_obj2 = get_matched(obj1, obj2, mask1, mask2)

    for (i, row_obj1), (j, row_obj2) in tqdm(zip(matched_obj1.iterrows(), matched_obj2.iterrows()),
                                             total=len(matched_obj1)):
        # Obj1 DATA
        filename_obj1 = obj1.matches[obj2.name].images[i].filters[filter].file_location
        d, med, mad, _ = collect_statistics(filename_obj1)
        gc = aplpy.FITSFigure(d, figsize=(4,4))
        gc.show_colorscale(vmin = med-15*mad,
                           vmax = med+15*mad,
                           interpolation = 'none',
                           cmap = 'Greys')
        del d

        gc.recenter(row_obj1[obj1.cols.ra.label],
                    row_obj1[obj1.cols.dec.label],
                    radius=set_image_radii(obj2, unit=u.deg, stretch=stretch))

        # Obj2 DATA
        filename_obj2 = obj2.matches[obj1.name].images[i].file_location
        print (i, j, filename_obj2)
        dd, dd_med, dd_mad, _  = collect_statistics(filename_obj2)
        bmaj, bmin, bpa = get_beam_info(filename_obj2)
        gc.show_contour(dd,
                        colors=['#7fcdbb', '#edf8b1', '#c7e9b4', '#ffffd9'],
                        labels=[r'$\mu+3\sigma$', r'$\mu+5\sigma$', r'$\mu+10\sigma$', r'$\mu+15\sigma$'],
                        alpha=0.75,
                        levels=[dd_med + 3*dd_mad,
                                dd_med + 5*dd_mad,
                                dd_med + 10*dd_mad,
                                dd_med + 15*dd_mad]
                        )
        del dd
        # gc._ax1.legend()

        try:
            gc.show_ellipses(row_obj1[obj1.cols.ra.label],
                             row_obj1[obj1.cols.dec.label],
                             width=obj1.semi_major(i, to_unit=u.deg) * 2,
                             height=obj1.semi_minor(i, to_unit=u.deg) * 2,
                             angle=-row_obj1[obj1.cols.pa.label],
                             edgecolor=pink_rgb,
                             linestyle='-',
                             facecolor='none',
                             zorder=100)
        except AttributeError:
            pass

        gc.add_beam(major=bmaj, minor=bmin, angle=bpa)
        gc.beam.set_edgecolor('black')
        gc.beam.set_facecolor('none')
        gc.beam.set_hatch('+')

        gc.show_markers(row_obj1[obj1.cols.ra.label],
                        row_obj1[obj1.cols.dec.label],
                        marker='*',
                        facecolor='#edf8b1',
                        edgecolor='#edf8b1',
                        alpha=0.75,
                        zorder=100)

        gc.show_markers(row_obj2[obj2.cols.ra.label],
                        row_obj2[obj2.cols.dec.label],
                        marker='+',
                        alpha=0.75,
                        edgecolor=blue_rgb, 
                        facecolor=blue_rgb)

        gc.show_ellipses(row_obj2[obj2.cols.ra.label],
                         row_obj2[obj2.cols.dec.label],
                         width=(row_obj2[obj2.cols.ra_err.label] * obj2.cols.ra_err.unit).to(u.deg) * 2,
                         height=(row_obj2[obj2.cols.dec_err.label] * obj2.cols.dec_err.unit).to(u.deg) * 2,
                         angle=0,
                         edgecolor=blue_rgb,
                         linestyle='--',
                         facecolor='none')

        gc.show_ellipses(row_obj2[obj2.cols.ra.label],
                         row_obj2[obj2.cols.dec.label],
                         width=obj2.semi_major(j, to_unit=u.deg) * 2,
                         height=obj2.semi_minor(j, to_unit=u.deg) * 2,
                         angle=-row_obj2[obj2.cols.pa.label],
                         edgecolor=blue_rgb,
                         linestyle='-',
                         facecolor='none')


        # Find all other galaxies in image frame
        # try:
        #     image_radii = set_image_radii(obj1, obj2, i, j, u.arcsec)
        #     # CASE clu
        #     if d2d_obj1 is None:
        #         local_d2d = obj1.as_SkyCoord(i).separation(obj1.as_SkyCoord())
        #     else:
        #         local_d2d = d2d_obj1[i]['d2d']
        #     test_mask = (local_d2d <= image_radii * u.arcsec)
        #
        #     if len(test_mask) > 1 or d2d_obj1 is None:
        #         for k, match in obj1.df.iloc[
        #             test_mask if d2d_obj1 is None else d2d_obj1[i]['filtered_idx'][test_mask]
        #         ].iterrows():
        #             if k != i:
        #                 gc.show_ellipses(match[obj1.cols.ra.label],
        #                                  match[obj1.cols.dec.label],
        #                                  width=obj1.semi_major(k, to_unit=u.deg) * 2,
        #                                  height=obj1.semi_minor(k, to_unit=u.deg) * 2,
        #                                  angle=-match[obj1.cols.pa.label] if not pd.isnull(match[obj1.cols.pa.label]) else 0,
        #                                  edgecolor='darkgrey',
        #                                  alpha=0.65,
        #                                  linestyle='dotted',
        #                                  facecolor='none')
        # except AttributeError:
        #     pass
        #
        # # CASE lotss
        # if d2d_obj2 is None:
        #     local_d2d = obj1.as_SkyCoord(i).separation(obj2.as_SkyCoord())
        # else:
        #     local_d2d = d2d_obj2[i]['d2d']
        # test_mask = (local_d2d <= image_radii * u.arcsec)
        # if len(test_mask) > 1 or d2d_obj2 is None:
        #     for k, match in obj2.df.iloc[
        #         test_mask if d2d_obj2 is None else d2d_obj2[i]['filtered_idx'][test_mask]
        #     ].iterrows():
        #         if k != j:
        #             gc.show_ellipses(match[obj2.cols.ra.label],
        #                              match[obj2.cols.dec.label],
        #                              width=obj2.semi_major(k, to_unit=u.deg) * 2,
        #                              height=obj2.semi_minor(k, to_unit=u.deg) * 2,
        #                              angle=-match[obj2.cols.pa.label],
        #                              edgecolor='grey',
        #                              alpha=0.7,
        #                              facecolor='none')

        # gc.set_title('LoTSS df id: %d' %j)
        gc.set_title('LoTSS df id: {}, {}{}'.format(
            j,
            row_obj1[obj1.cols.name.label],
            ':Type {}'.format(row_obj1[obj1.cols.obj_type.label]) if obj1.prop('obj_type') is not None else ''
        ))

        gc.set_theme('publication')

        outpath = check_folder_exists_or_create('images/matches/{}'.format(filter) if extra_path is None else 'images/matches/{}/{}'.format(extra_path, filter),
                                                return_folder=True)
        # gc.save("%s/%s_%s.pdf" % (outpath, j, get_name(row_obj1[obj1.cols.ra.label],
        #                                                row_obj1[obj1.cols.dec.label])), dpi=250)
        outfile = "%s/%s_%s.pdf" % (outpath, row_obj1[obj1.cols.name.label], j)
        print ('Saving to {}'.format(outfile))
        gc.save(outfile, dpi=250)

        gc.close()

def plot_obj1_img_obj2_contour(obj1, obj2,
                          mask1=None, mask2=None,
                          d2d_obj1=None, d2d_obj2=None,
                          stretch=4,
                          filter='r', extra_path=None):
    # PanSTARRS r filter fits most of CLU's 4 filters (H-alpha 656, 663, 672, 681)

    matched_obj1 = obj1.df.iloc[mask1] if mask1 is not None else obj1.df

    for i, row_obj1 in tqdm(matched_obj1.iterrows(), total=len(matched_obj1)):
        try:
            # Obj1 DATA
            filename_obj1 = obj1.matches[obj2.name].images[i].filters[filter].file_location
            d, med, mad, min = collect_statistics(filename_obj1)
            gc = aplpy.FITSFigure(d, figsize=(4,4))
            gc.show_colorscale(vmin = min,
                               vmax = med+15*mad,
                               interpolation = 'none',
                               cmap = 'Greys')
            del d

            gc.recenter(row_obj1[obj1.cols.ra.label],
                        row_obj1[obj1.cols.dec.label],
                        radius=set_image_radii(obj2, unit=u.deg, stretch=stretch))

            # Obj2 DATA
            filename_obj2 = obj2.matches[obj1.name].images[i].file_location

            dd, dd_med, dd_mad, _  = collect_statistics(filename_obj2)
            bmaj, bmin, bpa = get_beam_info(filename_obj2)
            gc.show_contour(dd,
                            colors=['#7fcdbb', '#edf8b1', '#c7e9b4', '#ffffd9'],
                            labels=[r'$\mu+3\sigma$', r'$\mu+5\sigma$', r'$\mu+10\sigma$', r'$\mu+15\sigma$'],
                            alpha=0.75,
                            levels=[dd_med + 3*dd_mad,
                                    dd_med + 5*dd_mad,
                                    dd_med + 10*dd_mad,
                                    dd_med + 15*dd_mad]
                            )
            del dd

            gc.add_beam(major=bmaj, minor=bmin, angle=bpa)
            gc.beam.set_edgecolor('black')
            gc.beam.set_facecolor('none')
            gc.beam.set_hatch('+')

            gc.show_ellipses(row_obj1[obj1.cols.ra.label],
                             row_obj1[obj1.cols.dec.label],
                             width=(6 * u.arcsec).to(u.deg).value,
                             height=(6 * u.arcsec).to(u.deg).value,
                             angle=0,
                             edgecolor=pink_rgb,
                             linestyle='-',
                             facecolor='none',
                             zorder=100)

            # gc.set_title('LoTSS df id: %d' %j)
            gc.set_title('{}{}'.format(
                    row_obj1[obj1.cols.tns.name.label],
                    ':Type {}'.format(row_obj1[obj1.cols.obj_type.label]) #if obj1.prop('obj_type') is not None else ''
                ))

            gc.set_theme('publication')

            outpath = check_folder_exists_or_create('images/matches/{}/{}'.format(filter, stretch) if extra_path is None else 'images/matches/{}/{}/{}'.format(extra_path, filter, stretch),
                                                    return_folder=True)

            outfile = "{}/{}_{}.pdf".format(outpath,
                                            row_obj1[obj1.cols.tns.name.label].replace(' ', ''),
                                            stretch)
            print ('Saving to {}'.format(outfile))
            gc.save(outfile, dpi=250)

            gc.close()
        except KeyError:
            pass

def plot_lotss(obj1, obj2, mask1=None, mask2=None, d2d_obj1=None, d2d_obj2=None, extra_path=None):
    # PanSTARRS r filter fits most of CLU's 4 filters (H-alpha 656, 663, 672, 681)
    matched_obj1, matched_obj2 = get_matched(obj1, obj2, mask1, mask2)

    for (i, row_obj1), (j, row_obj2) in tqdm(zip(matched_obj1.iterrows(), matched_obj2.iterrows()),
                                             total=len(matched_obj1)):
        # Obj2 DATA
        filename_obj2 = obj2.matches[obj1.name].images[i].file_location
        d, med, mad, _  = collect_statistics(filename_obj2)
        bmaj, bmin, bpa = get_beam_info(filename_obj2)
        gc = aplpy.FITSFigure(d, figsize=(4,4))
        gc.show_colorscale(vmin = med-15*mad,
                           vmax = med+15*mad,
                           interpolation = 'none',
                           cmap = 'RdBu_r')

        gc.show_contour(d,
                        colors=['#7fcdbb', '#edf8b1', '#c7e9b4', '#ffffd9'],
                        labels=[r'$\mu+3\sigma$', r'$\mu+5\sigma$', r'$\mu+10\sigma$', r'$\mu+15\sigma$'],
                        alpha=0.75,
                        levels=[med + 3*mad,
                                med + 5*mad,
                                med + 10*mad,
                                med + 15*mad])
        del d

        gc.recenter(row_obj1[obj1.cols.ra.label],
                    row_obj1[obj1.cols.dec.label],
                    radius=set_image_radii(obj1, obj2, i, j, unit=u.deg))

        gc.show_ellipses(row_obj1[obj1.cols.ra.label],
                         row_obj1[obj1.cols.dec.label],
                         width=obj1.semi_major(i, to_unit=u.deg) * 2,
                         height=obj1.semi_minor(i, to_unit=u.deg) * 2,
                         angle=-row_obj1[obj1.cols.pa.label],
                         edgecolor=pink_rgb,
                         linestyle='-',
                         facecolor='none')

        gc.add_beam(major=bmaj, minor=bmin, angle=bpa)
        gc.beam.set_edgecolor('black')
        gc.beam.set_facecolor('none')
        gc.beam.set_hatch('+')

        gc.show_markers(row_obj2[obj2.cols.ra.label],
                        row_obj2[obj2.cols.dec.label],
                        marker='+',
                        alpha=0.75,
                        edgecolor=blue_rgb,
                        facecolor=blue_rgb)

        gc.show_ellipses(row_obj2[obj2.cols.ra.label],
                         row_obj2[obj2.cols.dec.label],
                         width=(row_obj2[obj2.cols.ra_err.label] * obj2.cols.ra_err.unit).to(u.deg) * 2,
                         height=(row_obj2[obj2.cols.dec_err.label] * obj2.cols.dec_err.unit).to(u.deg) * 2,
                         angle=0,
                         edgecolor=blue_rgb,
                         linestyle='--',
                         facecolor='none')

        gc.show_ellipses(row_obj2[obj2.cols.ra.label],
                         row_obj2[obj2.cols.dec.label],
                         width=obj2.semi_major(j, to_unit=u.deg) * 2,
                         height=obj2.semi_minor(j, to_unit=u.deg) * 2,
                         angle=-row_obj2[obj2.cols.pa.label],
                         edgecolor=blue_rgb,
                         linestyle='-',
                         facecolor='none')


        # Find all other galaxies in image frame
        image_radii = set_image_radii(obj1, obj2, i, j, u.arcsec)
        # CASE clu
        if d2d_obj1 is None:
            local_d2d = obj1.as_SkyCoord(i).separation(obj1.as_SkyCoord())
        else:
            local_d2d = d2d_obj1[i]['d2d']
        test_mask = (local_d2d <= image_radii * u.arcsec)

        if len(test_mask) > 1 or d2d_obj1 is None:
            for k, match in obj1.df.iloc[
                test_mask if d2d_obj1 is None else d2d_obj1[i]['filtered_idx'][test_mask]
            ].iterrows():
                if k != i:
                    gc.show_ellipses(match[obj1.cols.ra.label],
                                     match[obj1.cols.dec.label],
                                     width=obj1.semi_major(k, to_unit=u.deg) * 2,
                                     height=obj1.semi_minor(k, to_unit=u.deg) * 2,
                                     angle=-match[obj1.cols.pa.label] if not pd.isnull(match[obj1.cols.pa.label]) else 0,
                                     edgecolor='darkgrey',
                                     alpha=0.65,
                                     linestyle='dotted',
                                     facecolor='none')
        # CASE lotss
        if d2d_obj2 is None:
            local_d2d = obj1.as_SkyCoord(i).separation(obj2.as_SkyCoord())
        else:
            local_d2d = d2d_obj2[i]['d2d']
        test_mask = (local_d2d <= image_radii * u.arcsec)
        if len(test_mask) > 1 or d2d_obj2 is None:
            for k, match in obj2.df.iloc[
                test_mask if d2d_obj2 is None else d2d_obj2[i]['filtered_idx'][test_mask]
            ].iterrows():
                if k != j:
                    gc.show_ellipses(match[obj2.cols.ra.label],
                                     match[obj2.cols.dec.label],
                                     width=obj2.semi_major(k, to_unit=u.deg) * 2,
                                     height=obj2.semi_minor(k, to_unit=u.deg) * 2,
                                     angle=-match[obj2.cols.pa.label],
                                     edgecolor='grey',
                                     alpha=0.7,
                                     facecolor='none')

        gc.set_title('LoTSS df id: %d' % j)

        gc.set_theme('publication')
        outpath = check_folder_exists_or_create('images/matches/lotss' if extra_path is None else 'images/matches/{}/lotss'.format(extra_path),
                                                return_folder=True)
        outfile = "%s/%s_%s.pdf" % (outpath, j, get_name(row_obj1[obj1.cols.ra.label],
                                                       row_obj1[obj1.cols.dec.label]))
        print ('Saving to {}'.format(outfile))
        gc.save(outfile, dpi=250)
        gc.close()

def plot_lotss_asiago(obj1, obj2, mask1=None, mask2=None, d2d_obj2=None, extra_path=None):
    # PanSTARRS r filter fits most of CLU's 4 filters (H-alpha 656, 663, 672, 681)
    matched_obj1, matched_obj2 = get_matched(obj1, obj2, mask1, mask2)

    for (i, row_obj1), (j, row_obj2) in tqdm(zip(matched_obj1.iterrows(), matched_obj2.iterrows()),
                                             total=len(matched_obj1)):
        # Obj2 DATA
        filename_obj2 = obj2.matches[obj1.name].images[i].file_location

        d, med, mad, _  = collect_statistics(filename_obj2)
        bmaj, bmin, bpa = get_beam_info(filename_obj2)
        gc = aplpy.FITSFigure(d, figsize=(4,4))
        gc.show_colorscale(vmin = med-15*mad,
                           vmax = med+15*mad,
                           interpolation = 'none',
                           cmap = 'RdBu_r')

        gc.add_colorbar()
        gc.colorbar.set_axis_label_text('Flux (Jy/beam)')

        gc.show_contour(d,
                        colors=['#7fcdbb', '#edf8b1', '#c7e9b4', '#ffffd9'],
                        labels=[r'$\mu+3\sigma$', r'$\mu+5\sigma$', r'$\mu+10\sigma$', r'$\mu+15\sigma$'],
                        alpha=0.75,
                        levels=[med + 3*mad,
                                med + 5*mad,
                                med + 10*mad,
                                med + 15*mad])
        del d

        gc.add_beam(major=bmaj, minor=bmin, angle=bpa)
        gc.beam.set_edgecolor('black')
        gc.beam.set_facecolor('none')
        gc.beam.set_hatch('+')

        gc.show_markers(row_obj1[obj1.cols.ra.label],
                        row_obj1[obj1.cols.dec.label],
                        marker='*',
                        facecolor='#edf8b1',
                        edgecolor='#edf8b1',
                        alpha=0.75,
                        zorder=100)

        gc.recenter(row_obj1[obj1.cols.ra.label],
                    row_obj1[obj1.cols.dec.label],
                    radius=set_image_radii(obj1, obj2, i, j, unit=u.deg))

        gc.show_markers(row_obj2[obj2.cols.ra.label],
                        row_obj2[obj2.cols.dec.label],
                        marker='+',
                        alpha=1,
                        edgecolor=blue_rgb,
                        facecolor=blue_rgb)
        #
        # gc.show_ellipses(row_obj2[obj2.cols.ra.label],
        #                  row_obj2[obj2.cols.dec.label],
        #                  width=(row_obj2[obj2.cols.ra_err.label] * obj2.cols.ra_err.unit).to(u.deg) * 2,
        #                  height=(row_obj2[obj2.cols.dec_err.label] * obj2.cols.dec_err.unit).to(u.deg) * 2,
        #                  angle=0,
        #                  edgecolor=blue_rgb,
        #                  linestyle='--',
        #                  facecolor='none')
        #
        gc.show_ellipses(row_obj2[obj2.cols.ra.label],
                         row_obj2[obj2.cols.dec.label],
                         width=obj2.semi_major(j, to_unit=u.deg) * 2,
                         height=obj2.semi_minor(j, to_unit=u.deg) * 2,
                         angle=-row_obj2[obj2.cols.pa.label],
                         edgecolor=blue_rgb,
                         linestyle='-',
                         facecolor='none')
        #
        #
        # # Find all other galaxies in image frame
        # image_radii = set_image_radii(obj1, obj2, i, j, u.arcsec)
        #
        # # CASE lotss
        # if d2d_obj2 is None:
        #     local_d2d = obj1.as_SkyCoord(i).separation(obj2.as_SkyCoord())
        # else:
        #     local_d2d = d2d_obj2[i]['d2d']
        # test_mask = (local_d2d <= image_radii * u.arcsec)
        # if len(test_mask) > 1 or d2d_obj2 is None:
        #     for k, match in obj2.df.iloc[
        #         test_mask if d2d_obj2 is None else d2d_obj2[i]['filtered_idx'][test_mask]
        #     ].iterrows():
        #         if k != j:
        #             gc.show_ellipses(match[obj2.cols.ra.label],
        #                              match[obj2.cols.dec.label],
        #                              width=obj2.semi_major(k, to_unit=u.deg) * 2,
        #                              height=obj2.semi_minor(k, to_unit=u.deg) * 2,
        #                              angle=-match[obj2.cols.pa.label],
        #                              edgecolor='grey',
        #                              alpha=0.7,
        #                              facecolor='none')

        # gc.set_title('LoTSS df id: {}, {}:Type {}'.format(j, row_obj1['Name'], row_obj1['Obj. Type']))
        gc.set_title('LoTSS df id: {}, {}{}'.format(
            j,
            row_obj1[obj1.cols.name.label],
            ':Type {}'.format(row_obj1[obj1.cols.obj_type.label]) if obj1.prop('obj_type') is not None else ''
        ))
        gc.set_theme('publication')
        outpath = check_folder_exists_or_create('images/matches/lotss' if extra_path is None else 'images/matches/{}/lotss'.format(extra_path),
                                                return_folder=True)
        outfile = "%s/%s_%s.pdf" % (outpath, j, get_name(row_obj1[obj1.cols.ra.label],
                                                       row_obj1[obj1.cols.dec.label]))
        print ('Saving to {}'.format(outfile))
        gc.save(outfile, dpi=250)
        gc.close()


def plot_source_on_lotss(obj1, obj2, mask1=None, core_collapse_only=False, stretch=4, extra_path=None, verbose=True):
    # PanSTARRS r filter fits most of CLU's 4 filters (H-alpha 656, 663, 672, 681)

    matched_obj1 = obj1.df.iloc[mask1] if mask1 is not None else obj1.df
    if core_collapse_only:
        # This is a temporary use case that should be done some other, generic, way
        matched_obj1 = matched_obj1.loc[
            (matched_obj1[obj1.cols.tns.obj_type.label] != 'SN Ia') &
            (matched_obj1['Redshift'] <= 0.005) &
            (~matched_obj1['Redshift'].isna())
        ]

    for i, row_obj1 in tqdm(matched_obj1.iterrows(), total=len(matched_obj1)):
        try:
            # Obj2 DATA
            filename_obj2 = obj2.matches[obj1.name].images[i].file_location

            d, med, mad, min  = collect_statistics(filename_obj2)
            bmaj, bmin, bpa = get_beam_info(filename_obj2)
            gc = aplpy.FITSFigure(d, figsize=(4,4))
            gc.show_colorscale(vmin = min,
                               vmax = med+15*mad,
                               interpolation = 'none',
                               cmap = 'RdBu_r')

            gc.add_colorbar()
            gc.colorbar.set_axis_label_text('Flux (Jy/beam)')

            gc.show_contour(d,
                            colors=['#7fcdbb', '#edf8b1', '#c7e9b4', '#ffffd9'],
                            labels=[r'$\mu+3\sigma$', r'$\mu+5\sigma$', r'$\mu+10\sigma$', r'$\mu+15\sigma$'],
                            alpha=0.75,
                            levels=[med + 3*mad,
                                    med + 5*mad,
                                    med + 10*mad,
                                    med + 15*mad])
            del d

            gc.add_beam(major=bmaj, minor=bmin, angle=bpa)
            gc.beam.set_edgecolor('black')
            gc.beam.set_facecolor('none')
            gc.beam.set_hatch('+')

            # gc.show_markers(row_obj1[obj1.cols.ra.label],
            #                 row_obj1[obj1.cols.dec.label],
            #                 marker='*',
            #                 facecolor='#edf8b1',
            #                 edgecolor='#edf8b1',
            #                 alpha=0.75,
            #                 zorder=100)

            gc.show_ellipses(row_obj1[obj1.cols.ra.label],
                             row_obj1[obj1.cols.dec.label],
                             width=(5 * u.arcsec).to(u.deg).value,
                             height=(5 * u.arcsec).to(u.deg).value,
                             angle=0,
                             # edgecolor='#edf8b1',
                             edgecolor=pink_rgb,
                             linestyle='-',
                             facecolor='none',
                             zorder=100)

            gc.recenter(row_obj1[obj1.cols.ra.label],
                        row_obj1[obj1.cols.dec.label],
                        radius=set_image_radii(obj2, unit=u.deg, stretch=stretch))

            gc.set_title('{}{}'.format(
                row_obj1[obj1.cols.tns.name.label],
                ':Type {}'.format(row_obj1[obj1.cols.obj_type.label]) #if obj1.prop('obj_type') is not None else ''
            ))
            gc.set_theme('publication')
            outpath = check_folder_exists_or_create('images/matches/lotss' if extra_path is None else 'images/matches/{}/lotss'.format(extra_path),
                                                    return_folder=True)
            outfile = "{}/{}_{}.pdf".format(outpath,
                                            row_obj1[obj1.cols.tns.name.label].replace(' ', ''),
                                            stretch)
            if verbose:
                print ('Saving to {}'.format(outfile))
            gc.save(outfile, dpi=250)
            gc.close()
        except FileNotFoundError:
            # print ('Missing: {}'.format(row_obj1[obj1.cols.tns.name.label]))
            pass
        # except:
        #     print ('Missing/outside field: {}'.format(row_obj1[obj1.cols.tns.name.label]))
        #     # print (row_obj1[[obj1.cols.ra.label, obj1.cols.dec.label]])
