import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import colors

from .utils import set_size

from ..matchmaker.utils import (
    check_folder_exists_or_create,
    check_strictly_positive,
    sort_dict_by_value,
    intersection
)
from ..matchmaker.measures import (
    luminosity_distance,
    w_hz_to_erg_s_hz,
    erg_s_hz_to_w_hz,
)
from ..matchmaker.model.lofar import power_law, sigma as l_sfr_sigma
from ..matchmaker.dataset import sources
from ..matchmaker.model import wise


def mscatter(x, y, ax=None, m=None, **kw):
    # Found on stack overflow. TODO: find ref.
    import matplotlib.markers as mmarkers
    if not ax: ax=plt.gca()
    sc = ax.scatter(x,y,**kw)
    if (m is not None) and (len(m)==len(x)):
        paths = []
        for marker in m:
            if isinstance(marker, mmarkers.MarkerStyle):
                marker_obj = marker
            else:
                marker_obj = mmarkers.MarkerStyle(marker)
            path = marker_obj.get_path().transformed(
                        marker_obj.get_transform())
            paths.append(path)
        sc.set_paths(paths)
    return sc

def diagnostics(obj1, obj2,
                save=False, path='images/matches',
                chance=None, use_a=True,
                show_colorbar=True,
                show_legend=True,
                z_range=None):
    to_plot = [
        # 'plot_flux',
        'wise_color-color',
        'lum_mstar',
        'lum_sfr',
        # 'lum_ssfr',
        # 'mag_mstar',
        # 'ssfr_mag',
        # 'mstar_ssfr',
        # 'lum_ssfr_mstar',
        # 'ssfr_lum_sfr',
        # 'ssfr_lum_m',
        'lum_offset'
    ]
    J, I = len(to_plot), 1
    fig, _ax = plt.subplots(J,I, figsize=(8, 4*J))
    size=15
    params = {'legend.fontsize': 'large',
              # 'figure.figsize': (20,8),
              'axes.labelsize': size,
              'axes.titlesize': size,
              'xtick.labelsize': size*0.75,
              'ytick.labelsize': size*0.75,
              'axes.titlepad': 25}
    plt.rcParams.update(params)

    # DATA PREP
    FRBs = [sources.frb121102, sources.frb190520b]
    matched_obj1 = obj1.df.iloc[obj1.matches[obj2.name].mask]
    matched_obj2 = obj2.df.iloc[obj2.matches[obj1.name].filtered_idx]
    matched = obj2.matches[obj1.name]

    lum_w_hz = obj2.luminosity_distance(obj2.matches[obj1.name].filtered_idx)

    # FORMATTING UTILS
    size = 15
    c = matched_obj1[obj1.cols.z.label]
    if z_range is None:
        norm = colors.LogNorm()
    else:
        norm = colors.LogNorm(vmin=check_strictly_positive(z_range[0]),
                              vmax=check_strictly_positive(z_range[1]))

    # Check which sources are likely AGNs
    matched_w1w2 = matched_obj1['w1mpro'] - matched_obj1['w2mpro']
    matched_w2w3 = matched_obj1['w2mpro'] - matched_obj1['w3mpro']

    # for i, idx in enumerate(matched.filtered_idx):
    #     if wise.is_AGN(matched_w1w2.values[i], matched_w2w3.values[i]):
    #         print (matched_obj1['sfr_fuv'].values[i],
    #                matched_w1w2.values[i],
    #                matched_w2w3.values[i],
    #                wise.is_AGN(matched_w1w2.values[i], matched_w2w3.values[i]))

    # markers = ['D' if i in matched.filtered_compact_idx else 'o' \
    #            for i, idx in enumerate(matched.filtered_idx)]

    markers = ['s' if i in matched.with_spectra else 'D' \
                   if i in matched.filtered_compact_idx else 'o' \
               for i, idx in enumerate(matched.filtered_idx)]

    edge_colors = ['red' if wise.is_AGN(matched_w1w2.values[i], matched_w2w3.values[i]) else \
                       'blue' if i in matched.with_spectra else 'black' \
        if i in matched.filtered_compact_idx else 'none' for i, idx in enumerate(matched.filtered_idx)]

    # TODO: use zorder to put diamonds in front of the other markers
    zorder = np.argsort([3 if wise.is_AGN(matched_w1w2.values[i], matched_w2w3.values[i]) else 2 \
        if i in matched.filtered_compact_idx else 1 for i, idx in enumerate(matched.filtered_idx)])

    i, j = 0, 0

    # for p, ax in zip(to_plot, grid):
    for p in to_plot:
        if J > 1 and I > 1:
            ax = _ax[j, i]
        else:
            ax = _ax[i]

        if p == 'plot_flux':
            ax.set_title('Expected chance association: %d out of %d' % (chance.expected, chance.total))
            scat = mscatter(matched_obj1['mstar'],
                            matched_obj2['Peak_flux'],
                            s=size,
                            c=c,
                            m=markers,
                            norm=norm,
                            ax=ax)
            ax.axhline(0.8)
            ax.axvline(sources.m_lmc, linestyle='--')
            ax.set_xlabel(r'Total stellar mass ($M_\odot$)')
            ax.set_ylabel('$S_{150~MHz}~(mJy)$')
            ax.set_xscale('log')

        if p == 'lum_mstar':
            scat = mscatter(matched_obj1['mstar'],
                            lum_w_hz,
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            sec_ax = ax.secondary_yaxis('right', functions=(w_hz_to_erg_s_hz, erg_s_hz_to_w_hz))
            sec_ax.set_ylabel('$L_{150~MHz}~(erg~s^{-1}~Hz^{-1})$')

            for frb in FRBs:
                ax.scatter(frb['Msun'],
                           frb['lum'],
                           marker='*',
                           c=frb['color'],
                           s=frb['s'],
                           label=frb['label'])
            # ax.axhline(sources.l_ptf10hgi, c='red', label='PTF10hgi (6 GHz)')
            ax.axvline(sources.m_lmc, linestyle='--')
            ax.set_xlabel(r'Total stellar mass ($M_\odot$)')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})$')
            ax.set_xscale('log')
            ax.set_yscale('log')
            if show_legend:
                ax.legend()

        if p == 'lum_ssfr':
            scat = mscatter(matched_obj1['sfr_fuv']/matched_obj1['mstar'],
                            lum_w_hz,
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            sec_ax = ax.secondary_yaxis('right', functions=(w_hz_to_erg_s_hz, erg_s_hz_to_w_hz))
            sec_ax.set_ylabel('$L_{150~MHz}~(erg~s^{-1}~Hz^{-1})$')

            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(r'sSFR')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})$')

        if p == 'lum_sfr':
            scat = mscatter(matched_obj1['sfr_fuv'],
                            lum_w_hz,
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            sec_ax = ax.secondary_yaxis('right', functions=(w_hz_to_erg_s_hz, erg_s_hz_to_w_hz))
            sec_ax.set_ylabel('$L_{150~MHz}~(erg~s^{-1}~Hz^{-1})$')

            for frb in FRBs:
                ax.scatter(frb['sfr'],
                           frb['lum'],
                           marker='*',
                           c=frb['color'],
                           s=frb['s'],
                           label=frb['label'])
            # ax.axhline(sources.l_ptf10hgi, c='red', label='PTF10hgi')

            ax.set_xscale('log')
            ax.set_yscale('log')

            # Add LOFAR's L-SFR relation (Gürkan et al.)
            xlim = ax.get_xlim()
            sfrs = np.linspace(xlim[0], xlim[1], 2)
            ax.plot(sfrs, power_law(sfrs), c='grey')
            ax.plot(sfrs, power_law(sfrs) * l_sfr_sigma, c='grey', linestyle='--')
            ax.plot(sfrs, power_law(sfrs) * (1.5 * l_sfr_sigma), c='grey', linestyle=':')
            ax.plot(sfrs, power_law(sfrs) / l_sfr_sigma, c='grey', linestyle='--')
            ax.set_xlim(xlim)

            ax.set_xlabel(r'SFR ($M_\odot$/yr)')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})$')


        if p == 'wise_color-color':
            scat = mscatter(matched_obj1['w2mpro'] - matched_obj1['w3mpro'],
                            matched_obj1['w1mpro'] - matched_obj1['w2mpro'],
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            # AGN box corners (from Jarret+2011)
            box = wise.AGN_color_color_box
            maxed = 1.7 if matched_w1w2.max() >= 1.7 else ax.get_ylim()[1]
            if maxed != 1.7:
                ylim = ax.get_ylim()
            ax.plot(
                [box[2][0], box[0][0]], [box[2][1], maxed],
                linestyle='--',
                color='black'
            )
            if maxed == 1.7:
                ax.plot(
                    [box[0][0], box[1][0]], [box[0][1], maxed],
                    linestyle='--',
                    color='black'
                )
            ax.plot(
                [box[2][0], box[3][0]], [box[2][1], box[3][1]],
                linestyle='--',
                color='black'
            )
            ax.plot(
                [box[3][0], box[1][0]], [box[3][1], box[1][1]],
                linestyle='--',
                color='black'
            )
            if maxed != 1.7:
                ax.set_ylim(ylim)
            ax.set_xlabel('W2-W3 (mag)')
            ax.set_ylabel('W1-W2 (mag)')

        if p == 'mag_mstar':
            scat = ax.scatter(matched_obj1['mstar'],
                              matched_obj1['modelmag_g'] - matched_obj1['modelmag_r'],
                              s=size,
                              label='SDSS',
                              c=c,
                              norm=norm)
            ax.axvline(sources.m_lmc, linestyle='--')
            ax.set_xlabel(r'Total stellar mass ($M_\odot$)')
            ax.set_ylabel('SDSS g-r (mag)')
            ax.set_xscale('log')

        if p == 'ssfr_mag':
            scat = ax.scatter(matched_obj1['modelmag_g'] - matched_obj1['modelmag_r'],
                              np.log10(matched_obj1['sfr_fuv'] / matched_obj1['mstar']),
                              s=size,
                              c=c,
                              norm=norm)
            ax.set_ylabel('log sSFR')
            ax.set_xlabel('SDSS g-r (mag)')
            fig.colorbar(scat, label='$\logz$', ax=ax)

        if p == 'mstar_ssfr':
            scat = ax.scatter(matched_obj1['mstar'],
                              np.log10(matched_obj1['sfr_fuv'] / matched_obj1['mstar']),
                              s=size,
                              c=c,
                              norm=norm)
            ax.set_ylabel('log sSFR')
            ax.set_xlabel(r'Sellar mass ($M_\odot$)')

        if p == 'lum_ssfr_mstar':
            scat = mscatter((matched_obj1['sfr_fuv']/matched_obj1['mstar'])/matched_obj1['mstar'],
                            lum_w_hz,
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            ax.set_xlabel(r'$\log_{10} (sSFR/M_\odot)$')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})$')
            ax.set_xscale('log')
            ax.set_yscale('log')

        if p == 'ssfr_lum_sfr':
            scat = mscatter((matched_obj1['sfr_fuv']/matched_obj1['mstar']),
                            lum_w_hz / matched_obj1['sfr_fuv'],
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            ax.set_xlabel(r'$\log_{10} (sSFR/M_\odot)$')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})/SFR$')
            ax.set_xscale('log')
            ax.set_yscale('log')
            fig.colorbar(scat, label=r'$z$')

        if p == 'ssfr_lum_m':
            scat = mscatter((matched_obj1['sfr_fuv']),
                            lum_w_hz / (matched_obj1['mstar']/(1e10 * u.solMass).value)**0.44,
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            for frb in FRBs:
                ax.scatter(frb['sfr'],
                           frb['lum'].to(u.W * u.Hz**-1).value / (frb['Msun']/(1e10 * u.solMass).value)**0.44,
                           marker='*',
                           c=frb['color'],
                           s=frb['s'],
                           label=frb['label'])
            # ax.axhline(sources.l_ptf10hgi.to(u.W * u.Hz**-1).value, c='red', label='PTF10hgi')

            ax.set_xlabel(r'$SFR (M_\odot~yr^{-1})$')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})/(M/10^{10}M_{\odot})^{0.44}$')
            ax.set_xscale('log')
            ax.set_yscale('log')

        if p == 'lum_offset':
            match = obj2.matches[obj1.name]

            scat = mscatter(match.sep3d.to(u.kpc), # Need to compute offset (Eftekhari quotes it in parsec)
                            lum_w_hz,
                            s=size,
                            c=c,
                            m=markers,
                            edgecolors=edge_colors,
                            norm=norm,
                            ax=ax)

            sec_ax = ax.secondary_yaxis('right', functions=(w_hz_to_erg_s_hz, erg_s_hz_to_w_hz))
            sec_ax.set_ylabel('$L_{150~MHz}~(erg~s^{-1}~Hz^{-1})$')

            ax.set_ylim(4e17, 4e23)
            ax.set_xlabel(r'Offset (kpc)')
            ax.set_ylabel('$L_{150~MHz}~(W~Hz^{-1})$')
            ax.set_xscale('log')
            ax.set_yscale('log')

        if J > 1 and I > 1:
            if i == I-1:
                i = 0
                j += 1
            else:
                i += 1
        else:
            i += 1


    fig.subplots_adjust(top=0.99,
                        left=0.1, right=0.77,
                        bottom=0.05,
                        wspace=0.02, hspace=0.2)
    if show_colorbar:
        cb_ax = fig.add_axes([
            0.87,  # left
            0.1,  # bottom
            0.02, # width
            0.8   # height
        ])
        cbar = fig.colorbar(scat,
                            cax=cb_ax,
                            label=r'$z$',
                            shrink=0.6)

    if save:
        check_folder_exists_or_create(path)
        filename = '{}/diagnostic_{}'.format(path, 'a' if use_a else 'no_a')
        print('Saving diagnostics to: {}'.format(filename))
        plt.savefig('{}.pdf'.format(filename))
        # plt.savefig('{}.png'.format(filename), dpi=250)

def spectrum(datasets:dict, reference_dataset:str, source_dataset:str, reference_base_mask=None,
             freq_unit_plot=u.Hz, flux_unit_plot=u.uJy,
             save=False, path='images/matches',
             individual=False):
    # Too many ifs for my liking... this could likely be simplified.

    ref_masks = {}
    lens = {}
    common = None
    try:
        datasets[reference_dataset].matches[source_dataset].radio_spectra
    except:
        print ('init radio_spectra dict')
        datasets[reference_dataset].matches[source_dataset].radio_spectra = {}

    for name in datasets.keys():
        if name != reference_dataset:
            if common is None:
                common = datasets[reference_dataset].matches[name].filtered_idx
                if reference_base_mask is not None:
                    common = reference_base_mask[common]
            else:
                mask = datasets[reference_dataset].matches[name].filtered_idx
                if reference_base_mask is not None:
                    mask = reference_base_mask[mask]
                common = intersection(common, mask)

    if len(common) > 0:
        freq_flux = []
        for i in common:
            freqs = []
            fluxes = []
            fluxes_err = []
            for name in datasets.keys():
                obj = datasets[name]
                freqs.append(obj.central_frequency.to(freq_unit_plot).value)
                if name == reference_dataset:
                    fluxes.append(
                        (obj.df.iloc[i][obj.cols.measure.total_flux.label] * obj.cols.measure.total_flux.unit).to(flux_unit_plot).value
                    )
                    try:
                        fluxes_err.append(
                            (obj.df.iloc[i][obj.cols.measure.total_flux_err.label] * obj.cols.measure.total_flux_err.unit).to(flux_unit_plot).value
                        )
                    except AttributeError:
                        fluxes_err.append(0)
                else:
                    mask = datasets[reference_dataset].matches[name].mask_idx
                    ref_mask = datasets[reference_dataset].matches[name].filtered_idx
                    if reference_base_mask is not None:
                        ref_mask = reference_base_mask[ref_mask]
                    ref_idx = np.where(ref_mask == i)[0][0]
                    fluxes.append(
                        (obj.df.iloc[mask[ref_idx]][obj.cols.measure.total_flux.label] * obj.cols.measure.total_flux.unit).to(flux_unit_plot).value
                    )
                    try:
                        fluxes_err.append(
                            (obj.df.iloc[mask[ref_idx]][obj.cols.measure.total_flux_err.label] * obj.cols.measure.total_flux_err.unit).to(flux_unit_plot).value
                        )
                    except AttributeError:
                        fluxes_err.append(0)
            freq_flux.append([np.array(freqs), np.array(fluxes), np.array(fluxes_err), i])

        if not individual:
            fig, ax = plt.subplots(1, 1, figsize=set_size(width=900, vertical=False))
            for freq, flux, flux_err, i in freq_flux:
                _sorted = np.argsort(freq)
                ax.errorbar(freq[_sorted], flux[_sorted], yerr=flux_err, fmt='o-', label='lotss df id:{}'.format(i))
            ax.set_xlabel(r'Frequency ({})'.format(freq_unit_plot.to_string('latex_inline')))
            ax.set_ylabel(r'Total flux ({})'.format(flux_unit_plot.to_string('latex_inline')))
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.legend()

            size=15
            params = {'legend.fontsize': 'large',
                      # 'figure.figsize': (20,8),
                      'axes.labelsize': size,
                      'axes.titlesize': size,
                      'xtick.labelsize': size*0.95,
                      'ytick.labelsize': size*0.95,
                      'axes.titlepad': 25}
            plt.rcParams.update(params)

            if save:
                check_folder_exists_or_create(path)
                filename = '{}/spectra.pdf'.format(path)
                plt.savefig(filename)
        else:
            for freq, flux, flux_err, i in freq_flux:
                fig, ax = plt.subplots(1, 1, figsize=set_size(width=900, vertical=False))
                _sorted = np.argsort(freq)
                ax.errorbar(freq[_sorted], flux[_sorted], yerr=flux_err, fmt='o-', label='lotss df id:{}'.format(i))
                ax.set_xlabel(r'Frequency ({})'.format(freq_unit_plot.to_string('latex_inline')))
                ax.set_ylabel(r'Total flux ({})'.format(flux_unit_plot.to_string('latex_inline')))
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.legend()

                size=15
                params = {'legend.fontsize': 'large',
                          # 'figure.figsize': (20,8),
                          'axes.labelsize': size,
                          'axes.titlesize': size,
                          'xtick.labelsize': size*0.95,
                          'ytick.labelsize': size*0.95,
                          'axes.titlepad': 25}
                plt.rcParams.update(params)

                if save:
                    check_folder_exists_or_create(path)
                    filename = '{}/spectra_{}.pdf'.format(path, i)
                    print ('saving to {}'.format(filename))
                    if i not in datasets[reference_dataset].matches[source_dataset].radio_spectra.keys():
                        plt.savefig(filename)
                        datasets[reference_dataset].matches[source_dataset].radio_spectra[i] = {
                            'freq': freq[_sorted],
                            'flux': flux[_sorted],
                            'flux_err': flux_err[_sorted]
                        }

    return common
    # else:
    #     print ('No common match between all sets')

    # return ref_indices, ref_masks


