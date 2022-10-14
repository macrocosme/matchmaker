from ..utils import set_size
import matplotlib.pyplot as plt
from astropy.io import fits
from specutils import Spectrum1D
import wget

from ...matchmaker.utils import check_folder_exists_or_create

def plot_sdss_spectrum(spectrum, idx, comment=None, save=False, path='spectra/matches', use_a=False, fig_width=900, vertical=False):
    # Read spectral data
    spec = Spectrum1D.read(spectrum[0], format="SDSS-III/IV spec")

    # Plot and set axes labels
    fig, ax = plt.subplots(1,1, figsize=set_size(width=fig_width, vertical=vertical))
    ax.plot(spec.spectral_axis, spec.flux)
    ax.set_xlabel(r'Wavelength ({})'.format(spec.spectral_axis.unit.to_string('latex_inline')))
    ax.set_ylabel(r'Flux ({})'.format(spec.flux.unit.to_string('latex_inline')))
    if comment is not None:
        ax.set_title(comment)
    # ax.set_title('idx' if comment is None else '{} {}'.format(idx, comment))

    if save:
        check_folder_exists_or_create(path)
        filename = '{}/{}_{}'.format(path, idx,'a' if use_a else 'no_a')
        print('Saving spectrum to: {}'.format(filename))
        plt.savefig('{}.pdf'.format(filename))

def plot_lamost_spectrum(file_url, idx, comment=None, save=False, path='spectra/matches', use_a=False, fig_width=900, vertical=False):
    d = wget.download(file_url, 'tmp/')
    hdu = fits.open(d)
    wavelength = [10**hdu[0].header['COEFF0'] + i * 10**hdu[0].header['COEFF1'] for i in range(hdu[0].header['NAXIS1'])]
    fig, ax = plt.subplots(1, 1, figsize=set_size(width=fig_width, vertical=vertical))
    ax.plot(wavelength, hdu[0].data[0])
    ax.set_xlabel(r'Wavelength ($\AA$)')
    ax.set_ylabel(r'Flux ')

    if save:
        check_folder_exists_or_create(path)
        filename = '{}/{}_lamost_{}'.format(path, idx,'a' if use_a else 'no_a')
        print('Saving spectrum to: {}'.format(filename))
        plt.savefig('{}.pdf'.format(filename))

