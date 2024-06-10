import os
import sys
import numpy as np

from tqdm import tqdm
import copy

# Add module path to system path
module_paths = ['..']
for module_path in module_paths:
    if os.path.abspath(os.path.join(module_path)) not in sys.path:
        sys.path.insert(0, module_path)
        
import astropy.units as u
from matchmaker.matchmaker.crossmatch import crossmatch
from matchmaker.matchmaker.utils import save, load


def apply_rand_shift(ra, dec, shift):
    # ra_offset = (((np.random.rand()-1)*shift) * u.arcmin).to(u.deg)
    # dec_offset = (((np.random.rand()-1)*shift) * u.arcmin).to(u.deg)
    ra_offset = (np.random.uniform(-shift, shift) * u.arcmin).to(u.deg)
    dec_offset = (np.random.uniform(-shift, shift) * u.arcmin).to(u.deg)
    
    # ra_ = (ra + np.random.normal(scale=scale))
    ra_ = (ra + ra_offset)
    ra_[np.where(ra_ > 360)[0]] = ra_[np.where(ra_ > 360)[0]] - 360
    ra_[np.where(ra_ < 0)[0]] = 360 + ra_[np.where(ra_ < 0)[0]]

    # dec_ = (dec + np.random.normal(scale=20))
    dec_ = dec + dec_offset
    dec_[np.where(dec_ > 90)[0]] = -(dec_[np.where(dec_ > 90)[0]] - 90)
    dec_[np.where(dec_ < -90)[0]] = -(dec_[np.where(dec_ < -90)[0]] + 90)
    return ra_, dec_, ra_offset, dec_offset



# Main
def run_mc(state_prefix_input, state_prefix, theta, nsigma = 3, shift=10, base_folder='states'):
    datasets = load('datasets', state_prefix=state_prefix_input)

    n_cands, n_cands_all, matches, ra_offsets, dec_offsets = [], [], [], [], []
    x1, x2, y1, y2 = [], [], [], []

    mc_datasets = {}
    ra_label, dec_label = datasets['lotss'].cols.ra.label, datasets['lotss'].cols.dec.label

    """
    I want to break this loop into chunks
    """
    for m in tqdm(range(1000)):
        mc_datasets = copy.deepcopy(datasets)
        
        # Make input smaller if possible
        del mc_datasets['clu_all']
        try:        
            del mc_datasets['sdss_matched']
            del mc_datasets['sdss_matched_multiple']
            del mc_datasets['vlass']
            del mc_datasets['first']
            del mc_datasets['nvss']
            del mc_datasets['lolss']
            del mc_datasets['spectral_indices']
            del mc_datasets['reines']
        except: 
            pass
        
        mc_datasets['lotss'].df.drop(['index_df_original', 'Peak_flux', 'E_Peak_flux', 'Maj',
                                    'E_Maj', 'Min', 'E_Min', 'DC_Maj', 'E_DC_Maj', 'DC_Min', 'E_DC_Min',
                                    'PA', 'E_PA', 'DC_PA', 'E_DC_PA', 'Isl_rms', 'S_Code', 'Mosaic_ID',
                                    'Number_Pointings', 'Masked_Fraction'], axis=1)
        mc_datasets['clu'].df.drop([
            'index_df_clu', 'cluid', 'id_other', 'a', 'b2a', 'pa', 'type_ned', 'name_galex', 'ra_galex', 'dec_galex', 'fuv', 'fuverr', 
            'nuv', 'nuverr', 'name_sdss', 'ra_sdss', 'dec_sdss', 'modelmag_u', 'modelmagerr_u', 'modelmag_g', 'modelmagerr_g', 
            'modelmag_r', 'modelmagerr_r', 'modelmag_i', 'modelmagerr_i', 'modelmag_z', 'modelmagerr_z', 'name_ps1', 'ra_ps1', 
            'dec_ps1', 'kronmag_g', 'kronmagerr_g', 'kronmag_r', 'kronmagerr_r', 'kronmag_i', 'kronmagerr_i', 'kronmag_z', 
            'kronmagerr_z', 'kronmag_y', 'kronmagerr_y', 'name_2mass', 'ra_2mass', 'dec_2mass', 'r_k20fe', 'j_m_k20fe', 
            'j_msig_k20fe', 'j_flg_k20fe', 'h_m_k20fe', 'h_msig_k20fe', 'h_flg_k20fe', 'k_m_k20fe', 'k_msig_k20fe', 'k_flg_k20fe', 
            'name_wise', 'ra_wise', 'dec_wise', 'w4mpro', 'w4sigmpro', 'w4snr', 'm21', 'm21err', 'name_cluha', 'ra_cluha', 
            'dec_cluha', 'maxcsig', 'cluhamag', 'cluhamagerr', 'sfr_ha', 'sfr_haerr', 'btc', 'btcerr', 'b_r25', 'b_r25err', 
            'magb', 'magberr', 'lum_b', 'lum_berr', 'source', 'btc_source', 'size_source', 'dm_source', 'z_source', 'flags'        
        ], axis=1)
        

        # Start shifted cross-matching process
        ra, dec = mc_datasets['lotss'].df[ra_label], mc_datasets['lotss'].df[dec_label]
        mc_datasets['lotss'].df[ra_label], mc_datasets['lotss'].df[dec_label], ra_offset, dec_offset = apply_rand_shift(ra, dec, shift)

        crossmatch(mc_datasets['clu'],
                mc_datasets['lotss'],
                sep=theta,
                unit_sep=u.arcsec,
                ellipse_filter=False)

        mc_datasets['lotss'].set_distance(mc_datasets['clu'])

        match = mc_datasets['lotss'].matches['clu']
        mc_datasets['lotss'].select_lum_sfr_outliers(mc_datasets['clu'], 
                                                    mask1=match.mask_idx, 
                                                    mask2=match.filtered_idx,
                                                    n_sigma=nsigma, 
                                                    verbose = False)

        # mask1 = match.mask_idx[match.filtered_l_sfr_idx]
        mask2 = match.filtered_idx[match.filtered_l_sfr_idx]

        n_cands_all.append(len(copy.deepcopy(match.filtered_idx)))
        n_cands.append(len(copy.deepcopy(mc_datasets['lotss'].df.iloc[mask2])))
        matches.append(match)
        ra_offsets.append(copy.deepcopy(ra_offset))
        dec_offsets.append(copy.deepcopy(dec_offset))
        
        #L-SFR
        #all
        x1.extend(copy.deepcopy(mc_datasets['clu'].df.iloc[match.mask_idx]['sfr_fuv'].values))
        y1.extend(copy.deepcopy(mc_datasets['lotss'].luminosity_distance(mask=match.filtered_idx)))
        #outliers
        x2.extend(copy.deepcopy(mc_datasets['clu'].df.iloc[match.mask_idx[match.filtered_l_sfr_idx]]['sfr_fuv'].values))
        y2.extend(copy.deepcopy(mc_datasets['lotss'].luminosity_distance(mask=match.filtered_idx[match.filtered_l_sfr_idx])))
        
        match.monte_carlo = {}
        match.monte_carlo['ra_offset'] = copy.deepcopy(ra_offset)
        match.monte_carlo['ra_offset'] = copy.deepcopy(dec_offset)
        
        del mc_datasets['clu']
        del mc_datasets['lotss'].df

        save('datasets_{}_{}'.format(m, theta), mc_datasets, state_prefix='mc', folder=f'{base_folder}/monte_carlo')
        del mc_datasets

    # After all realisations are processed, save results
    ncands = np.array(n_cands)
    ncands_all = np.array(n_cands_all)

    # Saves information for 
    save('ncands_{}'.format(theta), [ncands, ncands_all, matches, ra_offsets, dec_offsets], state_prefix=state_prefix, folder=f'{base_folder}')
    save('x1x2y1y2_{}'.format(theta), [np.array(x1), np.array(x2), np.array(y1), np.array(y2)], state_prefix='mc', folder='states/monte_carlo/mass_only')


"""
# Example of master script:
datasets = {}

# Input data 
_theta = 6
state_prefix_input = 'paper_mass_flux_const_{}_arcsec_single_only'.format(_theta)

# Run the following script on various values of theta (e.g. 1 to 10 arcsec)
for theta in np.arange(1, 11, 1):
    # Processed data
    state_prefix = 'monte_carlo_mass_only_1000_shift_{}_arcmin'.format(theta)

    # Compute 1000 realisations using a cross-matching limit of `theta` arcsec
    run_mc(state_prefix_input, state_prefix, theta, base_folder='states')
"""
