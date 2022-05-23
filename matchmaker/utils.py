import os
import pickle
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

sort_dict_by_value = lambda x, reverse=False: dict(sorted(x.items(),
                                                          key=lambda item: item[1],
                                                          reverse=reverse))
share_items = lambda a, b: not set(a).isdisjoint(b)
intersection = lambda a, b: np.sort(list(set(a).intersection(b)))

# Functions
def check_underscore(string):
    if string != '':
        if string[-1] != '_':
            string += '_'
    return string

def check_slash(string):
    if string != '':
        if string[-1] != '/':
            string += '/'
    return string

def check_folder_exists_or_create(folder, return_folder=True):
    if not os.path.exists(folder):
        os.makedirs(folder)
    if return_folder:
        return folder

def check_file_exists(file):
    return os.path.isfile(file)

def remove_file_from_folder(file_location):
    os.remove(file_location)

def remove_file_if_exists(filename, return_filename=False):
    if os.path.exists(filename):
        os.remove(filename)
    return filename

def save(variable:str, data, protocol=pickle.HIGHEST_PROTOCOL, state_prefix='', folder='states/'):
    check_folder_exists_or_create(folder)

    if state_prefix != '':
        with open(check_slash(folder) + check_underscore(state_prefix) + variable + '.pickle', 'wb') as f:
            pickle.dump(data, f, protocol)
    else:
        with open(check_slash(folder) + variable + '.pickle', 'wb') as f:
            pickle.dump(data, f, protocol)

def load(variable, state_prefix='', folder='states/'):
    if state_prefix != '':
        if os.path.exists(folder + check_underscore(state_prefix) + variable + '.pickle'):
            with open(folder + check_underscore(state_prefix) + variable + '.pickle', 'rb') as f:
                loaded = pickle.load(f)
            return loaded
        else:
            return None
    else:
        if os.path.exists(folder + variable + '.pickle'):
            with open(folder + variable + '.pickle', 'rb') as f:
                loaded = pickle.load(f)
            return loaded
        else:
            return None

def load_fits_hdul(file_location):
    with fits.open(file_location) as hdul:
        fits_file = hdul
    return fits_file

def load_fits_as_dataframe(file_location):
    with fits.open(file_location) as hdul:
        df = Table(hdul[1].data).to_pandas()
    return df

def load_apj_table(file_location):
    data = Table.read(file_location, format="ascii.cds")
    return data.to_pandas()

def load_fits_as_table(file_location):
    with fits.open(file_location) as hdul:
        table = Table(hdul[1].data)
    return table

def get_matched(fn, thres=0.95):
    df = load_fits_as_dataframe(fn)
    matched = df.loc[
        (df['p_any'] > thres) & (df['p_i'] > thres)
    ]
    return matched

def check_strictly_positive(v):
    if v <= 0:
        return 10**-4
    else:
        return v

def load_dataset_by_name(to_be_loaded=['lotss_dr2', 'reines', 'psrcat']):
    #assert t in ['lotss_dr1', 'lotss_dr2', 'lotass', 'sanidas', 'reines', 'clu_dwarfs', 'chime', 'psrcat']

    dfs = {}
    for dataset in to_be_loaded:
        if dataset == 'lotss_dr2':
            df_lotss_dr2 = load_fits_as_dataframe('data/LoTSS/LoTSS_DR2_v100.srl.fits')
            dfs[dataset] = df_lotss_dr2

        if dataset == 'clu_dwarfs':
            df_clu = pd.read_csv('data/CLU/CLU_20190708_marshalFormat.csv')
            df_clu_dwarfs = df_clu.loc[
            #     (df_clu['mstar'].apply(np.log10) > 7) &
                (df_clu['mstar'] <= 3e9)
            #     (df_clu['cluhamag'] < 14)
            ]
            dfs[dataset] = df_clu_dwarfs

        if dataset == 'chime':
            df_chime = pd.read_csv('data/chimefrbcat1.csv')
            dfs[dataset] = df_chime

    return dfs

def get_matched(obj1, obj2, mask1=None, mask2=None):
    if mask1 is not None and mask2 is not None:
        matched_obj1 = obj1.df.iloc[mask1]
        matched_obj2 = obj2.df.iloc[mask2]
    else:
        matched_obj1 = obj1.df.iloc[obj2.matches[obj1.name].mask_idx]
        matched_obj2 = obj2.df.iloc[obj2.matches[obj1.name].filtered_idx]

    return matched_obj1, matched_obj2
