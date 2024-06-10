import os
import pickle
import bz2
import _pickle as cPickle
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

sort_dict_by_value = lambda x, reverse=False: dict(sorted(x.items(),
                                                          key=lambda item: item[1],
                                                          reverse=reverse))
share_items = lambda a, b: not set(a).isdisjoint(b)
intersection = lambda a, b: np.sort(list(set(a).intersection(b)))
difference = lambda a, b: np.sort(list(set(a).difference(b)))
union = lambda a, b: np.sort(list(set(a).union(b)))

# Functions
def check_underscore(string):
    """
    Checks if the given string ends with an underscore. If not, appends an underscore to the string.

    Args:
        string (str): The input string to check.

    Returns:
        str: The modified string with an underscore at the end, if necessary.
    """
    if string != '':
        if string[-1] != '_':
            string += '_'
    return string

def check_slash(string):
    """
    Adds a trailing slash to the input string if it doesn't already have one.

    Args:
        string (str): The input string to check.

    Returns:
        str: The input string with a trailing slash added, if necessary.
    """
    if string != '':
        if string[-1] != '/':
            string += '/'
    return string

def check_folder_exists_or_create(folder, return_folder=True):
    """
    Checks if a folder exists and creates it if it doesn't.
    
    Args:
        folder (str): The path of the folder to check/create.
        return_folder (bool, optional): Whether to return the folder path after creating it. 
                                        Defaults to True.
    
    Returns:
        str: The folder path if `return_folder` is True, otherwise None.
    """
    if not os.path.exists(folder):
        os.makedirs(folder)
    if return_folder:
        return folder

def check_file_exists(file):
    """
    Check if a file exists.

    Args:
        file (str): The path to the file.

    Returns:
        bool: True if the file exists, False otherwise.
    """
    return os.path.isfile(file)

def remove_file_from_folder(file_location):
    """
    Removes a file from the specified folder.

    Args:
        file_location (str): The location of the file to be removed.

    Raises:
        FileNotFoundError: If the file does not exist at the specified location.
        PermissionError: If the user does not have permission to remove the file.

    """
    os.remove(file_location)

def remove_file_if_exists(filename, return_filename=False):
    """
    Removes the specified file if it exists.

    Args:
        filename (str): The path of the file to be removed.
        return_filename (bool, optional): If True, the filename will be returned. 
            Defaults to False.

    Returns:
        str: The filename if `return_filename` is True, otherwise None.
    """
    if os.path.exists(filename):
        os.remove(filename)
    return filename if return_filename else None

def save(variable:str, data, protocol=pickle.HIGHEST_PROTOCOL, state_prefix='', folder='states/'):
    """
    Save the given data to a pickle file.

    Args:
        variable (str): The name of the variable being saved.
        data: The data to be saved.
        protocol (int, optional): The pickle protocol to use. Defaults to pickle.HIGHEST_PROTOCOL.
        state_prefix (str, optional): The prefix to be added to the filename. Defaults to ''.
        folder (str, optional): The folder where the pickle file will be saved. Defaults to 'states/'.
    """
    check_folder_exists_or_create(folder)

    if state_prefix != '':
        with open(check_slash(folder) + check_underscore(state_prefix) + variable + '.pickle', 'wb') as f:
            pickle.dump(data, f, protocol)
    else:
        with open(check_slash(folder) + variable + '.pickle', 'wb') as f:
            pickle.dump(data, f, protocol)

def load(variable, state_prefix='', folder='states/'):
    """
    Load a variable from a pickle file.

    Args:
        variable (str): The name of the variable to load.
        state_prefix (str, optional): The prefix to use for the state file. Defaults to ''.
        folder (str, optional): The folder where the state files are stored. Defaults to 'states/'.

    Returns:
        The loaded variable if the file exists, otherwise None.
    """
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

def save_compressed(variable:str, data, state_prefix='', folder='states/'):
    """
    Save the given data to a compressed file.

    Args:
        variable (str): The name of the variable being saved.
        data: The data to be saved.
        state_prefix (str, optional): The prefix to be added to the filename. Defaults to ''.
        folder (str, optional): The folder where the file will be saved. Defaults to 'states/'.
    """
    check_folder_exists_or_create(folder)

    if state_prefix != '':
        with bz2.BZ2File(check_slash(folder) + check_underscore(state_prefix) + variable + '.pbz2', 'w') as f:
            cPickle.dump(data, f)
    else:
        with bz2.BZ2File(check_slash(folder) + variable + '.pickle', 'wb') as f:
            cPickle.dump(data, f)

def load_compressed(variable, state_prefix='', folder='states/'):
    """
    Load a compressed variable from a file.

    Args:
        variable (str): The name of the variable to load.
        state_prefix (str, optional): The prefix to use for the state file. Defaults to ''.
        folder (str, optional): The folder where the state files are stored. Defaults to 'states/'.

    Returns:
        The loaded variable if it exists, otherwise None.
    """
    if state_prefix != '':
        if os.path.exists(folder + check_underscore(state_prefix) + variable + '.pbz2'):
            with bz2.BZ2File(folder + check_underscore(state_prefix) + variable + '.pbz2', 'rb') as f:
                loaded = cPickle.load(f)
            return loaded
        else:
            return None
    else:
        if os.path.exists(folder + variable + '.pbz2'):
            with bz2.BZ2File(folder + variable + '.pbz2', 'rb') as f:
                loaded = cPickle.load(f)
            return loaded
        else:
            return None

def load_fits_hdul(file_location):
    """
    Load a FITS file from the given file location.

    Parameters:
    file_location (str): The path to the FITS file.

    Returns:
    HDUList: The FITS file as an HDUList object.

    """
    with fits.open(file_location) as hdul:
        fits_file = hdul
    return fits_file

def load_fits_as_dataframe(file_location):
    """
    Load a FITS file and convert it into a pandas DataFrame.

    Parameters:
    - file_location (str): The path to the FITS file.

    Returns:
    - df (pandas.DataFrame): The FITS data as a DataFrame.
    """
    with fits.open(file_location) as hdul:
        df = Table(hdul[1].data).to_pandas()
    return df

def load_apj_table(file_location):
    """
    Load an APJ table from the given file location.

    Parameters:
    - file_location (str): The path to the file containing the APJ table.

    Returns:
    - pandas.DataFrame: The loaded APJ table as a pandas DataFrame.
    """
    data = Table.read(file_location, format="ascii.cds")
    return data.to_pandas()

def load_fits_as_table(file_location):
    """
    Load a FITS file as a table.

    Parameters:
        file_location (str): The path to the FITS file.

    Returns:
        astropy.table.Table: The FITS file data as a table.
    """
    with fits.open(file_location) as hdul:
        table = Table(hdul[1].data)
    return table

def get_matched(fn, thres=0.95):
    """
    Get matched rows from a DataFrame based on given threshold values.

    Parameters:
    fn (str): The file name or path of the FITS file.
    thres (float, optional): The threshold value for matching. Defaults to 0.95.

    Returns:
    pandas.DataFrame: The matched rows from the DataFrame.
    """
    df = load_fits_as_dataframe(fn)
    matched = df.loc[
        (df['p_any'] > thres) & (df['p_i'] > thres)
    ]
    return matched

def check_strictly_positive(v):
    """
    Checks if a value is strictly positive.

    Parameters:
    v (float): The value to be checked.

    Returns:
    float: The original value if it is strictly positive, otherwise returns 10**-4.
    """
    if v <= 0:
        return 10**-4
    else:
        return v

def get_matched(obj1, obj2, mask1=None, mask2=None):
    """
    Get matched objects from obj1 and obj2 based on provided masks.

    Args:
        obj1: The first object.
        obj2: The second object.
        mask1: Optional mask for obj1.
        mask2: Optional mask for obj2.

    Returns:
        Tuple: A tuple containing the matched objects from obj1 and obj2.
    """
    if mask1 is not None and mask2 is not None:
        matched_obj1 = obj1.df.iloc[mask1]
        matched_obj2 = obj2.df.iloc[mask2]
    else:
        matched_obj1 = obj1.df.iloc[obj2.matches[obj1.name].mask_idx]
        matched_obj2 = obj2.df.iloc[obj2.matches[obj1.name].filtered_idx]

    return matched_obj1, matched_obj2

def eval_flux_given_spectral_index(freq1, freq2, flux1, alpha):
    """
    Calculate the flux at a given frequency, given the spectral index.

    Parameters:
    freq1 (float): The reference frequency.
    freq2 (float): The target frequency.
    flux1 (float): The flux at the reference frequency.
    alpha (float): The spectral index.

    Returns:
    float: The calculated flux at the target frequency.
    """
    return flux1 * (freq2/freq1)**alpha
