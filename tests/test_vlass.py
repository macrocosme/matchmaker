"""VLASS undefined (-99) positional errors."""
import numpy as np
import pandas as pd

from matchmaker.dataset.vlass import Vlass


def test_vlass_undefined_errors_keep_position(tmp_path):
    p = tmp_path / 'vlass.csv'
    pd.DataFrame({'RA': [10.0, 20.0], 'DEC': [5.0, -5.0], 'E_RA': [1e-4, -99], 'E_DEC': [1e-4, -99],
                  'Total_flux': [2.0, 3.0]}).to_csv(p, index=False)
    v = Vlass(load_data=False)
    v.file_location = str(p)
    v.load_data()
    assert list(v.df.RA) == [10.0, 20.0] and list(v.df.Total_flux) == [2.0, 3.0]
    assert np.isnan(v.df.E_RA[1]) and np.isnan(v.df.E_DEC[1])


# --- chance association -------------------------------------------------------------------
