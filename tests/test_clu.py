"""Clu.semi_minor with array masks."""
import numpy as np
import pandas as pd

from matchmaker.dataset.clu import Clu


def test_clu_semi_minor_accepts_array_mask():
    c = Clu(load_data=False)
    c.df = pd.DataFrame({'a': [1.0, 2.0, 3.0], 'b2a': [0.5, np.nan, 0.2]})
    assert np.allclose(c.semi_minor(mask=np.array([0, 2])), [0.5, 0.6])
