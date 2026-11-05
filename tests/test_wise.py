"""Vectorised WISE AGN box (Jarrett+2011)."""
import numpy as np
import pandas as pd

from matchmaker.dataset.clu import Clu
from matchmaker.model import wise


def test_wise_agn_mask_matches_scalar_definition():
    rng = np.random.default_rng(0)
    w1w2, w2w3 = rng.uniform(-0.5, 2.5, 500), rng.uniform(0, 6, 500)
    vec = wise.agn_mask(w1w2, w2w3)
    assert list(vec) == [bool(wise.is_AGN(a, b)) for a, b in zip(w1w2, w2w3)]
    assert not wise.agn_mask([np.nan], [3.0])[0]


def test_clu_non_agn_mask_without_matplotlib():
    c = Clu(load_data=False)
    c.df = pd.DataFrame({'w1mpro': [10.0, 10.0], 'w2mpro': [9.0, 9.9], 'w3mpro': [6.0, 9.0]})
    assert list(c.non_agn_mask()) == [False, True]    # (1.0, 3.0) inside the box; (0.1, 0.9) outside
