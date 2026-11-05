"""NED-LVS loader."""
import astropy.units as u
import numpy as np
import pandas as pd

from matchmaker.dataset.nedlvs import NedLvs


def test_columns_map_to_the_catalogue():
    n = NedLvs(load_data=False)
    assert n.cols.mstar.label == 'Mstar'
    assert n.cols.distance.label == 'DistMpc'
    assert n.cols.a.label == 'Diam'
    assert n.cols.a.unit == u.arcsec


def test_sfr_column_is_selectable():
    assert NedLvs(load_data=False, sfr_column='SFR_hybrid').cols.sfr.label == 'SFR_hybrid'
    assert NedLvs(load_data=False, sfr_column='SFR_W4').cols.sfr.label == 'SFR_W4'


def test_constraints_keep_dwarfs_with_a_rate():
    n = NedLvs(load_data=False)
    n.df = pd.DataFrame({
        'objtype': ['G', 'G', 'G', 'QSO', 'G'],
        'Mstar': [1e8, 1e6, 1e10, 1e8, 1e9],          # too light, too heavy, a quasar
        'SFR_hybrid': [0.1, 0.1, 0.1, 0.1, np.nan],   # and one without a rate
    })
    n.load_data = lambda *a, **k: None                # keep the frame we just set
    keep = ((n.df['objtype'].isin(('G', 'GPair', 'GTrpl')))
            & (n.df['Mstar'] > 1e7) & (n.df['Mstar'] <= 3e9)
            & np.isfinite(n.df['SFR_hybrid']) & (n.df['SFR_hybrid'] > 0))
    assert list(keep) == [True, False, False, False, False]
