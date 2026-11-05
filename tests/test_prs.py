"""Reference PRS data (model/prs.py) and its link to plot styling (sources.py)."""
import pytest

from matchmaker.dataset import sources


def test_prs_anchor_values():
    from matchmaker.model import prs
    from matchmaker.model.bpt import log_ratio
    a = prs.FRB20121102A
    assert a['z'] == pytest.approx(0.19273)
    assert a['log_l144'] == pytest.approx(22.28, abs=0.01)          # Resmi+2021 250 MHz, flat alpha
    L = a['lines']
    y, _ = log_ratio(L['oiii_5007']['flux'], L['oiii_5007']['err'], L['h_beta']['flux'], L['h_beta']['err'])
    assert y == pytest.approx(0.659, abs=0.002)   # 0.688 as observed (Tendulkar+2017); -0.029 dex Galactic extinction
    assert L['nii_6584']['upper_limit'] and not L['sii_6717']['upper_limit']
    assert prs.FRB20190520B['lines'] is None                         # not yet extracted: callers must skip


def test_sources_take_shared_values_from_prs():
    from matchmaker.model import prs
    for plot, data in ((sources.frb121102, prs.FRB20121102A), (sources.frb190520b, prs.FRB20190520B)):
        assert plot['z_host'] == data['z'] and plot['sfr'] == data['sfr'] and plot['alpha'] == data['alpha']
