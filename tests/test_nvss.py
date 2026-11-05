"""NVSS column labels."""
from matchmaker.dataset.nvss import Nvss


def test_nvss_residual_flux_is_not_labelled_peak():
    cols = Nvss(load_data=False).cols.measure
    assert not hasattr(cols, 'peak_flux')
    assert cols.residual_flux.label == 'resFlux'


# --- VLASS --------------------------------------------------------------------------------
