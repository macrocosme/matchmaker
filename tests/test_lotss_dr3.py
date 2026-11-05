"""Lotss(release='dr3') on a small catalogue with the DR3 column layout."""
import numpy as np
import pytest
from astropy.table import Table

pytest.importorskip('astropy_healpix')

from matchmaker.dataset.lotss import Lotss  # noqa: E402


@pytest.fixture
def dr3_file(tmp_path):
    rng = np.random.default_rng(1)
    n = 50_000
    ra = rng.uniform(150, 170, n)
    dec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(0)), np.sin(np.radians(30)), n)))
    # DR3 stores the 6" beam with float noise; reproduce it
    res = np.where(dec < 10, 9.0, rng.choice([6.000000000000001, 5.999999999999997], n, p=[.99, .01]))
    rms = np.full(n, 0.08)
    tot = rms * rng.uniform(3, 100, n)
    zeros = np.zeros(n)
    t = Table({'Source_Name': np.array(['ILTJx'] * n), 'RA': ra, 'E_RA': zeros + .3, 'DEC': dec, 'E_DEC': zeros + .3,
               'Total_flux': tot, 'E_Total_flux': .1 * tot, 'Peak_flux': .9 * tot, 'E_Peak_flux': rms,
               'Maj': zeros + 7, 'E_Maj': zeros, 'Min': zeros + 6, 'E_Min': zeros, 'PA': zeros, 'E_PA': zeros,
               'DC_Maj': zeros, 'E_DC_Maj': zeros, 'DC_Min': zeros, 'E_DC_Min': zeros, 'DC_PA': zeros,
               'E_DC_PA': zeros, 'Isl_Total_flux': tot, 'E_Isl_Total_flux': .1 * tot, 'Isl_rms': rms,
               'S_Code': rng.choice(['S', 'M'], n, p=[.8, .2]), 'HEALPIX': np.zeros(n, dtype=np.int64),
               'Resolution': res})
    p = tmp_path / 'dr3.fits'
    t.write(p)
    return str(p)


def _lotss(path, **kw):
    lot = Lotss(load_data=False, release='dr3')
    lot.file_location = path
    lot.load_data(**kw)
    return lot


def test_dr3_columns():
    lot = Lotss(load_data=False, release='dr3')
    assert lot.cols.radio.resolution.label == 'Resolution'
    assert not hasattr(lot.cols.radio, 'mosaic_id')
    assert lot.area is None and lot.file_location.endswith('LoTSS_DR3_v1.0.srl.fits')


def test_dr2_is_unchanged_default():
    lot = Lotss(load_data=False)
    assert lot.release == 'dr2' and lot.cols.radio.mosaic_id.label == 'Mosaic_ID'
    assert np.isclose(lot.area.value, 5635)


def test_dr3_footprint_built_before_cuts(dr3_file):
    lot = _lotss(dr3_file, constrain=True, smin=5.0, single_only=True)
    assert (lot.df.Total_flux >= 5.0).all() and (lot.df.S_Code == 'S').all()
    # footprint reflects every source, not only the bright single ones
    true = np.radians(20) * np.sin(np.radians(30)) * (180 / np.pi) ** 2
    assert true < lot.area.value < 1.2 * true
    assert lot.footprint.contains([160], [15])[0] and not lot.footprint.contains([160], [50])[0]
    assert lot.boxes is None


def test_dr3_resolution_per_source(dr3_file):
    lot = _lotss(dr3_file)
    r = lot.resolution()
    assert set(np.unique(r)) == {6.0, 9.0}          # float noise on the 6" beam is rounded away
    assert np.all(r[lot.df.DEC.values < 10] == 9.0)
    assert lot.footprint.area_where(6.0).value > 0 and lot.footprint.area_where(9.0).value > 0


def test_dr3_compact_only_refuses_dr2_envelope(dr3_file):
    with pytest.raises(NotImplementedError):
        _lotss(dr3_file, constrain=True, compact_only=True)


def test_constrain_without_cuts_no_longer_crashes(dr3_file):
    lot = _lotss(dr3_file, constrain=True)
    assert len(lot.df) == 50_000


def test_min_snr_cut(dr3_file):
    lot = _lotss(dr3_file, constrain=True, min_snr=9.0)
    assert (lot.df.Total_flux / lot.df.Isl_rms >= 9.0).all()
    assert 0 < len(lot.df) < 50_000
