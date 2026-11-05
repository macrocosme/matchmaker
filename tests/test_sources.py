"""Reference-object values."""
import pytest

from matchmaker.dataset import sources


def test_frb121102_host_redshift():
    assert sources.frb121102['z_host'] == pytest.approx(0.19273)


# --- BPT ----------------------------------------------------------------------------------
