"""BPT line ratios and demarcations."""
import math

import numpy as np

from matchmaker.model import bpt


def test_log_ratio_error_propagation():
    r, e = bpt.log_ratio(120.0, 4.0, 35.0, 3.0)
    assert np.isclose(r, np.log10(120 / 35))
    assert np.isclose(e, math.hypot(4 / 120, 3 / 35) / math.log(10))
    assert e > 0


def test_bpt_classification():
    cls = bpt.classify_nii([-1.0, -0.2, 0.3, np.nan], [0.0, 0.0, 1.0, 0.0])   # Kewley at x=-0.2 is 0.28
    assert list(cls) == ['star-forming', 'composite', 'AGN', '']


# --- spectral index -----------------------------------------------------------------------


def test_seyfert_liner_line_starts_on_the_starburst_line():
    """The Seyfert/LINER division only applies above Kewley+2001; it must not cross below it."""
    import numpy as np

    from matchmaker.model.bpt import (
        SII_SEYFERT_LINER_START,
        kewley01_sii,
        kewley06_sii_seyfert_liner,
        kewley06_sii_valid,
    )

    x0 = SII_SEYFERT_LINER_START
    assert np.isclose(kewley06_sii_seyfert_liner(x0), kewley01_sii(x0), atol=1e-9)

    x = np.linspace(x0 + 1e-6, 0.31, 50)          # inside the valid range, the line is above the curve
    assert np.all(kewley06_sii_seyfert_liner(x) >= kewley01_sii(x) - 1e-9)
    assert np.all(kewley06_sii_valid(x))

    below = np.linspace(-1.5, x0 - 0.01, 50)
    assert not np.any(kewley06_sii_valid(below))
