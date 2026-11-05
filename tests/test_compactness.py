"""Compactness envelope (matchmaker.model.compactness)."""
import numpy as np
import pytest

from matchmaker.model import compactness as cp


def test_envelope_matches_the_published_dr2_formula():
    snr = np.array([5.0, 50.0, 500.0])
    expected = 0.42 + 1.08 / (1 + (snr / 96.57) ** 2.49)
    assert np.allclose(cp.envelope(snr, *cp.ENVELOPE_DR2), expected)
    assert list(cp.is_compact([0.1, 3.0], [50.0, 50.0])) == [True, False]


def _catalogue(n=400_000, centre=0.25, width=2.0, resolved_fraction=0.3, seed=0):
    """Unresolved sources scattering about `centre` with sigma = width/snr, plus a resolved tail."""
    rng = np.random.default_rng(seed)
    snr = 10 ** rng.uniform(np.log10(5), np.log10(2000), n)
    r = rng.normal(centre, width / snr)
    resolved = rng.random(n) < resolved_fraction
    r[resolved] += rng.exponential(0.8, resolved.sum())      # only ever positive
    return r, snr, resolved


@pytest.mark.parametrize('method', ['tail', 'mode'])
def test_measured_envelope_recovers_the_unresolved_tail(method):
    r, snr, _ = _catalogue()
    centre, r999, counts, modes = cp.measure_envelope(r, snr, min_per_bin=200, method=method)
    truth = 0.25 + 3.09 * 2.0 / centre            # 0.1st percentile of a Gaussian is -3.09 sigma
    assert len(centre) >= 10
    assert np.allclose(modes, 0.25, atol=0.08)     # the mode, not the median, finds the centre
    assert np.abs(r999 - truth).mean() < 0.05 and np.abs(r999 - truth).max() < 0.2


def test_fit_recovers_the_envelope_and_flags_resolved_sources():
    r, snr, resolved = _catalogue()
    params, diag = cp.fit_envelope(r, snr, min_per_bin=200)
    assert params[0] == pytest.approx(0.25, abs=0.05)         # asymptote at high S/N
    snr_probe = np.array([10.0, 100.0, 1000.0])
    assert np.allclose(cp.envelope(snr_probe, *params), 0.25 + 3.09 * 2.0 / snr_probe, atol=0.15)
    compact = cp.is_compact(r, snr, params)
    # a smooth fit to per-bin 99.9% values keeps ~99% overall, not 99.9%
    assert compact[~resolved].mean() > 0.98
    assert compact[resolved].mean() < 0.3                     # and rejects most resolved ones


def test_tail_estimator_survives_a_contaminated_mode():
    """A second population of marginally-resolved sources captures the mode; the low tail stays clean.

    This is what DR3's 9" imaging shows above S/N ~ 50, where the mode jumps by ~0.1 dex.
    """
    rng = np.random.default_rng(4)
    n = 400_000
    snr = 10 ** rng.uniform(np.log10(5), np.log10(2000), n)
    r = rng.normal(0.25, 2.0 / snr)
    resolved = rng.random(n) < 0.6                      # majority, but the 10th percentile is still unresolved
    r[resolved] = rng.normal(0.60, 0.05, resolved.sum())    # narrow bump that dominates the histogram
    centre, e_tail, _, modes = cp.measure_envelope(r, snr, min_per_bin=200, method='tail')
    _, e_mode, _, _ = cp.measure_envelope(r, snr, min_per_bin=200, method='mode')
    truth = 0.25 + 3.09 * 2.0 / centre
    assert modes[:5].mean() > 0.5                        # at low S/N the mode is captured by the bump
    assert np.abs(e_tail - truth).mean() < np.abs(e_mode - truth).mean() / 2


def test_median_would_be_wrong():
    """The median centre is pulled up by resolved sources; the mode is not."""
    r, snr, _ = _catalogue(resolved_fraction=0.8)
    b = (snr > 100) & (snr < 300)
    assert np.median(r[b]) > 0.40 and cp._mode(r[b]) == pytest.approx(0.25, abs=0.05)


def test_fit_needs_enough_bins():
    with pytest.raises(ValueError):
        cp.fit_envelope(np.zeros(20), np.full(20, 10.0))
