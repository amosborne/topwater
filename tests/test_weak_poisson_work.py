"""Small independent-work controls; no production solver or FFT kernel."""
import numpy as np
from numpy.testing import assert_allclose
import pytest
from tests.experiments.free_surface.compact_p1_trace import CompactP1Trace
from tests.experiments.free_surface.weak_poisson_work import hat_pair, h1_seminorm_squared, poisson_work


def test_hat_seminorm_and_zero_polarization():
    hx, hy = .31, .23
    plus, zero = hat_pair(hx, hy, [0, 0])
    unit = CompactP1Trace(plus.x, plus.y, plus.values/2)
    expected = 2*(hx/hy+hy/hx)
    assert_allclose(h1_seminorm_squared(unit), expected, rtol=2e-15)
    assert_allclose(h1_seminorm_squared(plus), 4*expected, rtol=2e-15)
    result = poisson_work(zero, .02, absolute_tolerance=1e-5, max_evaluations=10)
    assert result['work'] == 0 and result['target_evaluations'] == 0
    assert result['inserted_continuum_interval'] == [0., 0.]


def test_positive_depth_work_scaling_and_monotonicity():
    plus, _ = hat_pair(.31, .23, [0, 0])
    unit = CompactP1Trace(plus.x, plus.y, plus.values/2)
    deep = poisson_work(unit, .08, absolute_tolerance=2e-5, max_evaluations=20000)
    shallow = poisson_work(unit, .04, absolute_tolerance=2e-5, max_evaluations=20000)
    scaled = poisson_work(plus, .08, absolute_tolerance=8e-5, max_evaluations=20000)
    assert 0 < deep['work'] < shallow['work']
    assert shallow['work'] <= deep['inserted_continuum_interval'][1]
    assert_allclose(scaled['work'], 4*deep['work'], atol=2e-8, rtol=2e-8)
    assert not deep['certified_total_error']
    with pytest.raises(RuntimeError, match='budget'):
        poisson_work(unit, .04, absolute_tolerance=1e-5, max_evaluations=2)
