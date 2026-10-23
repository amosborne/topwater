"""Actual GPU batches against direct signed-lag Toeplitz multiplication."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

cp=pytest.importorskip('cupy')
from tests.experiments.free_surface.compact_p1_weak import WeakPatch


def test_batch_axis_and_single_column_share_one_weak_fft_action():
    shape=(7,8);operator=WeakPatch(shape,.31,.23)
    iy,ix=np.indices(shape);flat_y=iy.ravel();flat_x=ix.ravel()
    direct=operator.kernel[flat_y[:,None]-flat_y[None,:]+shape[0]-1,
                           flat_x[:,None]-flat_x[None,:]+shape[1]-1]
    values=np.random.default_rng(17).standard_normal((*shape,4))
    expected=(direct@values.reshape(-1,4)).reshape(*shape,4)
    actual=operator.apply_host(values)
    assert_allclose(actual,expected,rtol=3e-12,atol=2e-14)
    assert_allclose(cp.asnumpy(operator.apply(cp.asarray(values[:,:,2]))),expected[:,:,2],rtol=3e-12,atol=2e-14)
