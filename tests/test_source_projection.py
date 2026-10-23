"""Analytic low/high mode separation, including vector-valued source columns."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')
from tests.experiments.free_surface.boat_source_projection import project_trace
from tests.experiments.free_surface.cuda_graph import cp


def test_projection_preserves_modes_and_rejects_high_mode_coarse_alias():
    size, samples = 17, 33
    x, y = np.meshgrid(2*np.pi*np.arange(samples)/samples,
                       2*np.pi*np.arange(samples)/samples)
    values = np.column_stack(((3+.7*np.cos(3*x-2*y)+.4*np.sin(10*x+y)).ravel(),
                              (2*np.sin(x+4*y)+.3*np.cos(9*y)).ravel()))
    xx, yy = np.meshgrid(2*np.pi*np.arange(size)/size, 2*np.pi*np.arange(size)/size)
    expected = np.column_stack(((3+.7*np.cos(3*xx-2*yy)).ravel(),
                                (2*np.sin(xx+4*yy)).ravel()))
    assert_allclose(cp.asnumpy(project_trace(values, samples, size)), expected, atol=2e-14, rtol=0)
    assert_allclose(cp.asnumpy(project_trace(values[:, 0], samples, size)), expected[:, 0], atol=2e-14, rtol=0)
