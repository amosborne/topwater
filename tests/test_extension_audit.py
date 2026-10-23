"""Known biharmonic fields through the actual reference waterplane hole."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial import ConvexHull

from topwater.hydrostatics import hydrostatics
from tests.reference_boat import simulation_inputs
from tests.experiments.free_surface.boat_extension_audit import biharmonic_fill


def test_extension_reproduces_nonharmonic_biharmonic_polynomials():
    inputs = simulation_inputs()
    wet = hydrostatics(inputs.body, inputs.initial_state, density=1025.)
    hull = ConvexHull(wet.waterline_segments[:, :, :2].reshape(-1, 2))
    size, length = 33, 16.
    origin = inputs.initial_state.position[:2]-length/2
    x, y = np.meshgrid(origin[0]+np.arange(size)*length/size,
                        origin[1]+np.arange(size)*length/size)
    xy = np.column_stack((x.ravel(), y.ravel()))
    hidden = np.flatnonzero(np.all(xy@hull.equations[:, :2].T+hull.equations[:, 2] <= 0, axis=1))
    boundary, weights, diagnostics = biharmonic_fill(size, hidden)
    fields = np.stack((1+.2*x-.1*y, x*x+y*y, x*x*y), axis=-1).reshape(-1, 3)
    assert_allclose(weights@fields[boundary], fields[hidden], atol=3e-12, rtol=3e-12)
    assert diagnostics['rank'] == len(hidden)
    assert diagnostics['stationarity_relative_residual'] < 1e-12
