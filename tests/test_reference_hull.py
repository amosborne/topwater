import numpy as np
from numpy.testing import assert_allclose

from topwater import HullMesh
from .reference_hull import (
    DISPLACEMENT, FULL_VOLUME, INERTIA, MASS, SEMI_AXES, WATER_DENSITY,
    geometry_metrics, mesh_arrays,
)


def test_canonical_hull_converges_to_exact_geometry_and_inertia():
    errors = []
    for refinement in (2, 3, 4):
        vertices, faces = mesh_arrays(refinement)
        HullMesh(vertices, faces)
        # Closed-form references: V=4*pi*a*b*c/3, A=pi*a*b,
        # z_CB=-3*c/8 and I_xx=m*(b^2+c^2)/5, cyclically.
        metrics = geometry_metrics(vertices, faces)
        assert_allclose(vertices.max(axis=0), [3, 1, 0.5])
        assert_allclose(vertices.min(axis=0), [-3, -1, -0.5])
        assert_allclose(metrics["centroid"], [0, 0, 0], atol=1e-14)
        assert_allclose(metrics["buoyancy_centroid"][:2], [0, 0], atol=1e-14)
        assert_allclose(metrics["displacement"], metrics["volume"] / 2, rtol=1e-14)
        errors.append([
            abs(metrics["volume"] / FULL_VOLUME - 1),
            abs(metrics["displacement"] / DISPLACEMENT - 1),
            abs(metrics["buoyancy_centroid"][2] / (-3 * SEMI_AXES[2] / 8) - 1),
            abs(metrics["waterplane_area"] / (np.pi * 3 * 1) - 1),
            *np.abs(metrics["inertia_diagonal"] / INERTIA - 1),
        ])
    errors = np.array(errors)
    # Piecewise-planar approximation: second-order convergence as edge lengths
    # halve, before any fluid solver or empirical damping is introduced.
    assert np.all(errors[1:] < errors[:-1])
    assert np.all(errors[:-1] / errors[1:] > 3)
    assert np.all(errors[-1] < 0.02)
    assert_allclose(MASS / WATER_DENSITY, np.pi)  # Exact reference displacement.
