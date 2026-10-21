import numpy as np
from numpy.testing import assert_allclose

from topwater import HullMesh
from .reference_hull import (
    CG, DISPLACEMENT, FULL_VOLUME, INERTIA, MASS, RIM_HEIGHT, SEMI_AXES, WATER_DENSITY,
    geometry_metrics, input_data, mesh_arrays, opening_boundaries,
)


def test_canonical_hull_converges_to_exact_geometry_and_inertia():
    errors = []
    for refinement in (2, 3, 4):
        vertices, faces = mesh_arrays(refinement)
        HullMesh(vertices, faces)
        # Closed-form cross-section integrals of the truncated ellipsoid;
        # the submerged half remains unchanged by the above-water cut.
        metrics = geometry_metrics(vertices, faces)
        assert_allclose(vertices.max(axis=0), [3, 1, 0.25])
        assert_allclose(vertices.min(axis=0), [-3, -1, -0.5])
        assert_allclose(metrics["centroid"][:2], [0, 0], atol=1e-14)
        assert_allclose(metrics["buoyancy_centroid"][:2], [0, 0], atol=1e-14)
        errors.append([
            abs(metrics["volume"] / FULL_VOLUME - 1),
            abs(metrics["centroid"][2] / CG[2] - 1),
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


def test_rim_and_artificial_cap_are_explicit():
    vertices, faces = mesh_arrays()
    rim, cap = opening_boundaries(vertices, faces)
    assert len(rim) >= 4
    assert len(cap) > 0
    assert_allclose(vertices[rim, 2], RIM_HEIGHT)
    assert_allclose(np.sum((vertices[rim] / SEMI_AXES)**2, axis=1), 1)
    assert_allclose(vertices[faces[cap], 2], RIM_HEIGHT)
    # Every consecutive rim segment bounds exactly one cap and one hull face.
    edges = np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owners = np.tile(np.arange(len(faces)), 3)
    is_cap = np.zeros(len(faces), dtype=bool)
    is_cap[cap] = True
    for start, end in zip(rim, np.roll(rim, -1)):
        adjacent = owners[np.all(edges == sorted([start, end]), axis=1)]
        assert len(adjacent) == 2
        assert np.count_nonzero(is_cap[adjacent]) == 1


def test_loading_and_pose_use_the_truncated_hull_cg():
    data = input_data()
    assert_allclose(CG, [0, 0, -1 / 16], atol=1e-15)
    assert_allclose(FULL_VOLUME, 27 * np.pi / 16)
    assert_allclose(data["body"]["cg"], CG)
    # Pose stores world CG; using the same coordinates keeps mesh waterline z=0.
    assert_allclose(data["initial_state"]["position"], CG)
    assert_allclose(INERTIA / MASS, [63 / 256, 2491 / 1280, 17 / 8])
