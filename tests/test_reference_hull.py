import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.testing import assert_allclose

from topwater import HullMesh
from .reference_hull import (
    AREAL_DENSITY, CG, DISPLACEMENT, FULL_VOLUME, LOADING, MASS, PAYLOAD_MASS, PAYLOAD_POSITION,
    RIM_HEIGHT, SEMI_AXES, WATER_DENSITY, WATERLINE_HEIGHT,
    geometry_metrics, input_data, mesh_arrays, mesh_loading, opening_boundaries,
    reference_loading,
)


def test_canonical_hull_converges_to_exact_geometry_and_inertia():
    errors = []
    for refinement in (2, 3, 4):
        vertices, faces = mesh_arrays(refinement)
        HullMesh(vertices, faces)
        # Geometry is checked against exact half-ellipsoid integrals. Shell
        # loading is separately checked against parametric surface quadrature.
        metrics = geometry_metrics(vertices, faces)
        assert_allclose(vertices.max(axis=0), [1.75, 0.45, 0])
        assert_allclose(vertices.min(axis=0), [-1.75, -0.45, -0.35])
        assert_allclose(metrics["centroid"][:2], [0, 0], atol=1e-14)
        loading = mesh_loading(vertices, faces)
        assert_allclose(loading["cg"][:2], [0, 0], atol=1e-14)
        off_diagonal = loading["inertia_cg"] - np.diag(np.diag(loading["inertia_cg"]))
        assert_allclose(off_diagonal, 0, atol=1e-12)
        errors.append([
            abs(metrics["volume"] / FULL_VOLUME - 1),
            abs(metrics["centroid"][2] / (-3 * SEMI_AXES[2] / 8) - 1),
            abs(loading["shell_mass"] / LOADING["shell_mass"] - 1),
            abs(loading["shell_cg"][2] / LOADING["shell_cg"][2] - 1),
            abs(loading["cg"][2] - CG[2]) / SEMI_AXES[2],
            *np.abs(np.diag(loading["inertia_cg"]) / np.diag(LOADING["inertia_cg"]) - 1),
        ])
    errors = np.array(errors)
    # Piecewise-planar approximation: second-order convergence as edge lengths
    # halve, before any fluid solver or empirical damping is introduced.
    assert np.all(errors[1:] < errors[:-1])
    assert np.all(errors[:-1] / errors[1:] > 3)
    assert np.all(errors[-1] < 0.02)


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


def test_independent_quadrature_reference_is_resolved():
    refined = reference_loading(order=128)
    for name in ("shell_mass", "shell_cg", "mass", "cg", "inertia_cg"):
        assert_allclose(refined[name], LOADING[name], rtol=1e-10, atol=1e-12)


def test_inertia_by_direct_integration_about_combined_cg():
    # Different parameterization (mu=-z/c and full azimuth), no symmetry
    # cancellation assumptions, and no shared parallel-axis helper.
    nodes, weights = leggauss(128)
    mu, phi = np.meshgrid((nodes+1)/2, (nodes+1)*np.pi, indexing="ij")
    a, b, c = SEMI_AXES
    radius = np.sqrt(1-mu*mu)
    points = np.stack([a*radius*np.cos(phi), b*radius*np.sin(phi), -c*mu], axis=-1)
    jacobian = np.sqrt((b*c*radius*np.cos(phi))**2
                       + (a*c*radius*np.sin(phi))**2 + (a*b*mu)**2)
    dm = AREAL_DENSITY * jacobian * np.outer(weights/2, weights*np.pi)
    centred = points - CG
    squared_radius = np.sum(centred**2, axis=-1)
    tensor = (squared_radius[:, :, None, None] * np.eye(3)
              - centred[:, :, :, None] * centred[:, :, None, :])
    inertia = np.sum(dm[:, :, None, None] * tensor, axis=(0, 1))
    offset = PAYLOAD_POSITION - CG
    inertia += PAYLOAD_MASS * (np.dot(offset, offset) * np.eye(3) - np.outer(offset, offset))
    assert_allclose(inertia, LOADING["inertia_cg"], rtol=1e-9, atol=1e-10)


def test_loading_and_pose_use_the_combined_shell_and_payload_cg():
    data = input_data()
    assert_allclose(MASS, LOADING["shell_mass"] + PAYLOAD_MASS)
    assert_allclose(MASS * CG, LOADING["shell_mass"] * LOADING["shell_cg"]
                    + PAYLOAD_MASS * PAYLOAD_POSITION)
    assert_allclose(data["body"]["cg"], CG)
    # Initial position stores world CG. The body-coordinate reference waterline
    # must map to world z=0, not the mesh origin or the opening rim.
    world_waterline = (np.array([0, 0, WATERLINE_HEIGHT]) - CG
                       + data["initial_state"]["position"])
    assert_allclose(world_waterline, [0, 0, 0], atol=1e-15)
    assert -SEMI_AXES[2] < WATERLINE_HEIGHT < RIM_HEIGHT


def test_reference_loading_has_positive_initial_roll_stability():
    # Independent closed-form waterplane moment and submerged first moment.
    # This checks the chosen reference loading, not a production buoyancy solver.
    a, b, c = SEMI_AXES
    h = WATERLINE_HEIGHT
    waterplane_moment = np.pi * a * b**3 * (1 - (h/c)**2)**2 / 4
    buoyancy_z = np.pi * a * b * (h*h/2 - h**4/(4*c*c) - c*c/4) / DISPLACEMENT
    metacentric_height = buoyancy_z + waterplane_moment / DISPLACEMENT - CG[2]
    assert metacentric_height > 0
    assert_allclose(DISPLACEMENT * WATER_DENSITY, MASS)
