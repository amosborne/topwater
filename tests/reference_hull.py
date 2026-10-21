"""One canonical numerical hull, shared by tests, timings and figures.

The reference is a triaxial ellipsoid, not a measured boat. Refinement projects
subdivided triangles onto the same exact surface; the equator stays at z=0.
This is test infrastructure, not a production mesh generator or hydrostatics API.
"""

import numpy as np


SEMI_AXES = np.array([3.0, 1.0, 0.5])  # metres; bow/port/up
WATER_DENSITY = 1000.0  # kg/m^3 for this reference case only
FULL_VOLUME = 4 * np.pi * np.prod(SEMI_AXES) / 3
DISPLACEMENT = FULL_VOLUME / 2
MASS = WATER_DENSITY * DISPLACEMENT
INERTIA = MASS / 5 * (np.sum(SEMI_AXES**2) - SEMI_AXES**2)


def mesh_arrays(refinement=2):
    """Return contiguous vertices/faces; there are 8*4**refinement panels."""
    if isinstance(refinement, bool) or not isinstance(refinement, int) or refinement < 0:
        raise ValueError("refinement must be a nonnegative integer")
    vertices = np.array([
        [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
    ], dtype=float)
    faces = np.array([
        [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
        [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5],
    ], dtype=np.int64)
    for _ in range(refinement):
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        unique, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        midpoints = vertices[unique].mean(axis=1)
        midpoints /= np.linalg.norm(midpoints, axis=1)[:, None]
        ab, bc, ca = (inverse + len(vertices)).reshape(3, -1)
        a, b, c = faces.T
        faces = np.vstack([
            np.column_stack([a, ab, ca]), np.column_stack([ab, b, bc]),
            np.column_stack([ca, bc, c]), np.column_stack([ab, bc, ca]),
        ])
        vertices = np.vstack([vertices, midpoints])
    return np.ascontiguousarray(vertices * SEMI_AXES), np.ascontiguousarray(faces)


def input_data(refinement=2):
    """Analytical calm-water reference pose plus explicit synthetic wave inputs.

    CG is the ellipsoid centre; a uniform body density of 500 kg/m^3 supplies
    the declared mass and inertia. Mass balances the exact half-ellipsoid's
    displacement; the finite mesh slightly under-displaces. No equilibrium
    has been solved. Initial velocities are zero. No random phases are generated
    here; the seed is illustrative provenance metadata.
    """
    vertices, faces = mesh_arrays(refinement)
    return {
        "format_version": 1,
        "body": {
            "mesh": {"vertices": vertices.tolist(), "faces": faces.tolist()},
            "mass": float(MASS), "cg": [0, 0, 0],
            "inertia_cg": np.diag(INERTIA).tolist(),
        },
        "waves": {
            "amplitudes": [0.1, 0.05], "angular_frequencies": [1.2, 2.4],
            "directions": [0, np.pi / 2], "phases": [0.3, -0.7], "seed": 42,
            "provenance": {"source": "synthetic reference case",
                           "phase_generation": "explicit; seed is illustrative metadata"},
        },
        "initial_state": {
            "position": [0, 0, 0], "orientation": [1, 0, 0, 0],
            "linear_velocity": [0, 0, 0], "angular_velocity": [0, 0, 0],
        },
    }


def surface_integrals(vertices, faces):
    """Independent divergence-theorem reference integrals of planar triangles.

    Uses z*n_z for volume (the input validator uses tetrahedral determinants).
    First and diagonal second moments use exact triangle polynomial averages.
    The lower half may omit its z=0 cap, whose contributions here are zero.
    """
    triangles = vertices[faces]
    p, q, r = triangles.transpose(1, 0, 2)
    normal_area = np.cross(q - p, r - p) / 2
    volume = np.sum(normal_area[:, 2] * triangles[:, :, 2].mean(axis=1))
    average_squared = (p*p + q*q + r*r + p*q + p*r + q*r) / 6
    first_moments = np.sum(normal_area * average_squared, axis=0) / 2
    average_cubed = (
        p**3 + q**3 + r**3 + p*q*(p+q) + p*r*(p+r) + q*r*(q+r) + p*q*r
    ) / 10
    second_moments = np.sum(normal_area * average_cubed, axis=0) / 3
    return volume, first_moments / volume, second_moments


def geometry_metrics(vertices, faces):
    """Evaluate full/half geometry and waterplane for the canonical pose."""
    volume, centroid, second = surface_integrals(vertices, faces)
    lower = faces[np.all(vertices[faces, 2] <= 0, axis=1)]
    displacement, buoyancy_centroid, _ = surface_integrals(vertices, lower)
    equator = vertices[vertices[:, 2] == 0]
    angles = np.arctan2(equator[:, 1], equator[:, 0])
    x, y = equator[np.argsort(angles), :2].T
    waterplane_area = np.sum(x * np.roll(y, -1) - y * np.roll(x, -1)) / 2
    inertia = (MASS / FULL_VOLUME) * (np.sum(second) - second)
    return dict(volume=volume, centroid=centroid, displacement=displacement,
                buoyancy_centroid=buoyancy_centroid, waterplane_area=waterplane_area,
                inertia_diagonal=inertia)
