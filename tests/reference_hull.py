"""One canonical numerical hull, shared by tests, timings and figures.

The reference is an ellipsoid truncated above z=0, not a measured boat.
The opening rim is at z=0.25 m. A marked artificial cap closes the geometry;
it is not an impermeable deck. Refinement preserves that rim and the equator.
This is test infrastructure, not a production mesh generator or hydrostatics API.
"""

import numpy as np


SEMI_AXES = np.array([3.0, 1.0, 0.5])  # metres; bow/port/up
RIM_HEIGHT = 0.25  # metres above the reference waterline, in body coordinates
WATER_DENSITY = 1000.0  # kg/m^3 for this reference case only
# Integrate elliptical cross sections A(z)=pi*a*b*(1-z^2/c^2) from -c to h.
# The second-moment integrals below follow the same independent exact domain.
a, b, c = SEMI_AXES
u = RIM_HEIGHT / c
f0 = u - u**3 / 3 + 2 / 3
f1 = u**2 / 2 - u**4 / 4 - 1 / 4
f2 = u - 2 * u**3 / 3 + u**5 / 5 + 8 / 15
fz2 = u**3 / 3 - u**5 / 5 + 2 / 15
FULL_VOLUME = np.pi * a * b * c * f0
CG = np.array([0.0, 0.0, c * f1 / f0])
DISPLACEMENT = 2 * np.pi * a * b * c / 3
MASS = WATER_DENSITY * DISPLACEMENT
BODY_DENSITY = MASS / FULL_VOLUME  # Synthetic uniform loading, not water inside.
mean_squared = np.array([a*a*f2/(4*f0), b*b*f2/(4*f0), c*c*fz2/f0]) - CG**2
INERTIA = MASS * (np.sum(mean_squared) - mean_squared)


def mesh_arrays(refinement=2):
    """Return contiguous vertices/faces; there are 16*4**refinement panels."""
    if isinstance(refinement, bool) or not isinstance(refinement, int) or refinement < 0:
        raise ValueError("refinement must be a nonnegative integer")
    radius = np.sqrt(1 - u*u)
    vertices = np.array([
        [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1],
        [radius, 0, u], [-radius, 0, u], [0, radius, u], [0, -radius, u], [0, 0, u],
    ])
    faces = np.array([
        [2, 0, 4], [1, 2, 4], [3, 1, 4], [0, 3, 4],
        [0, 2, 7], [0, 7, 5], [2, 1, 6], [2, 6, 7],
        [1, 3, 8], [1, 8, 6], [3, 0, 5], [3, 5, 8],
        [9, 5, 7], [9, 7, 6], [9, 6, 8], [9, 8, 5],
    ], dtype=np.int64)
    cap_faces = np.all(vertices[faces, 2] == u, axis=1)
    for _ in range(refinement):
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        unique, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        midpoints = vertices[unique].mean(axis=1)
        curved = np.zeros(len(unique), dtype=bool)
        np.logical_or.at(curved, inverse, np.tile(~cap_faces, 3))
        rim = curved & np.all(vertices[unique, 2] == u, axis=1)
        midpoints[curved] /= np.linalg.norm(midpoints[curved], axis=1)[:, None]
        midpoints[rim, :2] *= radius / np.linalg.norm(midpoints[rim, :2], axis=1)[:, None]
        midpoints[rim, 2] = u
        ab, bc, ca = (inverse + len(vertices)).reshape(3, -1)
        a, b, c = faces.T
        faces = np.vstack([
            np.column_stack([a, ab, ca]), np.column_stack([ab, b, bc]),
            np.column_stack([ca, bc, c]), np.column_stack([ab, bc, ca]),
        ])
        vertices = np.vstack([vertices, midpoints])
        cap_faces = np.tile(cap_faces, 4)
    return np.ascontiguousarray(vertices * SEMI_AXES), np.ascontiguousarray(faces)


def opening_boundaries(vertices, faces):
    """Identify the ordered open rim and artificial closure faces explicitly.

    Future moving-geometry checks must terminate at rim immersion, before an
    artificial cap can produce pressure loads. No flooding detector exists yet.
    """
    on_surface = np.isclose(np.sum((vertices / SEMI_AXES)**2, axis=1), 1, rtol=0, atol=1e-12)
    rim = np.flatnonzero((vertices[:, 2] == RIM_HEIGHT) & on_surface)
    angles = np.arctan2(vertices[rim, 1], vertices[rim, 0])
    cap = np.flatnonzero(np.all(vertices[faces, 2] == RIM_HEIGHT, axis=1))
    return rim[np.argsort(angles)], cap


def input_data(refinement=2):
    """Analytical calm-water reference pose plus explicit synthetic wave inputs.

    A uniform synthetic loading of the truncated envelope supplies the declared
    mass, CG and inertia. It does not represent flood water or a real canoe's
    mass distribution. Mass balances the exact lower half-ellipsoid's
    displacement; the finite mesh slightly under-displaces. No equilibrium
    has been solved. Initial velocities are zero. No random phases are generated
    here; the seed is illustrative provenance metadata.
    """
    vertices, faces = mesh_arrays(refinement)
    return {
        "format_version": 1,
        "body": {
            "mesh": {"vertices": vertices.tolist(), "faces": faces.tolist()},
            "mass": float(MASS), "cg": CG.tolist(),
            "inertia_cg": np.diag(INERTIA).tolist(),
        },
        "waves": {
            "amplitudes": [0.1, 0.05], "angular_frequencies": [1.2, 2.4],
            "directions": [0, np.pi / 2], "phases": [0.3, -0.7], "seed": 42,
            "provenance": {"source": "synthetic reference case",
                           "phase_generation": "explicit; seed is illustrative metadata"},
        },
        "initial_state": {
            "position": CG.tolist(), "orientation": [1, 0, 0, 0],
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
    about_cg = second - 2 * CG * centroid * volume + CG**2 * volume
    inertia = BODY_DENSITY * (np.sum(about_cg) - about_cg)
    return dict(volume=volume, centroid=centroid, displacement=displacement,
                buoyancy_centroid=buoyancy_centroid, waterplane_area=waterplane_area,
                inertia_diagonal=inertia)
