"""One canonical numerical hull, shared by tests, timings and figures.

The reference is an open lower half-ellipsoid, not a measured boat.
The opening rim is at body z=0. A marked artificial cap closes the geometry;
it is not an impermeable deck. Loading is a thin shell plus one point payload.
This is test infrastructure, not a production mesh generator or hydrostatics API.
"""

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq


SEMI_AXES = np.array([1.75, 0.45, 0.35])  # 3.5 m long, 0.9 m wide, 0.35 m deep
RIM_HEIGHT = 0.0  # body coordinates, not the loaded waterline
WATER_DENSITY = 1000.0  # kg/m^3 for this reference case only
SHELL_THICKNESS = 0.003  # m; thin-shell area model, no through-thickness inertia
LAMINATE_DENSITY = 1800.0  # kg/m^3; illustrative fiberglass/resin assumption
AREAL_DENSITY = SHELL_THICKNESS * LAMINATE_DENSITY
PAYLOAD_MASS = 80.0  # kg; a single adult represented as a point mass
PAYLOAD_POSITION = np.array([0.0, 0.0, 0.1])  # m; CG 0.1 m above rim at centre
FULL_VOLUME = 2 * np.pi * np.prod(SEMI_AXES) / 3


def _combine_loading(shell_mass, shell_first, shell_second):
    """Combine shell raw moments with the point payload, then shift to total CG."""
    mass = shell_mass + PAYLOAD_MASS
    cg = (shell_first + PAYLOAD_MASS * PAYLOAD_POSITION) / mass
    second = shell_second + PAYLOAD_MASS * np.outer(PAYLOAD_POSITION, PAYLOAD_POSITION)
    inertia_origin = np.trace(second) * np.eye(3) - second
    inertia_cg = inertia_origin - mass * (np.dot(cg, cg) * np.eye(3) - np.outer(cg, cg))
    return dict(shell_mass=shell_mass, shell_cg=shell_first / shell_mass,
                mass=mass, cg=cg, inertia_cg=inertia_cg)


def reference_loading(order=64):
    """Independent parametric surface quadrature; includes no mesh triangles.

    Integrate one azimuth quadrant and reflect its even moments over all four.
    Symmetry makes the other first moments and all mixed second moments zero.
    """
    nodes, weights = leggauss(order)
    angles = (nodes + 1) * np.pi / 4
    theta, phi = np.meshgrid(angles, angles, indexing="ij")
    a, b, c = SEMI_AXES
    st, ct = np.sin(theta), np.cos(theta)
    sp, cp = np.sin(phi), np.cos(phi)
    points = np.stack([a*st*cp, b*st*sp, -c*ct], axis=-1)
    jacobian = st * np.sqrt((b*c*st*cp)**2 + (a*c*st*sp)**2 + (a*b*ct)**2)
    dm = 4 * np.outer(weights, weights) * (np.pi / 4)**2 * jacobian * AREAL_DENSITY
    shell_mass = float(np.sum(dm))
    first = np.array([0.0, 0.0, np.sum(dm * points[:, :, 2])])
    second = np.diag(np.sum(dm[:, :, None] * points**2, axis=(0, 1)))
    return _combine_loading(shell_mass, first, second)


LOADING = reference_loading()
SHELL_MASS = LOADING["shell_mass"]
MASS = LOADING["mass"]
CG = LOADING["cg"]
INERTIA = np.diag(LOADING["inertia_cg"])
DISPLACEMENT = MASS / WATER_DENSITY


def displaced_volume(level):
    """Exact lower-ellipsoid volume below a body-coordinate level in [-c, 0]."""
    a, b, c = SEMI_AXES
    return np.pi * a * b * (level - level**3 / (3*c*c) + 2*c/3)


# An analytical shape-specific draft estimate, not a production buoyancy solve.
WATERLINE_HEIGHT = brentq(
    lambda level: displaced_volume(level) - DISPLACEMENT,
    -SEMI_AXES[2], RIM_HEIGHT, xtol=1e-14,
)


def mesh_arrays(refinement=2):
    """Return contiguous vertices/faces; there are 8*4**refinement panels."""
    if isinstance(refinement, bool) or not isinstance(refinement, int) or refinement < 0:
        raise ValueError("refinement must be a nonnegative integer")
    vertices = np.array([
        [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, -1],
        [0, 0, 0],
    ], dtype=float)
    faces = np.array([
        [2, 0, 4], [1, 2, 4], [3, 1, 4], [0, 3, 4],
        [5, 0, 2], [5, 2, 1], [5, 1, 3], [5, 3, 0],
    ], dtype=np.int64)
    cap_faces = np.all(vertices[faces, 2] == 0, axis=1)
    for _ in range(refinement):
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        unique, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        midpoints = vertices[unique].mean(axis=1)
        curved = np.zeros(len(unique), dtype=bool)
        np.logical_or.at(curved, inverse, np.tile(~cap_faces, 3))
        midpoints[curved] /= np.linalg.norm(midpoints[curved], axis=1)[:, None]
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

    Rim immersion terminates calm-water evaluation before an artificial cap can
    produce pressure loads. Evolving-wave contact detection remains future work.
    """
    on_surface = np.isclose(np.sum((vertices / SEMI_AXES)**2, axis=1), 1, rtol=0, atol=1e-12)
    rim = np.flatnonzero((vertices[:, 2] == RIM_HEIGHT) & on_surface)
    angles = np.arctan2(vertices[rim, 1], vertices[rim, 0])
    cap = np.flatnonzero(np.all(vertices[faces, 2] == RIM_HEIGHT, axis=1))
    return rim[np.argsort(angles)], cap


def input_data(refinement=2):
    """Analytical calm-water reference pose plus explicit synthetic wave inputs.

    Integrated thin-shell mass plus one centred point payload supplies mass, CG
    and inertia. The artificial cap has no mass. The exact ellipsoid's displaced
    volume gives a shape-specific reference draft; a finite mesh only approximates
    that pose. No wave solver exists. Initial velocities are zero. No random
    phases are generated here; the seed is illustrative provenance metadata.
    """
    vertices, faces = mesh_arrays(refinement)
    rim, cap = opening_boundaries(vertices, faces)
    return {
        "format_version": 2,
        "body": {
            "mesh": {"vertices": vertices.tolist(), "faces": faces.tolist(),
                     "opening_rim": rim.tolist(), "artificial_faces": cap.tolist()},
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
            "position": (CG - [0, 0, WATERLINE_HEIGHT]).tolist(),
            "orientation": [1, 0, 0, 0],
            "linear_velocity": [0, 0, 0], "angular_velocity": [0, 0, 0],
        },
    }


def surface_integrals(vertices, faces):
    """Independent divergence-theorem reference integrals of planar triangles.

    Uses z*n_z for volume (the input validator uses tetrahedral determinants).
    First moments use exact triangle polynomial averages.
    """
    triangles = vertices[faces]
    p, q, r = triangles.transpose(1, 0, 2)
    normal_area = np.cross(q - p, r - p) / 2
    volume = np.sum(normal_area[:, 2] * triangles[:, :, 2].mean(axis=1))
    average_squared = (p*p + q*q + r*r + p*q + p*r + q*r) / 6
    first_moments = np.sum(normal_area * average_squared, axis=0) / 2
    return volume, first_moments / volume


def geometry_metrics(vertices, faces):
    """Full-envelope volume/centroid; no submerged clipping is implemented."""
    volume, centroid = surface_integrals(vertices, faces)
    return dict(volume=volume, centroid=centroid)


def mesh_loading(vertices, faces):
    """Integrate thin-shell triangle mass, CG and full inertia; cap has no mass."""
    _, cap = opening_boundaries(vertices, faces)
    triangles = vertices[np.delete(faces, cap, axis=0)]
    p, q, r = triangles.transpose(1, 0, 2)
    dm = AREAL_DENSITY * np.linalg.norm(np.cross(q-p, r-p), axis=1) / 2
    summed = p + q + r
    # Exact average r*r^T on each planar triangle, not a centroid point lump.
    quadratic = (np.einsum("ni,nj->nij", summed, summed)
                 + np.einsum("nki,nkj->nij", triangles, triangles)) / 12
    return _combine_loading(float(np.sum(dm)), np.sum(dm[:, None] * summed / 3, axis=0),
                            np.sum(dm[:, None, None] * quadratic, axis=0))
