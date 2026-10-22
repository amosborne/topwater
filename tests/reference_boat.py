"""Authored 24 ft open V-hull and round loading; not a production boat model.

The aft three quarters are prismatic; the forward quarter tapers to the raised
stem. All facets are planar, so subdivision changes only numerical resolution.
The artificial sheer closure has neither shell mass nor pressure-bearing area.
"""
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation

from topwater import HullMesh, InitialState, RigidBody
from topwater.hydrostatics import equilibrate_heave


LENGTH = 24*.3048
BEAM = 9*.3048
DEPTH = 1.
DEADRISE_DEGREES = 20.
CHINE_BEAM_FRACTION = .8
SHELL_MASS = 1600.
PAYLOAD_MASS = 600.
MOTOR_MASS = 300.
PAYLOAD_POSITION = np.array([0., 0., .8])
MOTOR_POSITION = np.array([-LENGTH/2, 0., .7])
WATER_DENSITY = 1025.
GRAVITY = 9.81


def mesh_arrays(refinement=0):
    """Return vertices, outward triangles, ordered rim and artificial face IDs.

    Coordinates are x forward, y port, z up; keel z=0 and sheer z=1 m.
    Every refinement splits each triangle into four without moving the hull.
    """
    if isinstance(refinement, bool) or not isinstance(refinement, int) or refinement < 0:
        raise ValueError('nonnegative integer refinement required')
    chine = CHINE_BEAM_FRACTION*BEAM/2
    z = chine*np.tan(np.deg2rad(DEADRISE_DEGREES))
    ring = np.array([[0., 0.], [chine, z], [BEAM/2, DEPTH],
                     [-BEAM/2, DEPTH], [-chine, z]])
    vertices = np.vstack([np.column_stack((np.full(5, x), ring))
                          for x in (-LENGTH/2, LENGTH/4)] + [[[LENGTH/2, 0., DEPTH]]])
    faces, cap = [], []
    for j in range(5):
        a, b = j, (j+1) % 5
        faces.extend([[a, b, b+5], [a, b+5, a+5], [a+5, b+5, 10]])
        cap.extend([j == 2]*3)
    # End-view ring is counterclockwise about +x; stern normal is -x.
    faces.extend([[0, j+1, j] for j in (1, 2, 3)])
    cap.extend([False]*3)
    faces, cap = np.asarray(faces, dtype=np.int64), np.asarray(cap)
    for _ in range(refinement):
        edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
        unique, inverse = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True)
        ab, bc, ca = (inverse+len(vertices)).reshape(3, -1)
        a, b, c = faces.T
        faces = np.vstack([np.column_stack(v) for v in
                           ((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
        vertices = np.vstack((vertices, vertices[unique].mean(axis=1)))
        cap = np.tile(cap, 4)
    boundary_edges = np.concatenate([faces[cap][:, [0, 1]], faces[cap][:, [1, 2]],
                                     faces[cap][:, [2, 0]]])
    _, inverse, counts = np.unique(np.sort(boundary_edges, axis=1), axis=0,
                                   return_inverse=True, return_counts=True)
    rim_edges = boundary_edges[counts[inverse] == 1]
    successor = dict(rim_edges)
    rim = [min(successor)]
    while successor[rim[-1]] != rim[0]:
        rim.append(successor[rim[-1]])
    return vertices, faces, np.asarray(rim), np.flatnonzero(cap)


def shell_loading(vertices, physical_faces):
    """Exact uniform-area triangular shell moments, normalized to 1600 kg.

    Input shell mass is an explicit budget, not inferred laminate thickness or
    density. The point loads have no intrinsic inertia; their offsets do.
    """
    triangles = vertices[physical_faces]
    area = np.linalg.norm(np.cross(triangles[:, 1]-triangles[:, 0],
                                   triangles[:, 2]-triangles[:, 0]), axis=1)/2
    dm = SHELL_MASS*area/area.sum()
    sums = triangles.sum(axis=1)
    second = (np.einsum('fi,fj->fij', sums, sums)
              +np.einsum('fai,faj->fij', triangles, triangles))/12
    first = (dm[:, None]*sums/3).sum(axis=0)
    raw_second = (dm[:, None, None]*second).sum(axis=0)
    mass = SHELL_MASS+PAYLOAD_MASS+MOTOR_MASS
    for m, p in ((PAYLOAD_MASS, PAYLOAD_POSITION), (MOTOR_MASS, MOTOR_POSITION)):
        first += m*p
        raw_second += m*np.outer(p, p)
    cg = first/mass
    central = raw_second-mass*np.outer(cg, cg)
    return dict(mass=mass, cg=cg, inertia_cg=np.trace(central)*np.eye(3)-central,
                shell_area=float(area.sum()), shell_cg=(dm[:, None]*sums/3).sum(axis=0)/SHELL_MASS)


def reference_body(refinement=0):
    vertices, faces, rim, cap = mesh_arrays(refinement)
    mesh = HullMesh(vertices, faces, rim, cap)
    loading = shell_loading(mesh.vertices, mesh.faces[mesh.physical_faces])
    return RigidBody(mesh, loading['mass'], loading['cg'], loading['inertia_cg'])


@lru_cache(None)
def equilibrium(refinement=0):
    """Independent flat-water heave and pitch balance of this actual mesh.

    Positive trim is bow up (negative world-y rotation). Port/starboard symmetry
    supplies roll balance; this routine does not assume roll stability.
    """
    body = reference_body(refinement)

    def at(trim):
        orientation = np.roll(Rotation.from_rotvec([0., -trim, 0.]).as_quat(), 1)
        return equilibrate_heave(body, InitialState(orientation=orientation),
                                 density=WATER_DENSITY, gravity=GRAVITY)

    trim = brentq(lambda angle: at(angle)[1].moment[1], -.1, .1, xtol=1e-13)
    return at(trim)
