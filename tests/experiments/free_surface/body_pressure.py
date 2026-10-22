"""Bernoulli pressure and instantaneous six-motion reaction on a wet hull.

The caller owns the fluid solve, moving contact line and ALE time derivative.
Faces wind out of the fluid, into the body. No hydrostatic or coefficient load
is added separately: pressure below includes gravity and acceleration once.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PressureResult:
    acceleration: np.ndarray
    hydrodynamic_load: np.ndarray
    added_mass: np.ndarray
    known_pressure_load: np.ndarray
    diagnostics: dict


def _finite(value, shape, name):
    value = np.asarray(value, dtype=float)
    if value.shape != shape or not np.all(np.isfinite(value)):
        raise ValueError(f'{name} must be finite with shape {shape}')
    return value


def quadrature():
    """Interior degree-three rule shared with material-point field sampling."""
    return (np.array([[1/3, 1/3, 1/3], [.6, .2, .2], [.2, .6, .2], [.2, .2, .6]]),
            np.array([-27., 25., 25., 25.])/48)


def evaluate_body_pressure(points, faces, body_faces, center, mass, inertia_world,
                          velocity_world, potential, potential_rate_ale,
                          mesh_velocity, acceleration_potentials, *, density=1000.,
                          gravity=9.81, external_load=None):
    """Solve current world acceleration, ordered translation then rotation.

    ``potential_rate_ale`` follows the supplied moving nodes while holding the
    current world twist fixed. It may instead be supplied directly at the four
    interior body quadrature samples, shape (body-face count, 4), following their
    interpolated mesh velocity. This permits matched material sampling across
    remeshing without differencing intermittently dry waterline vertices.
    ``acceleration_potentials[:, j]`` is the fluid
    potential response to unit world velocity j at the current geometry and
    fixed surface potential. Thus phi_t = rate_ALE - w_mesh dot grad(phi)
    + acceleration_potentials @ acceleration. The caller must differentiate
    field solves on matching ALE labels and check its directional step size.

    Interior body velocities are recovered from the P1 tangential potential
    gradient and exact rigid-body normal velocity, not a boundary kernel limit.
    Four-point degree-three triangle quadrature integrates the resulting P1
    pressure/load products exactly. This does not certify fluid spatial accuracy.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError('points must be finite with shape (N, 3)')
    count = len(points)
    faces = np.asarray(faces)
    if (faces.ndim != 2 or faces.shape[1] != 3 or faces.dtype.kind not in 'iu'
            or np.any(faces < 0) or np.any(faces >= count)):
        raise ValueError('faces must contain valid triangle vertex indices')
    body_faces = np.asarray(body_faces)
    if (body_faces.ndim != 1 or body_faces.dtype.kind not in 'iu'
            or len(body_faces) == 0 or len(np.unique(body_faces)) != len(body_faces)
            or np.any(body_faces < 0) or np.any(body_faces >= len(faces))):
        raise ValueError('body_faces must contain unique wetted face indices')
    center = _finite(center, (3,), 'center')
    inertia_world = _finite(inertia_world, (3, 3), 'inertia_world')
    velocity_world = _finite(velocity_world, (6,), 'velocity_world')
    potential = _finite(potential, (count,), 'potential')
    rate_shape = np.shape(potential_rate_ale)
    if rate_shape not in ((count,), (len(body_faces), 4)):
        raise ValueError('potential_rate_ale must be nodal or body quadrature values')
    potential_rate_ale = _finite(potential_rate_ale, rate_shape, 'potential_rate_ale')
    mesh_velocity = _finite(mesh_velocity, (count, 3), 'mesh_velocity')
    acceleration_potentials = _finite(
        acceleration_potentials, (count, 6), 'acceleration_potentials')
    if (not np.isfinite(mass) or mass <= 0 or not np.isfinite(density) or density <= 0
            or not np.isfinite(gravity) or gravity < 0):
        raise ValueError('mass/density must be positive and gravity nonnegative')
    if not np.allclose(inertia_world, inertia_world.T, rtol=1e-12,
                       atol=1e-14*np.max(np.abs(inertia_world))):
        raise ValueError('inertia_world must be symmetric')
    np.linalg.cholesky(inertia_world)
    external_load = (np.zeros(6) if external_load is None else
                     _finite(external_load, (6,), 'external_load'))

    body = faces[body_faces]
    triangles = points[body]
    normal_area = np.cross(triangles[:, 1]-triangles[:, 0],
                           triangles[:, 2]-triangles[:, 0])
    twice_area = np.linalg.norm(normal_area, axis=1)
    if np.any(twice_area <= 0):
        raise ValueError('wetted pressure faces must be nondegenerate')
    normals = normal_area/twice_area[:, None]
    opposite = np.roll(triangles, -2, axis=1)-np.roll(triangles, -1, axis=1)
    gradient = np.cross(normals[:, None], opposite)/twice_area[:, None, None]
    tangent_velocity = np.einsum('fcj,fc->fj', gradient, potential[body])

    bary, weights = quadrature()
    samples = np.einsum('qc,fcj->fqj', bary, triangles)
    displacement = samples-center
    omega = velocity_world[3:]
    body_velocity = velocity_world[:3]+np.cross(omega, displacement)
    normal_velocity = np.einsum('fqj,fj->fq', body_velocity, normals)
    fluid_velocity = (tangent_velocity[:, None]
                      + normal_velocity[:, :, None]*normals[:, None])
    moving_velocity = np.einsum('qc,fcj->fqj', bary, mesh_velocity[body])
    material_rate = (np.einsum('qc,fc->fq', bary, potential_rate_ale[body])
                     if potential_rate_ale.ndim == 1 else potential_rate_ale)
    eulerian_rate = material_rate-np.sum(moving_velocity*fluid_velocity, axis=2)
    pressure = -density*(eulerian_rate+.5*np.sum(fluid_velocity**2, axis=2)
                         + gravity*samples[:, :, 2])
    modes = np.concatenate((np.broadcast_to(normals[:, None], samples.shape),
                            np.cross(displacement, normals[:, None])), axis=2)
    measure = twice_area[:, None]/2*weights
    known_load = np.einsum('fq,fq,fqj->j', measure, pressure, modes)
    auxiliary = np.einsum('qc,fck->fqk', bary, acceleration_potentials[body])
    added_mass = density*np.einsum('fq,fqj,fqk->jk', measure, modes, auxiliary)

    rigid_mass = np.zeros((6, 6))
    rigid_mass[:3, :3] = mass*np.eye(3)
    rigid_mass[3:, 3:] = inertia_world
    weight = np.array([0., 0., -mass*gravity, 0., 0., 0.])
    gyro = np.r_[np.zeros(3), np.cross(omega, inertia_world@omega)]
    right = known_load+weight+external_load-gyro
    acceleration = np.linalg.solve(rigid_mass+added_mass, right)
    hydro_load = known_load-added_mass@acceleration
    if not np.all(np.isfinite(acceleration)) or not np.all(np.isfinite(hydro_load)):
        raise FloatingPointError('nonfinite pressure/rigid-body acceleration')
    balance = rigid_mass@acceleration+gyro-hydro_load-weight-external_load
    scale = np.linalg.norm(added_mass)
    diagnostics = dict(
        rigid_body_balance_world=balance.tolist(),
        added_mass_reciprocity_relative=(float(np.linalg.norm(added_mass-added_mass.T)/scale)
                                        if scale else 0.),
        gyro_torque_world=gyro[3:].tolist(),
        wetted_area=float(twice_area.sum()/2),
    )
    return PressureResult(acceleration, hydro_load, added_mass, known_load, diagnostics)
