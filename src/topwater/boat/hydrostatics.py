"""Calm-water geometry and pressure; no BEM or wave evolution is assumed.

Pressure is rho*g*depth on the physical wetted hull. This is an independent
verification/initialization tool, not a second load added to dynamic Bernoulli.
"""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, root
from scipy.spatial.transform import Rotation
from ..constants import GRAVITY, WATER_DENSITY

from ..geometry.clipping import clip_flat_water
from .._validation import array
from ..geometry.transforms import rotation_matrix
from ..inputs import InitialState, RigidBody


class FloodingError(ValueError):
    """The open rim touches or crosses the specified water surface."""


class LossOfContactError(ValueError):
    """The hull has no wetted area enclosing positive displaced volume."""


class NoEquilibriumError(ValueError):
    """There is no surface-piercing heave balance before a failure boundary."""


@dataclass(frozen=True, eq=False)
class HydrostaticResult:
    """World-axis pressure load and its CG reference; triangles remain uncapped.

    stiffness is minus the load derivative under world translations/rotations
    about the current CG. It is symmetric at a force/moment-balanced pose.
    """

    reference_point: np.ndarray
    displaced_volume: float
    center_of_buoyancy: np.ndarray
    waterplane_area: float
    wetted_area: float
    force: np.ndarray
    moment: np.ndarray
    stiffness: np.ndarray
    wetted_triangles: np.ndarray
    source_faces: np.ndarray
    waterline_segments: np.ndarray
    rim_clearance: float | None


def _parameters(body, state, water_level, density, gravity):
    if not isinstance(body, RigidBody) or not isinstance(state, InitialState):
        raise ValueError('expected RigidBody and InitialState')
    level = array(water_level, (), 'water_level').item()
    rho = array(density, (), 'density').item()
    g = array(gravity, (), 'gravity').item()
    if rho <= 0 or g <= 0:
        raise ValueError('density and gravity must be positive')
    relative = (body.mesh.vertices-body.cg) @ rotation_matrix(state.orientation).T
    return relative, level, rho, g


def _waterplane(segments):
    """Exact polygon moments from oriented wet-hull edges, relative to CG x/y."""
    x, y = segments[:, 0, 0], segments[:, 0, 1]
    u, v = segments[:, 1, 0], segments[:, 1, 1]
    cross = -(x*v-u*y)  # Wet-hull boundary is opposite the upward plane closure.
    area = cross.sum()/2
    first = np.array([np.sum((x+u)*cross), np.sum((y+v)*cross)])/6
    second = np.array([[np.sum((y*y+y*v+v*v)*cross)/12,
                        np.sum((2*x*y+x*v+u*y+2*u*v)*cross)/24],
                       [0., np.sum((x*x+x*u+u*u)*cross)/12]])
    second[1, 0] = second[0, 1]
    return area, first, second


def hydrostatics(body, state=InitialState(), *, water_level=0., density=1000., gravity=9.81):
    """Clip the actual hull at flat water and integrate pressure about world CG.

    Rim contact is flooding onset; zero contact is a failed case. Artificial
    closure faces never enter clipping or pressure. Detection here is for flat
    water only; local evolved-surface events belong to the future wave solver.
    """
    relative, level, rho, g = _parameters(body, state, water_level, density, gravity)
    local_level = level-state.position[2]
    mesh = body.mesh
    clearance = None
    if len(mesh.opening_rim):
        clearance = float(np.min(relative[mesh.opening_rim, 2])-local_level)
        if clearance <= 0:
            raise FloodingError(f'opening rim reaches the water surface (clearance {clearance:g} m)')
    wet, source, lines, volume, first, force, moment, area = clip_flat_water(
        relative, mesh.faces, mesh.physical_faces, local_level, True)
    if len(wet) == 0 or volume <= 0:
        raise LossOfContactError('hull has lost water contact')
    if not np.isfinite(volume):
        raise ValueError('non-finite clipped volume')
    buoyancy = first/volume+np.array([0., 0., local_level])
    plane_area, plane_first, plane_second = _waterplane(lines)
    stiffness = np.zeros((6, 6))
    stiffness[2, 2] = plane_area
    stiffness[2, 3], stiffness[3, 2] = plane_first[1], plane_first[1]
    stiffness[2, 4], stiffness[4, 2] = -plane_first[0], -plane_first[0]
    stiffness[3, 3] = plane_second[0, 0]+volume*buoyancy[2]
    stiffness[4, 4] = plane_second[1, 1]+volume*buoyancy[2]
    stiffness[3, 4] = stiffness[4, 3] = -plane_second[0, 1]
    stiffness[3, 5], stiffness[4, 5] = -volume*buoyancy[0], -volume*buoyancy[1]
    force *= rho*g
    moment *= rho*g
    stiffness *= rho*g
    wet += state.position
    lines += state.position
    buoyancy += state.position
    for values in (force, moment, stiffness, wet, source, lines, buoyancy):
        values.setflags(write=False)
    return HydrostaticResult(state.position, volume, buoyancy, float(plane_area), area,
                             force, moment, stiffness, wet, source, lines, clearance)


def equilibrate_heave(body, state=InitialState(), *, water_level=0., density=1000., gravity=9.81):
    """Return (rest_state, hydrostatics) with heave force balance at fixed attitude.

    Horizontal CG and orientation are retained; returned velocities are zero.
    This does not solve roll/pitch equilibrium or assert stability.
    """
    relative, level, rho, g = _parameters(body, state, water_level, density, gravity)
    mesh = body.mesh
    scale = float(np.ptp(relative, axis=0).max())
    margin = 64*np.finfo(float).eps*scale
    if len(mesh.opening_rim):
        lower = -np.min(relative[mesh.opening_rim, 2])+margin
    else:
        lower = -np.max(relative[:, 2])-margin
    upper = -np.min(relative[:, 2])-margin
    target = body.mass/rho
    def volume(height):
        return clip_flat_water(relative, mesh.faces, mesh.physical_faces,
                               -height, False)[3]
    if lower >= upper or volume(lower) <= target:
        raise NoEquilibriumError('no heave balance before rim immersion or complete submersion')
    if volume(upper) >= target:
        raise NoEquilibriumError('heave balance is unresolved near loss of contact')
    height = brentq(lambda height: volume(height)-target, lower, upper,
                    xtol=1e-12*scale, rtol=1e-12)
    position = state.position.copy()
    position[2] = level+height
    rest = InitialState(position, state.orientation)
    return rest, hydrostatics(body, rest, water_level=level, density=rho, gravity=g)


def floating_equilibrium(body):
    """Balance heave/roll/pitch at fixed horizontal CG and heading.

    Heave is bracketed by actual dry/rim-contact limits. The two attitude roots
    use exact calm-water facet integrals; no simulated settling is performed.
    Stability is checked in the heave/roll/pitch subspace at the resulting root.
    """
    length = float(np.ptp(body.mesh.vertices, axis=0).max())
    scale = body.mass * GRAVITY * length

    def state_at(angles):
        orientation = Rotation.from_euler('xyz', [*angles, 0.]).as_quat(scalar_first=True)
        return equilibrate_heave(body, InitialState(orientation=orientation),
                                 density=WATER_DENSITY, gravity=GRAVITY)

    def residual(angles):
        _, h = state_at(angles)
        return h.moment[:2] / scale

    solution = root(residual, np.zeros(2), method='hybr', options={'xtol': 1e-10})
    if not solution.success and np.linalg.norm(residual(solution.x), ord=np.inf) > 1e-10:
        raise ValueError('equilibrium: could not balance roll and pitch before physical failure limits')
    state, h = state_at(solution.x)
    wrench = np.r_[h.force + [0., 0., -body.mass * GRAVITY], h.moment]
    scales = np.r_[np.full(3, body.mass * GRAVITY), np.full(3, scale)]
    if np.max(abs(wrench) / scales) > 1e-9:
        raise ValueError('equilibrium: residual force or moment exceeds the balance tolerance')
    dofs = [2, 3, 4]
    stiffness = h.stiffness[np.ix_(dofs, dofs)]
    units = np.diag([length, 1., 1.])
    scaled = units @ stiffness @ units
    if np.linalg.eigvalsh((scaled + scaled.T) / 2).min() <= 1e-10 * body.mass * GRAVITY * length:
        raise ValueError('equilibrium: heave/roll/pitch equilibrium is not stably restoring')
    return state, h
