"""Shared mean geometry and physical incoming train for the WEC-Sim workflow.

The packet is an initial spatial window propagated with signed-kappa dispersion,
not a time-dependent multiplier on body force. Its definition is independent of
the requested stopping time. The profile and bounding method reproduce the
accepted REF002 construction; no amplitude or phase is fitted to the response.
"""
from collections import deque
from dataclasses import dataclass
import itertools
import math
import time

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import root
from scipy.spatial.transform import Rotation

from .constants import GRAVITY, WATER_DENSITY
from .frames import rotation_matrix
from .hydrostatics import equilibrate_heave
from .inputs import InitialState, RigidBody, WaveComponents
from .sea_states import validate_envelope


TIME_STEP_S = .02
RADIATION_MEMORY_S = 10.
CALM_FRACTION = .001


@dataclass(frozen=True)
class PreparedCase:
    body: RigidBody
    equilibrium: InitialState
    hydrostatic: object
    mean_hull_vertices_m: np.ndarray
    physical_faces: np.ndarray
    wet_triangles_m: np.ndarray
    mean_cg_m: np.ndarray
    mass_matrix: np.ndarray
    world_from_native: np.ndarray
    native_waves: WaveComponents
    time_s: np.ndarray
    train: object


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


def prepare_case(body, waves, duration_s):
    """Construct native mean-boat geometry and a single global embedding."""
    duration_s = float(duration_s)
    if not np.isfinite(duration_s) or duration_s <= 0:
        raise ValueError('duration must be positive and finite')
    steps = int(round(duration_s / TIME_STEP_S))
    if steps < 1 or abs(steps * TIME_STEP_S - duration_s) > 1e-9:
        raise ValueError('WEC-Sim duration must be a positive multiple of the fixed 0.02 s time step')
    beam = float(np.ptp(body.mesh.vertices[:, 1]))
    validate_envelope(waves, beam_m=beam, gravity=GRAVITY)
    equilibrium, h = floating_equilibrium(body)
    energy = waves.amplitudes ** 2
    vector = np.array([energy @ np.cos(waves.directions), energy @ np.sin(waves.directions)])
    if energy.sum() == 0:
        yaw = 0.
    elif np.linalg.norm(vector) <= 64 * np.finfo(float).eps * energy.sum():
        raise ValueError('heading: sea has no resolved dominant incoming direction')
    else:
        yaw = float((np.arctan2(vector[1], vector[0]) + np.pi) % (2 * np.pi))
    embedding = Rotation.from_euler('z', yaw).as_matrix()
    native = WaveComponents(waves.amplitudes, waves.angular_frequencies,
                            (waves.directions - yaw) % (2 * np.pi), waves.phases)
    mean_rotation = rotation_matrix(equilibrium.orientation)
    vertices = (body.mesh.vertices - body.cg) @ mean_rotation.T + equilibrium.position
    mass = np.zeros((6, 6))
    mass[:3, :3] = body.mass * np.eye(3)
    mass[3:, 3:] = mean_rotation @ body.inertia_cg @ mean_rotation.T
    train = IncomingTrain(native, beam)
    return PreparedCase(body, equilibrium, h, vertices, body.mesh.faces[body.mesh.physical_faces],
                        h.wetted_triangles, equilibrium.position.copy(), mass, embedding,
                        native, np.arange(steps + 1) * TIME_STEP_S, train)


def spatial_window(tau):
    q = np.clip((abs(np.asarray(tau)) - 15.) / 10., 0., 1.)
    return 1 - 10 * q**3 + 15 * q**4 - 6 * q**5


def _window_spectrum(nu):
    x, w = np.polynomial.legendre.leggauss(256)
    tau = 20 + 5*x
    weights = 5*w*spatial_window(tau)/np.pi
    flat = np.asarray(nu).reshape(-1)
    result = np.empty_like(flat)
    for start in range(0, len(flat), 2048):
        part = flat[start:start+2048]
        result[start:start+2048] = 15/np.pi*np.sinc(15*part/np.pi) + np.cos(part[:, None]*tau) @ weights
    return result.reshape(np.shape(nu))


class IncomingTrain:
    """Accepted fixed40s/15s/25s incoming packet, with arbitrary carrier count.

    C[carrier,sign,frequency] multiplies exp(i*signed_k*s-i*omega*t).
    Frequencies, carrier amplitudes and phases remain unchanged. The finite
    packet has sidebands and is not the original stationary measured spectrum.
    """
    def __init__(self, waves, beam_m):
        maximum = 2*np.pi*min(.6, np.sqrt(GRAVITY/(4*np.pi*beam_m)))
        self.omega = np.r_[np.arange(0., maximum, .005), maximum]
        dw = np.diff(self.omega)
        weights = np.r_[dw[0]/2, (dw[:-1]+dw[1:])/2, dw[-1]/2]
        self.k = self.omega**2 / GRAVITY
        self.theta = waves.directions
        self.carriers = waves
        group = GRAVITY / (2*waves.angular_frequencies)
        k0 = waves.angular_frequencies**2 / GRAVITY
        nu = group[:, None, None] * (np.array([1., -1.])[None, :, None]*self.k - k0[:, None, None])
        jacobian = 2*group[:, None]*self.omega[None, :]/GRAVITY*weights
        self.C = (waves.amplitudes[:, None, None] * _window_spectrum(nu)
                  * np.exp(1j*(waves.phases[:, None, None] + 40*nu)) * jacobian[:, None, :])

    def elevation(self, xy, times):
        xy, times = np.atleast_2d(xy), np.atleast_1d(times)
        value = np.zeros((len(times), len(xy)))
        temporal = np.exp(-1j*times[:, None]*self.omega)
        for j in np.flatnonzero(self.carriers.amplitudes):
            s = xy @ np.array([np.cos(self.theta[j]), np.sin(self.theta[j])])
            spatial = (self.C[j, 0, :, None]*np.exp(1j*self.k[:, None]*s)
                       + self.C[j, 1, :, None]*np.exp(-1j*self.k[:, None]*s))
            value += (temporal @ spatial).real
        return value

    def excitation(self, omega, headings, transfer, stiffness, times):
        """Preserve accepted periodic interpolation and declared low-omega bridge.

        transfer has shape (frequency,heading,6), in dimensional native units
        with exp(-i*omega*t). The exact constant-pressure limit is K[:,Heave].
        This force approximation is not a prescription for surface fields.
        """
        omega, headings, transfer = np.asarray(omega), np.asarray(headings), np.asarray(transfer)
        if transfer.shape != (len(omega), len(headings), 6) or not np.isfinite(transfer).all():
            raise ValueError('excitation: incomplete frequency/heading transfer table')
        if abs(headings[0]) > 1e-12 or abs(headings[-1] - 2*np.pi) > 1e-12:
            raise ValueError('excitation: transfer headings must include the periodic seam')
        if not np.allclose(transfer[:, 0], transfer[:, -1], rtol=0, atol=1e-10):
            raise ValueError('excitation: periodic heading seam does not close')
        if abs(omega[0] - .1) > 1e-12 or omega[-1] < self.omega[-1]:
            raise ValueError('excitation: retained packet is outside the supported transfer table')
        interp = RegularGridInterpolator((omega, headings), transfer, bounds_error=True)
        w = np.broadcast_to(self.omega, self.C.shape)
        theta = np.mod(self.theta[:, None] + [0., np.pi], 2*np.pi)
        theta = np.broadcast_to(theta[:, :, None], self.C.shape)
        h = interp(np.stack([np.maximum(w, .1), theta], axis=-1))
        low = self.omega < .1
        limit = np.asarray(stiffness)[:, 2]
        h[:, :, low, :] = limit + (self.omega[low][None, None, :, None]/.1)**2*(h[:, :, low, :] - limit)
        spectrum = np.sum(self.C[:, :, :, None]*h, axis=(0, 1))
        force = (np.exp(-1j*np.asarray(times)[:, None]*self.omega) @ spectrum).real
        carrier_h = interp(np.c_[self.carriers.angular_frequencies, np.mod(self.theta, 2*np.pi)])
        scale = np.sum(self.carriers.amplitudes[:, None]*abs(carrier_h), axis=0)
        return force, scale

    def certify_startup(self, setup, force, carrier_scale, *, maximum_seconds=480.):
        """Full-hull Taylor bounds and engineering calm loads; fail inconclusive.

        These are startup conditions, not a proof of motion accuracy or linearity.
        They use exactly the immutable field later used for forcing and export.
        """
        amplitude = float(self.carriers.amplitudes.sum())
        if amplitude == 0:
            return
        deadline = time.monotonic() + maximum_seconds
        field = _FieldBounds(self)
        lo = setup.mean_hull_vertices_m.min(axis=0)
        hi = setup.mean_hull_vertices_m.max(axis=0)
        eta_limit = CALM_FRACTION*amplitude
        speed_limit = CALM_FRACTION*float(self.carriers.amplitudes @ self.carriers.angular_frequencies)
        field.certify(lo[:2], hi[:2], [1., 1.], eta_limit, deadline)
        field.certify(np.r_[lo[:2], lo[2]], np.r_[hi[:2], 0.], [1., 1., .2], speed_limit, deadline)
        first = np.asarray(force)[0]
        weight = setup.body.mass*GRAVITY
        length = float(np.ptp(setup.body.mesh.vertices, axis=0).max())
        if (np.linalg.norm(first[:3]) > CALM_FRACTION*weight
                or np.linalg.norm(first[3:]) > CALM_FRACTION*weight*length
                or np.any(abs(first) > CALM_FRACTION*np.asarray(carrier_scale))):
            raise ValueError('startup: initial six-component excitation exceeds the engineering calm limits')
        times = setup.time_s
        peak = np.zeros(len(times))
        temporal = np.exp(-1j*self.omega[:, None]*times)
        for j in np.flatnonzero(self.carriers.amplitudes):
            value = (self.C[j].sum(axis=0) @ temporal)
            peak += abs(value)**2
        arrival = np.sqrt(peak/np.sum(self.carriers.amplitudes**2)).max()
        if arrival < .25:
            raise ValueError('startup: incoming waves do not reach 25% envelope arrival during the requested duration')


class _FieldBounds:
    """Order-six combined-field box bounds, preserving cancellation at centers."""
    def __init__(self, train):
        values = train.C.reshape(-1)
        direction = np.broadcast_to(train.theta[:, None, None], train.C.shape).reshape(-1)
        signs = np.broadcast_to(np.array([1., -1.])[None, :, None], train.C.shape).reshape(-1)
        wave = np.broadcast_to(train.k, train.C.shape).reshape(-1)
        omega = np.broadcast_to(train.omega, train.C.shape).reshape(-1)
        active = abs(values) > 0
        self.C = values[active]
        ka, omega = wave[active], omega[active]
        kx, ky = signs[active]*ka*np.cos(direction[active]), signs[active]*ka*np.sin(direction[active])
        self.k = np.c_[kx, ky, ka]
        self.velocity = np.c_[omega*kx/ka, omega*ky/ka, -1j*omega]
        self.powers = [np.array([v**i for i in range(7)]) for v in (1j*kx, 1j*ky, ka)]
        self.weight_eta = abs(self.C)
        self.weight_velocity = np.sqrt(2)*omega*abs(self.C)

    def certify(self, lo, hi, width, limit, deadline):
        dimension = len(lo)
        axes = [np.linspace(a, b, int(np.ceil((b-a)/w))+1) for a, b, w in zip(lo, hi, width)]
        boxes = [np.array([[axes[d][i] for d, i in enumerate(index)],
                           [axes[d][i+1] for d, i in enumerate(index)]])
                 for index in itertools.product(*(range(len(axis)-1) for axis in axes))]
        pending = deque((box, 0) for box in boxes)
        while pending:
            group = [pending.popleft() for _ in range(min(8, len(pending)))]
            bounds, centers = self._evaluate([box for box, _ in group], dimension, deadline)
            for (box, level), bound, value in zip(group, bounds, centers):
                if value > limit:
                    raise ValueError('startup: initial field exceeds the full-hull calm limit')
                if bound <= limit:
                    continue
                if level >= 3:
                    raise ValueError('startup: full-hull calm bound remains unresolved at the assessment refinement limit')
                middle = box.mean(axis=0)
                for sides in itertools.product([0, 1], repeat=dimension):
                    child = np.array([[box[0, d] if s == 0 else middle[d] for d, s in enumerate(sides)],
                                      [middle[d] if s == 0 else box[1, d] for d, s in enumerate(sides)]])
                    pending.append((child, level+1))

    def _evaluate(self, boxes, dimension, deadline):
        boxes = np.asarray(boxes)
        center = boxes.mean(axis=1)
        half = (boxes[:, 1]-boxes[:, 0])/2
        phase = self.C*np.exp(1j*(center @ self.k[:, :2].T)) if dimension == 2 else self.C*np.exp(
            1j*(center[:, :2] @ self.k[:, :2].T) + center[:, 2, None]*self.k[:, 2])
        weights = self.weight_eta if dimension == 2 else self.weight_velocity
        total = np.zeros(len(boxes))
        center_value = None
        for alpha in itertools.product(range(7), repeat=dimension):
            if sum(alpha) > 6:
                continue
            if time.monotonic() > deadline:
                raise TimeoutError('startup: full-hull calm assessment exceeded its time budget')
            factor = np.ones(len(self.C), dtype=complex)
            width = np.ones(len(boxes))
            for axis, power in enumerate(alpha):
                factor *= self.powers[axis][power]
                width *= half[:, axis]**power/math.factorial(power)
            value = abs((phase @ factor).real) if dimension == 2 else np.linalg.norm(
                (phase @ (factor[:, None]*self.velocity)).real, axis=1)
            if sum(alpha) == 0:
                center_value = value.copy()
            total += value*width
        radius = half @ abs(self.k[:, :dimension]).T
        total += radius**7 @ weights / math.factorial(7)
        return total, center_value
