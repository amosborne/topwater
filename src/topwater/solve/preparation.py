"""Shared mean geometry, packet construction and physical incoming wave train.

Packets are constructed in time, decomposed into physical harmonics and then
propagated with dispersion. Both motion backends consume this same physical
train; no time-dependent force multiplier or fitted response is introduced.
"""
from collections import deque
from dataclasses import dataclass
import itertools
import math
import sys
import time

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial.transform import Rotation

from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.geometry.transforms import rotation_matrix
from topwater.boat.hydrostatics import floating_equilibrium
from topwater.inputs import InitialState, RigidBody, WaveComponents
from topwater.solve.packets import TIME_STEP_S, minimum_cutoff_packet


CALM_FRACTION = .001
from topwater.solve.settings import DEFAULTS

PACKET_ENVELOPE_TOLERANCE_M = DEFAULTS.packet_tolerance_m


def waves_from_sea(sea):
    return WaveComponents(sea.amplitudes_m, sea.angular_frequencies_rad_s,
                          sea.directions_rad, sea.phases_rad)


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
    world_vertical_offset_m: float = 0.


def prepare_case(body, waves, dwell_s, *, rest_s=10., ramp_s=5., time_step_s=TIME_STEP_S,
                 envelope_tolerance_m=PACKET_ENVELOPE_TOLERANCE_M):
    """Construct native mean-boat geometry and a single global embedding."""
    if not np.isfinite(time_step_s) or time_step_s <= 0:
        raise ValueError('time_step_s must be positive and finite')
    timing = np.asarray([dwell_s, rest_s, ramp_s], dtype=float)
    if (not np.isfinite(timing).all() or dwell_s <= 0 or rest_s < 0 or ramp_s <= 0
            or not np.allclose(timing/time_step_s, np.round(timing/time_step_s), rtol=0, atol=1e-8)):
        raise ValueError(f'dwell/ramp must be positive and rest nonnegative, aligned to {time_step_s:g} s')
    duration_s = float(dwell_s+2*rest_s+2*ramp_s)
    steps = int(round(duration_s / time_step_s))
    if steps < 1 or abs(steps * time_step_s - duration_s) > 1e-9:
        raise ValueError('Duration must be a positive multiple of time_step_s')
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
    train = IncomingTrain(native, dwell_s=dwell_s, rest_s=rest_s, ramp_s=ramp_s,
                          time_step_s=time_step_s, envelope_tolerance_m=envelope_tolerance_m)
    return PreparedCase(body, equilibrium, h, vertices, body.mesh.faces[body.mesh.physical_faces],
                        h.wetted_triangles, equilibrium.position.copy(), mass, embedding,
                        native, np.arange(steps + 1) * time_step_s, train, train.dc_mean_m)


class IncomingTrain:
    """Quintic packets, with a separate static sea-level reference offset.

    C[carrier,sign,frequency] multiplies exp(i*signed_k*s-i*omega*t).
    Carrier phases are directly at dwell midpoint in the physical convention
    cos(k.s-omega*(t-midpoint)+phase). DC shifts the native mean-water reference;
    it is never silently removed from the physical incoming surface.
    """
    def __init__(self, waves, *, dwell_s=100., rest_s=10., ramp_s=5.,
                 envelope_tolerance_m=PACKET_ENVELOPE_TOLERANCE_M, time_step_s=TIME_STEP_S):
        self.time_step_s = time_step_s
        self.period_s = float(dwell_s+2*rest_s+2*ramp_s)
        self.midpoint_s = self.period_s/2
        self.rest_s, self.ramp_s, self.dwell_s = float(rest_s), float(ramp_s), float(dwell_s)
        self.envelope_tolerance_m = float(envelope_tolerance_m)
        if not np.isfinite(envelope_tolerance_m) or envelope_tolerance_m <= 0:
            raise ValueError('Packet envelope tolerance must be positive finite metres')
        amplitude = float(waves.amplitudes.sum())
        rows, dc, diagnostics = [], [], []
        from topwater.solve.packet_basis import PacketBasis
        harmonics = PacketBasis()
        maximum_bin = 2
        for j, (a, w, phi) in enumerate(zip(waves.amplitudes, waves.angular_frequencies, waves.phases)):
            if a == 0:
                rows.append(np.zeros(0, complex)); dc.append(0.)
                continue
            tolerance = envelope_tolerance_m*float(a)/amplitude
            packet = minimum_cutoff_packet(float(w/(2*np.pi)), float(-phi), rest_s, ramp_s,
                dwell_s, amplitude_m=float(a), error_tolerance_m=tolerance, time_step_s=time_step_s, _harmonics=harmonics)
            count = packet.positive_components
            positive = np.conj(2*packet.coefficients[1:count+1])
            if len(packet.target.time_s) % 2 == 0 and count == len(packet.target.time_s)//2:
                positive[-1] *= .5  # Real-FFT Nyquist is a singleton, like DC.
            rows.append(positive)
            dc.append(packet.dc_m)
            maximum_bin = max(maximum_bin, count)
            diagnostics.append(dict(carrier=j, cutoff_hz=packet.cutoff_hz, positive_bins=count,
                allocated_tolerance_m=tolerance,
                certified_exceedance_m=packet.diagnostics['certified_envelope_exceedance_bound_m']))
        self.omega = 2*np.pi*np.arange(maximum_bin+1)/self.period_s
        self.k = self.omega**2 / GRAVITY
        self.theta = waves.directions
        self.carriers = waves
        self.C = np.zeros((len(rows), 2, len(self.omega)), complex)
        for j, values in enumerate(rows):
            self.C[j, 0, 1:len(values)+1] = values
        self.dc_by_carrier_m = np.asarray(dc, dtype=float)
        self.dc_mean_m = float(self.dc_by_carrier_m.sum())
        self.packet_diagnostics = diagnostics
        certified = sum(item['certified_exceedance_m'] for item in diagnostics)
        self.certified_origin_exceedance_m = float(certified)
        print(f'PACKET rest_s={rest_s:g} ramp_s={ramp_s:g} dwell_s={dwell_s:g} '
              f'total_s={self.period_s:g} midpoint_s={self.midpoint_s:g} dt_s={self.time_step_s:g} '
              f'positive_frequency_nodes={maximum_bin} cutoff_hz={maximum_bin/self.period_s:.9g} '
              f'active_carriers={len(diagnostics)} allocated_total_tolerance_m={envelope_tolerance_m:g} '
              f'certified_origin_exceedance_m={certified:.9g} dc_mean_m={self.dc_mean_m:.9g}',
              file=sys.stderr)

    def elevation(self, xy, times):
        xy, times = np.atleast_2d(xy), np.atleast_1d(times)
        temporal = np.exp(-1j*times[:, None]*self.omega)
        spatial = np.zeros((len(self.omega), len(xy)), dtype=complex)
        for j in np.flatnonzero(self.carriers.amplitudes):
            s = xy @ np.array([np.cos(self.theta[j]), np.sin(self.theta[j])])
            spatial += (self.C[j, 0, :, None]*np.exp(1j*self.k[:, None]*s)
                        + self.C[j, 1, :, None]*np.exp(-1j*self.k[:, None]*s))
        return (temporal @ spatial).real+float(getattr(self, 'dc_mean_m', 0.))

    def excitation(self, omega, headings, transfer, stiffness, times):
        """OFFLINE BASELINE ORACLE ONLY; never called by the production solve.

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

    def directional_elevation_records(self, times):
        """Physical elevations at native origin, keeping both signed branches.

        No hydrodynamic transfer or force calculation occurs here. Directions
        are propagation radians; absolute time and the original weighted C
        coefficients are retained. A genuinely zero sea has one zero record.
        """
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not len(times) or not np.isfinite(times).all():
            raise ValueError('directional elevation requires finite sample times')
        active = np.flatnonzero(self.carriers.amplitudes)
        if not len(active):
            return np.zeros(1), np.zeros((len(times), 1))
        directions = np.mod(self.theta[active, None] + [0., np.pi], 2*np.pi).reshape(-1)
        coefficients = self.C[active].reshape(-1, len(self.omega))
        elevation = (np.exp(-1j*times[:, None]*self.omega) @ coefficients.T).real
        elevation[:, ::2] += getattr(self, 'dc_by_carrier_m', np.zeros(len(self.theta)))[active]
        return directions, elevation

    def certify_incident_startup(self, setup, *, maximum_seconds=480.):
        """Full-hull incident eta/u and arrival bounds; fail inconclusive.

        These are startup conditions, not a proof of motion accuracy or linearity.
        The native WEC preparation separately certifies its own initial force.
        """
        amplitude = float(self.carriers.amplitudes.sum())
        if amplitude == 0:
            return
        deadline = time.monotonic() + maximum_seconds
        field = _FieldBounds(self)
        lo = setup.mean_hull_vertices_m.min(axis=0)
        hi = setup.mean_hull_vertices_m.max(axis=0)
        eta_limit = self.envelope_tolerance_m
        speed_limit = CALM_FRACTION*float(self.carriers.amplitudes @ self.carriers.angular_frequencies)
        field.certify(lo[:2], hi[:2], [1., 1.], eta_limit, deadline)
        field.certify(np.r_[lo[:2], lo[2]], np.r_[hi[:2], 0.], [1., 1., .2], speed_limit, deadline)
        times = setup.time_s
        peak = np.zeros(len(times))
        temporal = np.exp(-1j*self.omega[:, None]*times)
        for j in np.flatnonzero(self.carriers.amplitudes):
            value = (self.C[j].sum(axis=0) @ temporal)+self.dc_by_carrier_m[j]
            peak += abs(value)**2
        arrival = np.sqrt(peak/np.sum(self.carriers.amplitudes**2)).max()
        if arrival < .25:
            raise ValueError('startup: incoming waves do not reach 25% envelope arrival during the requested duration')


class _FieldBounds:
    """Order-six combined-field box bounds, preserving cancellation at centers."""
    def __init__(self, train):
        self.dc_mean_m = float(getattr(train, 'dc_mean_m', 0.))
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
            value = abs((phase @ factor).real+(self.dc_mean_m if sum(alpha) == 0 else 0.)) if dimension == 2 else np.linalg.norm(
                (phase @ (factor[:, None]*self.velocity)).real, axis=1)
            if sum(alpha) == 0:
                center_value = value.copy()
            total += value*width
        radius = half @ abs(self.k[:, :dimension]).T
        total += radius**7 @ weights / math.factorial(7)
        return total, center_value
