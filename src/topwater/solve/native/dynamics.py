"""Coupled Newton-Euler integration; no MATLAB or private model dependencies."""
import warnings

import numpy as np
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater.solve.native.excitation import Excitation
from topwater.solve.native.influence import BoatState, LoadEvaluation, finite_vector, readonly
from topwater.solve.native.radiation import history_force, radiation_kernel, stage_weights
from topwater.solve.settings import DEFAULTS


from topwater.solve.types import Motion


def effective_inertia(mass, added):
    """Continuous-time form of the current reference Newton-Euler equations.

The reference solid includes a scalar translational augmentation and the full
rotational added inertia. Keep its orientation/gyroscopic terms unchanged, but
move the residual acceleration load to the LHS instead of delaying it. At the
mean attitude, solid + residual is exactly M + Ainf.
"""
    solid = np.array(mass, copy=True)
    solid[:3, :3] += 2*np.trace(added[:3, :3])*np.eye(3)
    rotational = np.triu(added[3:, 3:])
    # Simscape's solid has symmetric products; retain genuine BEM asymmetry
    # in the residual rather than silently symmetrizing the complete operator.
    solid[3:, 3:] += rotational+np.triu(rotational, 1).T
    return solid, mass+added-solid


def quaternion_derivative(quaternion, angular_velocity):
    w, xyz = quaternion[0], quaternion[1:]
    # Left multiplication: angular velocity is resolved in the parent frame.
    return .5*np.r_[-angular_velocity@xyz,
                    w*angular_velocity+np.cross(angular_velocity, xyz)]


def integrate_motion(setup, hydro, *, influence=None, dt_s=DEFAULTS.time_step_s, excitation=None, kernel=None):
    """Integrate motion and optional device states on a uniform accepted clock.

dt_s is exposed for convergence checks and must exactly subdivide the case.
Radiation uses the same internal clock. The default workflow step is 10 ms.
"""
    end = float(setup.time_s[-1])
    if not np.isfinite(dt_s) or dt_s <= 0 or abs(round(end/dt_s)*dt_s-end) > 1e-9:
        raise ValueError('Native timestep must be finite, positive, and divide the case duration')
    time = np.arange(round(end/dt_s)+1)*dt_s
    excitation = Excitation.from_case(setup, hydro) if excitation is None else excitation
    kernel = radiation_kernel(hydro.radiation_omega, hydro.radiation_damping,
        dt_s=dt_s, duration_s=end)[1] if kernel is None else kernel
    operators = [stage_weights(kernel, dt_s, fraction) for fraction in (0., .5, 1.)]
    solid, residual = effective_inertia(hydro.M, hydro.Ainf)
    heading = np.asarray(setup.world_from_native)
    equilibrium = Rotation.from_quat(setup.equilibrium.orientation, scalar_first=True).as_matrix()
    offset = np.array([0., 0., setup.world_vertical_offset_m])
    state = np.zeros(13)
    state[2], state[3] = -setup.world_vertical_offset_m, 1.

    def snapshot(stamp, value):
        perturbation = Rotation.from_quat(value[3:7], scalar_first=True).as_matrix()
        return BoatState(float(stamp), readonly(setup.mean_cg_m+heading@value[:3]+offset),
                         readonly(Rotation.from_matrix(heading@perturbation@equilibrium).as_quat(scalar_first=True)),
                         readonly(heading@value[7:10]), readonly(heading@value[10:13]))

    def call(method, stamp, value, *args):
        try:
            return method(snapshot(stamp, value), *args)
        except Exception as error:
            name = type(influence).__qualname__
            raise RuntimeError(f'External model {name}.{method.__name__} failed at t={stamp:.9g} s: {error}') from error

    if influence is not None:
        raw = call(influence.initial_state, 0., state)
        try:
            raw = np.asarray(raw)
            if raw.ndim != 1:
                raise ValueError('initial_state must be a numeric vector')
            initial = finite_vector(raw, len(raw), 'initial_state')
        except (ValueError, TypeError) as error:
            raise ValueError(f'External model {type(influence).__qualname__} at t=0 s: {error}') from error
        state = np.r_[state, initial]
    nz = len(state)-13
    count = len(time)
    history = np.zeros((count, 6))
    states = np.zeros((count, len(state)))
    accelerations = np.zeros((count, 6))
    components = {key: np.zeros((count, 6)) for key in
                  ('excitation', 'radiation', 'restoring', 'additional', 'balance_residual')}

    def update(stamp, value):
        callback = getattr(influence, 'on_step', None)
        if callback is not None:
            call(callback, stamp, value, readonly(value[13:]))

    def derivative(stamp, value, past, current):
        rotation = Rotation.from_quat(value[3:7], scalar_first=True)
        R = rotation.as_matrix()
        with warnings.catch_warnings():
            warnings.simplefilter('error', UserWarning)
            angles = rotation.as_euler('xyz')
        q = np.r_[value[:3], angles]
        v, omega = value[7:13], value[10:13]
        radiation = past+current@v
        restoring = hydro.K@q
        force = excitation.at(stamp)
        additional, dz = np.zeros(6), np.zeros(nz)
        if influence is not None:
            load = call(influence.evaluate, stamp, value, readonly(value[13:]))
            try:
                if not isinstance(load, LoadEvaluation):
                    raise TypeError('evaluate must return LoadEvaluation')
                additional = np.r_[heading.T@finite_vector(load.force_world_N, 3, 'force_world_N'),
                                   heading.T@finite_vector(load.torque_world_Nm, 3, 'torque_world_Nm')]
                dz = finite_vector(load.state_derivative, nz, 'state_derivative')
            except (ValueError, TypeError) as error:
                raise ValueError(f'External model {type(influence).__qualname__} at t={stamp:.9g} s: {error}') from error
        inertia = R@solid[3:, 3:]@R.T
        matrix = residual.copy()
        matrix[:3, :3] += solid[:3, :3]
        matrix[3:, 3:] += inertia
        gyroscopic = np.r_[np.zeros(3), np.cross(omega, inertia@omega)]
        rhs = force-radiation-restoring+additional-gyroscopic
        acceleration = np.linalg.solve(matrix, rhs)
        change = np.r_[v[:3], quaternion_derivative(value[3:7], omega), acceleration, dz]
        values = (force, radiation, restoring, additional, matrix@acceleration-rhs)
        return change, values

    update(0., state)
    # One policy boundary, not a Windows threadpool scan at each small solve.
    with threadpool_limits(limits=1):
        for n, stamp in enumerate(time):
            history[n] = state[7:13]
            states[n] = state
            past0 = history_force(history, n, operators[0][0])
            k1, values = derivative(stamp, state, past0, operators[0][1])
            accelerations[n] = k1[7:13]
            for key, value in zip(components, values):
                components[key][n] = value
            if n == count-1:
                break
            half = history_force(history, n, operators[1][0])
            last = history_force(history, n, operators[2][0])
            k2, _ = derivative(stamp+dt_s/2, state+dt_s/2*k1, half, operators[1][1])
            k3, _ = derivative(stamp+dt_s/2, state+dt_s/2*k2, half, operators[1][1])
            k4, _ = derivative(stamp+dt_s, state+dt_s*k3, last, operators[2][1])
            state = state+dt_s/6*(k1+2*k2+2*k3+k4)
            state[3:7] /= np.linalg.norm(state[3:7])
            if not np.isfinite(state).all():
                raise ArithmeticError(f'Native state became nonfinite at t={time[n+1]:.9g} s')
            update(time[n+1], state)
    rotations = Rotation.from_quat(states[:, 3:7], scalar_first=True)
    q = np.c_[states[:, :3], rotations.as_euler('xyz')]
    orientation = Rotation.from_matrix(heading@rotations.as_matrix()@equilibrium).as_quat(scalar_first=True)
    position = states[:, :3]@heading.T+setup.mean_cg_m+offset
    return Motion(time, q, history, accelerations, position, orientation,
                  components['excitation'], states[:, 13:], components)
