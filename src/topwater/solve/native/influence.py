"""Public boundary for private force, actuator and feedback-control models."""
from dataclasses import dataclass, field
from typing import Protocol

import numpy as np


@dataclass(frozen=True)
class BoatState:
    """Read-only SI state: world east/north/up, body bow/port/up.

Orientation is active body-to-world (w,x,y,z); angular velocity is an actual
world-frame vector in rad/s, not a vector of Euler-angle derivatives.
"""
    time_s: float
    position_m: np.ndarray
    orientation_wxyz: np.ndarray
    velocity_m_s: np.ndarray
    angular_velocity_rad_s: np.ndarray


@dataclass(frozen=True)
class LoadEvaluation:
    """World force [N], world torque about current CG [N m], and dz/dt.

Device hardware mass must already be reconciled in the Boat mass properties.
The solver never infers additional mass from an Influence.
"""
    force_world_N: np.ndarray = field(default_factory=lambda: np.zeros(3))
    torque_world_Nm: np.ndarray = field(default_factory=lambda: np.zeros(3))
    state_derivative: np.ndarray = field(default_factory=lambda: np.zeros(0))


class Influence(Protocol):
    """Structural interface; no inheritance or plugin registration required.

evaluate must be pure: RK stages may be provisional or repeated. on_step is
optional and is the ONLY place to advance sampled controller memory/commands.
Both methods receive read-only snapshots. on_step has the same (boat, state)
arguments as evaluate and its return value is ignored. A slower controller can
update on integer multiples of the integration step; commands are held between
updates. Continuous device states must be returned through state_derivative.
Acceleration-dependent loads and algebraic device mass coupling are not part
of this explicit interface.
"""
    def initial_state(self, boat: BoatState) -> np.ndarray: ...
    def evaluate(self, boat: BoatState, state: np.ndarray) -> LoadEvaluation: ...


def readonly(value):
    """Own the snapshot: external code cannot modify the integrator's storage."""
    result = np.array(value, dtype=float, copy=True)
    result.flags.writeable = False
    return result


def finite_vector(value, size, name):
    raw = np.asarray(value)
    if not np.isrealobj(raw):
        raise ValueError(f'{name} must be real')
    result = np.asarray(raw, dtype=float)
    if result.shape != (size,) or not np.isfinite(result).all():
        raise ValueError(f'{name} must be a finite vector of shape ({size},)')
    return result
