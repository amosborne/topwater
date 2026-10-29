"""Shared motion results and backend preparation contracts."""
from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class Motion:
    time_s: np.ndarray
    q: np.ndarray
    velocity_native: np.ndarray
    acceleration_native: np.ndarray
    position_m: np.ndarray
    orientation_wxyz: np.ndarray
    excitation_native: np.ndarray
    influence_state: np.ndarray | None = None
    force_components: dict = field(default_factory=dict)


@dataclass(frozen=True)
class PreparedMotion:
    hydrodynamics: object
    inputs: object = None


class ScratchCleanupError(RuntimeError):
    """A complete run cannot publish while scratch cleanup has failed."""


MODES = ('Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw')
