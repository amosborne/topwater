"""Chosen numerical and output controls, with SI units and one default owner."""
from dataclasses import dataclass
import math

import numpy as np

DEFAULT_DWELL_S = 100.
VIDEO_SIZE_PX = (800, 800)
SURFACE_PIXELS_PER_SAMPLE = 10

@dataclass(frozen=True)
class OfflineRadiationSettings:
    """Full-complex OFFLINE FFT policy, requiring case-specific qualification.

    Defaults use exactly the saved native samples, without padding or resampling.
    An explicit period_s selects the historical padded policy; its optional
    future polynomial is numerical continuation, not predicted body motion.
    """
    period_s: float | None = None
    fft_dt_s: float | None = None
    future_s: float = 0.
    chunk_points: int = 128


@dataclass(frozen=True)
class SolveOptions:
    time_step_s: float = .01
    render_every: int = 10
    packet_tolerance_m: float = .01
    surface_padding_beams: float = 1.
    coefficient_step_rad_s: float = .1
    radiation_bandwidth_rad_s: float = 12.
    batch_points: int | None = None
    cpu_workers: int | None = None

    def __post_init__(self):
        for name in ('time_step_s', 'packet_tolerance_m',
                     'coefficient_step_rad_s',
                     'radiation_bandwidth_rad_s'):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f'{name} must be positive and finite')
        if not math.isfinite(self.surface_padding_beams) or self.surface_padding_beams < 0:
            raise ValueError('surface_padding_beams must be nonnegative and finite')
        for name in ('render_every', 'batch_points', 'cpu_workers'):
            value = getattr(self, name)
            if (value is not None or name == 'render_every') and (
                    not isinstance(value, int) or isinstance(value, bool) or value < 1):
                raise ValueError(f'{name} must be a positive integer')

    @property
    def surface_interval_s(self):
        return self.render_every*self.time_step_s

    def radiation_frequencies(self, packet_maximum):
        step = self.coefficient_step_rad_s
        last = max(2, math.ceil(max(self.radiation_bandwidth_rad_s, packet_maximum)/step))
        return np.arange(1, last+1)*step

    def radiation_settings(self, chunk_points):
        return OfflineRadiationSettings(
            chunk_points=self.batch_points if self.batch_points is not None else chunk_points)


DEFAULTS = SolveOptions()
