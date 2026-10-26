"""Six saved motion histories in the initial horizontal heading frame.

Only Result motion samples are consumed: no surface interpolation or solver
calculation belongs here. Absolute quaternions remain the orientation authority;
Euler curves are a display of departure from the initial attitude.
"""
import os
import math
from pathlib import Path
import tempfile
import warnings

import numpy as np
from scipy.spatial.transform import Rotation

from .workflow_data import Result


def _motion_ticks(values):
    """Seven symmetric nice decimal ticks, with at least 5% data headroom."""
    extent = float(np.max(np.abs(values)))
    if extent == 0.:
        step = 1.
    else:
        target = (extent/3)*1.05
        power = 10.**math.floor(math.log10(target))
        step = next(multiplier*power for multiplier in (1., 2., 2.5, 5., 10.)
                    if multiplier*power >= target)
    return np.arange(-3, 4, dtype=float)*step


def motion_histories(result: Result):
    """Return saved times and [surge, sway, heave, roll, pitch, yaw] curves.

    Translation is CG displacement from its initial location, expressed in a
    fixed inertial basis B: initial bow projected horizontally, port, world up.
    Rotation is B.T @ R(t) @ R(0).T @ B, displayed as unwrapped extrinsic XYZ
    Euler angles in degrees. Initial trim is removed. Euler angles are not
    angular velocities or uniquely defined at gimbal lock; that condition warns.
    """
    if not isinstance(result, Result):
        raise TypeError('motion plotting requires a Result')
    rotations = Rotation.from_quat(result.orientation_wxyz, scalar_first=True)
    initial = rotations[0]
    bow = initial.as_matrix()[:, 0]
    if np.linalg.norm(bow[:2]) <= 1e-12:
        raise ValueError('initial bow has no horizontal heading for the six-motion plot')
    heading = Rotation.from_euler('z', np.arctan2(bow[1], bow[0]))
    displacement = (result.position_m-result.position_m[0]) @ heading.as_matrix()
    relative = heading.inv() * rotations * initial.inv() * heading
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', UserWarning)
        angles = relative.as_euler('xyz')
    if caught:
        warnings.warn('Six-motion Euler labels are ambiguous at gimbal lock; '
                      'saved quaternions remain authoritative.', UserWarning, stacklevel=2)
    angles = np.rad2deg(np.unwrap(angles, axis=0))
    histories = np.column_stack((displacement, angles))
    if not np.isfinite(histories).all():
        raise ValueError('six-motion display produced nonfinite values')
    return result.time_s, histories


def save_motion_plot(result: Result, path):
    """Write a 2160 by 1170 PNG, atomically replacing a previous chart.

    All saved motion samples, including startup, are plotted. A private staging
    directory beside the destination directory keeps incomplete files outside
    canonical examples and permits same-volume atomic publication. A failed
    plot or replacement preserves the previous PNG. No global backend changes.
    """
    destination = Path(path)
    if destination.suffix.lower() != '.png':
        raise ValueError('motion plot output must name a .png file')
    times, histories = motion_histories(result)
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure(figsize=(12, 6.5), dpi=180, facecolor='white')
    FigureCanvasAgg(figure)
    axes = figure.subplots(2, 3, sharex=True, sharey='row')
    labels = ('Surge (m)', 'Sway (m)', 'Heave (m)',
              'Roll (deg)', 'Pitch (deg)', 'Yaw perturbation (deg)')
    for index, axis in enumerate(axes.flat):
        axis.plot(times, histories[:, index], color='#075c80', linewidth=1.1)
        axis.set_title(labels[index])
        axis.set_facecolor('white')
        axis.grid(True, alpha=.3)
        axis.set_xlim(times[0], times[-1])
        if index >= 3:
            axis.set_xlabel('Time from declared initial state (s)')
    for row, values in zip(axes, (histories[:, :3], histories[:, 3:])):
        ticks = _motion_ticks(values)
        for axis in row:
            axis.set_yticks(ticks)
            axis.set_ylim(ticks[0], ticks[-1])
    figure.suptitle('Complete six motion histories; initial attitude removed')
    figure.text(.5, .014,
                'CG displacement from initial position (m); fixed initial heading axes, z up. '
                'Unwrapped relative XYZ Euler angles (deg).',
                ha='center', fontsize=8)
    figure.tight_layout(rect=(0, .045, 1, .96))
    destination.parent.mkdir(parents=True, exist_ok=True)
    try:
        with tempfile.TemporaryDirectory(prefix='.topwater-motion-',
                                         dir=destination.absolute().parent.parent) as staging:
            temporary = Path(staging)/'result.png'
            figure.canvas.print_png(temporary)
            os.replace(temporary, destination)
    finally:
        figure.clear()
    return destination
