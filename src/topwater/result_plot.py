"""Six saved motion histories in the initial horizontal heading frame.

Motion uses saved Result samples; an optional original Sea supplies a static
spectrum beside them. No surface interpolation or solver calculation belongs
here. Absolute quaternions remain the orientation authority; Euler curves are a
display of departure from the initial attitude.
"""
import os
import math
from decimal import Decimal, ROUND_CEILING
from pathlib import Path
import tempfile
import warnings

import numpy as np
from scipy.spatial.transform import Rotation

from .workflow_data import Result, Sea


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


def _sea_spectrum(sea):
    """Original-carrier per-bin RMS elevations, wave height and default heading."""
    if not isinstance(sea, Sea):
        raise TypeError('sea-spectrum plotting requires the original Sea')
    frequency = sea.angular_frequencies_rad_s/(2*np.pi)
    # Keep world FROM bins fixed; rotate only the polar axes for presentation.
    angle = np.mod(sea.directions_rad+np.pi, 2*np.pi)
    angle_edges = np.linspace(0., 2*np.pi, 25)
    frequency_bins = max(20, int(np.ceil(float(frequency.max(initial=0.))/.025)))
    frequency_edges = np.arange(frequency_bins+1)*.025
    with np.errstate(over='ignore', invalid='ignore'):
        energy = sea.amplitudes_m**2
        total = float(np.sum(energy))
    if not np.isfinite(total):
        raise ValueError('sea carrier amplitudes exceed the finite spectrum display range')
    vector = np.array([energy @ np.cos(sea.directions_rad),
                       energy @ np.sin(sea.directions_rad)])
    # Match prepare_case's default bow heading without running case preparation.
    if total == 0.:
        default_heading = 0.
    elif np.linalg.norm(vector) <= 64*np.finfo(float).eps*total:
        raise ValueError('original sea has no unique default boat heading')
    else:
        default_heading = float((np.arctan2(vector[1], vector[0])+np.pi) % (2*np.pi))
    power, _, _ = np.histogram2d(frequency, angle,
        bins=(frequency_edges, angle_edges), weights=energy)
    return angle_edges, frequency_edges, np.sqrt(.5*power), float(4*np.sqrt(total/2)), default_heading


def _rms_color_maximum(rms_elevation):
    """Round outward to two significant decimal digits, never clip RMS data."""
    peak = float(np.max(rms_elevation, initial=0.))
    if peak == 0.:
        return 1.
    decimal_peak = Decimal(str(peak))
    quantum = Decimal(1).scaleb(decimal_peak.adjusted()-1)
    return float(decimal_peak.quantize(quantum, rounding=ROUND_CEILING))


def _add_sea_spectrum(figure, spectrum, actual_heading):
    """Add the right-hand polar panel and RMS-elevation colorbar to a wide figure."""
    from matplotlib.ticker import MaxNLocator
    angle, frequency, rms_elevation, _, default_heading = spectrum
    axes = figure.add_axes((12.20/17.5, 1.15/6.5, 5.1/17.5, 5.1/6.5),
                           projection='polar', label='sea spectrum')
    axes.set_facecolor('white')
    axes.grid(False)
    color_maximum = _rms_color_maximum(rms_elevation)
    mesh = axes.pcolormesh(angle, frequency, rms_elevation, cmap='cividis', vmin=0.,
        vmax=color_maximum, shading='flat', rasterized=True)
    axes.set_theta_offset(np.pi/2-default_heading)
    axes.set_theta_direction(1)
    axes.set_xticks([])
    axes.xaxis.grid(False)
    axes.set_ylim(0., frequency[-1])
    rings = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]).tick_values(0., frequency[-1])
    rings = rings[(rings > 0.) & (rings < frequency[-1]-1e-12)]
    axes.set_yticks(rings, [f'{value:g}'+(' Hz' if index == 0 else '')
                           for index, value in enumerate(rings)])
    axes.set_rlabel_position(np.rad2deg(-np.pi/2-axes.get_theta_offset()))
    axes.tick_params(axis='y', labelcolor='black', labelsize=10)
    for label in axes.get_yticklabels():
        label.set_fontfamily('sans-serif')
    axes.grid(axis='y', color='#d0d0d0', alpha=.8, linewidth=.65)
    axes.plot([actual_heading, actual_heading], [0., frequency[-1]], color='black',
              linestyle=':', linewidth=1.5, label='initial boat heading', zorder=3)
    colorbar_axes = figure.add_axes((12.20/17.5, .65/6.5, 5.1/17.5, .16/6.5))
    bar = figure.colorbar(mesh, cax=colorbar_axes, orientation='horizontal',
                         label='Bin RMS elevation (m)')
    ticks = MaxNLocator(nbins=6, steps=[1, 2, 2.5, 5, 10]).tick_values(0., color_maximum)
    interior = ticks[(ticks > color_maximum*1e-12) &
                     (ticks < color_maximum*(1-1e-12))]
    bar.set_ticks(np.r_[0., interior, color_maximum])


def save_motion_plot(result: Result, path, *, sea=None):
    """Write a six-motion PNG, optionally with the original sea beside it.

    Without Sea the chart remains 2160 by 1170 pixels. Supplying Sea adds a
    spectrum column, expanding width to 3150 pixels without shrinking the six
    motion panels. No carrier data is read from or written into Result.

    All saved motion samples, including startup, are plotted. A private staging
    directory beside the destination directory keeps incomplete files outside
    canonical examples and permits same-volume atomic publication. A failed
    plot or replacement preserves the previous PNG. No global backend changes.
    """
    destination = Path(path)
    if destination.suffix.lower() != '.png':
        raise ValueError('motion plot output must name a .png file')
    times, histories = motion_histories(result)
    spectrum = None if sea is None else _sea_spectrum(sea)
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    width = 12. if spectrum is None else 17.5
    motion_fraction = 12/width
    figure = Figure(figsize=(width, 6.5), dpi=180, facecolor='white')
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
            axis.set_xlabel('Time (s)')
    for row, values in zip(axes, (histories[:, :3], histories[:, 3:])):
        ticks = _motion_ticks(values)
        for axis in row:
            axis.set_yticks(ticks)
            axis.set_ylim(ticks[0], ticks[-1])
    figure.tight_layout(rect=(0, 0, motion_fraction, 1))
    if spectrum is not None:
        bow = Rotation.from_quat(result.orientation_wxyz[0], scalar_first=True).as_matrix()[:, 0]
        _add_sea_spectrum(figure, spectrum, float(np.arctan2(bow[1], bow[0])))
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
