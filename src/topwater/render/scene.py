"""Producer-neutral mesh animation using shared native CPU/GPU scene drawing.

Inputs are already prepared world geometry, never simulation parameters or loads.
One sequence has fixed triangular connectivity; another hull is another input.
"""
from dataclasses import dataclass, field
import json

import numpy as np
from scipy.spatial import ConvexHull, QhullError, cKDTree


@dataclass
class SceneSequence:
    """World meshes in metres, right handed with z up, at increasing times.

    ``*_points`` have shape (time, vertex, 3); ``*_faces`` are (triangle, 3)
    integer indices. Hull faces may leave a gunwale open, but must lie on the
    convex envelope. Rendering projects these faces and their actual waterline.
    Water faces should be consistently oriented toward the air.
    """
    time_s: np.ndarray
    hull_points: np.ndarray
    hull_faces: np.ndarray
    water_points: np.ndarray
    water_faces: np.ndarray
    metadata: dict = field(default_factory=dict)
    telemetry: dict = field(default_factory=dict)

    def __post_init__(self):
        self.time_s = np.asarray(self.time_s, dtype=float)
        if self.time_s.ndim != 1 or not len(self.time_s) or not np.isfinite(self.time_s).all():
            raise ValueError('time_s must be a nonempty finite vector')
        if np.any(np.diff(self.time_s) <= 0):
            raise ValueError('time_s must increase strictly')
        for name in ('hull', 'water'):
            points = np.asarray(getattr(self, name + '_points'), dtype=float)
            raw_faces = np.asarray(getattr(self, name + '_faces'))
            if points.ndim != 3 or points.shape[0] != len(self.time_s) or points.shape[2] != 3 or not np.isfinite(points).all():
                raise ValueError(f'{name}_points must be finite (time, vertex, 3)')
            if raw_faces.ndim != 2 or raw_faces.shape[1] != 3 or not len(raw_faces) or not np.issubdtype(raw_faces.dtype, np.integer):
                raise ValueError(f'{name}_faces must be integer (triangle, 3)')
            if raw_faces.min() < 0 or raw_faces.max() >= points.shape[1]:
                raise ValueError(f'{name}_faces index outside vertex array')
            setattr(self, name + '_points', points)
            setattr(self, name + '_faces', raw_faces.astype(np.int64, copy=False))
        if self.metadata.get('units', 'm') != 'm' or self.metadata.get('coordinate_system', 'right-handed-z-up') != 'right-handed-z-up':
            raise ValueError('Normalize geometry to metres and right-handed-z-up before rendering')
        for name, values in self.telemetry.items():
            values = np.asarray(values, dtype=float)
            if values.shape != self.time_s.shape or not np.isfinite(values).all():
                raise ValueError(f'telemetry {name} must be finite with one value per time')
            self.telemetry[name] = values

    def save(self, path):
        """Save prepared geometry once; style changes can reuse it without physics."""
        np.savez(path, time_s=self.time_s, hull_points=self.hull_points,
                 hull_faces=self.hull_faces, water_points=self.water_points,
                 water_faces=self.water_faces, metadata=json.dumps(self.metadata),
                 **{'telemetry__' + name: values for name, values in self.telemetry.items()})

    @classmethod
    def load(cls, path):
        with np.load(path, allow_pickle=False) as data:
            return cls(**{key: data[key] for key in ('time_s', 'hull_points', 'hull_faces', 'water_points', 'water_faces')},
                       metadata=json.loads(str(data['metadata'])),
                       telemetry={key.removeprefix('telemetry__'): data[key] for key in data.files if key.startswith('telemetry__')})


def convex_envelope(points, faces):
    """Validate convex boundary faces and return outward oriented closure faces.

    An internal face is rejected rather than silently replacing a nonconvex boat
    by its convex hull. Open boundaries are permitted for the visual hull.
    """
    points = np.asarray(points, dtype=float)
    faces = np.asarray(faces, dtype=int)
    scale = np.linalg.norm(np.ptp(points, axis=0))
    try:
        envelope = ConvexHull(points)
    except QhullError as error:
        raise ValueError('Hull must enclose a three-dimensional convex volume') from error
    if np.any(np.linalg.norm(np.cross(points[faces[:, 1]] - points[faces[:, 0]], points[faces[:, 2]] - points[faces[:, 0]]), axis=1) <= 1e-12 * scale**2):
        raise ValueError('Hull contains degenerate triangles')
    oriented = _outward(points, faces)
    triangle = points[oriented]
    normal = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    # A boundary triangle's outward plane must match a convex-envelope plane.
    # Query by normal, then test its vertices against that support plane. This
    # avoids an O(faces * envelope_planes) dense allocation on refined meshes.
    _, nearest = cKDTree(envelope.equations[:, :3]).query(normal)
    planes = envelope.equations[nearest]
    residual = abs(np.einsum('fvi,fi->fv', triangle, planes[:, :3]) + planes[:, 3, None])
    if np.any(residual.max(axis=1) > 1e-8 * max(scale, 1e-12)):
        raise ValueError('Hull faces must lie on a convex boundary')
    closure = envelope.simplices.copy()
    return _outward(points, closure), oriented


def _outward(points, faces):
    result = faces.copy()
    triangle = points[result]
    normal = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    flip = np.einsum('ij,ij->i', normal, triangle.mean(axis=1) - points.mean(axis=0)) < 0
    result[flip] = result[flip][:, [0, 2, 1]]
    return result


@dataclass
class RenderSettings:
    """Fixed top-down presentation settings, independent of the producer."""
    framing: float = 1.45


def render_sequence(scene, output_path, settings=None):
    """Render the prepared sequence through shared CPU/GPU scene drawing."""
    from pathlib import Path
    from .frames import VideoFrames
    from .workflow import render_result
    settings = RenderSettings() if settings is None else settings
    output=Path(output_path)
    if output.suffix.lower()!='.mp4':
        raise ValueError('Output must name the complete MP4 file, including .mp4')
    if output.exists():
        raise FileExistsError(f'Output already exists: {output}')
    return render_result(VideoFrames.from_scene(scene,settings),output,settings)
