"""Small numerical NPZ contracts for boats, realized seas and final results.

No pickle, solver markers, generic metadata, or implementation objects are
stored. Units/frames are fixed by this module. Saves publish a complete archive
atomically; validation happens before any destination is replaced.
"""

from dataclasses import dataclass
import os
from pathlib import Path
import tempfile
from zipfile import BadZipFile, ZipFile
import zlib

import numpy as np

from ._validation import array
from .boat_geometry import physical_hull
from .inputs import RigidBody, WaveComponents


def _integer(value, shape, name):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iu' or raw.ndim != len(shape) or any(
            want is not None and got != want for got, want in zip(raw.shape, shape)):
        raise ValueError(f'{name} must be an integer array with shape {shape}')
    if raw.dtype.kind == 'u' and np.any(raw > np.iinfo(np.int64).max):
        raise ValueError(f'{name} exceeds signed 64-bit index range')
    result = np.array(raw, dtype=np.int64, order='C', copy=True)
    result.setflags(write=False)
    return result


def _times(value, name):
    times = array(value, (None,), name)
    if len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0):
        raise ValueError(f'{name} must start at zero and contain at least two strictly increasing times')
    return times


def _read_archive(cls, path):
    path = Path(path)
    with ZipFile(path) as archive:
        members = archive.namelist()
        expected = {name+'.npy' for name in cls.members}
        if len(members) != len(expected) or set(members) != expected:
            raise ValueError(f'{cls.__name__} archive must contain exactly {", ".join(cls.members)}')
    with np.load(path, allow_pickle=False) as data:
        return cls(**{name: data[name] for name in cls.members})


def _load(cls, path):
    try:
        return _read_archive(cls, path)
    except (BadZipFile, EOFError, zlib.error) as error:
        raise ValueError(f'invalid or truncated {cls.__name__} NPZ archive: {path}') from error


def _save(value, path, overwrite):
    path = Path(path)
    if path.suffix.lower() != '.npz':
        raise ValueError('output must name a complete .npz file')
    if not isinstance(overwrite, bool):
        raise TypeError('overwrite must be a bool')
    if path.exists() and not overwrite:
        raise FileExistsError(f'output already exists: {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    # Prefer the OS temp directory so example folders contain only final files.
    # Cross-filesystem publication needs a same-volume final staging file.
    staging = None if os.stat(path.parent).st_dev == os.stat(tempfile.gettempdir()).st_dev else path.parent
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='wb', dir=staging,
                                         prefix='.'+path.stem+'-', suffix='.npz', delete=False) as stream:
            temporary = Path(stream.name)
            np.savez(stream, **{name: getattr(value, name) for name in value.members})
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary, path)
        elif os.name == 'nt':
            # Windows rename is atomic and refuses an existing destination.
            # Hard-link creation may be unavailable in a restricted workspace.
            os.rename(temporary, path)
        else:
            # A hard link publishes without the check-then-rename overwrite race.
            os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


class _Archive:
    @classmethod
    def load(cls, path):
        return _load(cls, path)

    def save(self, path, *, overwrite=False):
        _save(self, path, overwrite)


@dataclass(frozen=True, eq=False)
class Boat(_Archive):
    """Uncapped physical hull in body bow/port/up axes, all quantities SI.

    cg_m is expressed in the same coordinates as vertices_m. The full inertia
    tensor is about this CG in body axes. One opening is inferred;
    users do not supply artificial caps or a separate opening declaration.
    """

    vertices_m: np.ndarray
    faces: np.ndarray
    mass_kg: float
    cg_m: np.ndarray
    inertia_cg_kg_m2: np.ndarray

    members = ('vertices_m', 'faces', 'mass_kg', 'cg_m', 'inertia_cg_kg_m2')

    def __post_init__(self):
        mesh = physical_hull(self.vertices_m, self.faces)
        body = RigidBody(mesh, self.mass_kg, self.cg_m, self.inertia_cg_kg_m2)
        # Preserve exactly the physical input, not the internally generated cap.
        object.__setattr__(self, 'vertices_m', array(self.vertices_m, (None, 3), 'vertices_m'))
        object.__setattr__(self, 'faces', _integer(self.faces, (None, 3), 'faces'))
        object.__setattr__(self, 'mass_kg', body.mass)
        object.__setattr__(self, 'cg_m', body.cg)
        object.__setattr__(self, 'inertia_cg_kg_m2', body.inertia_cg)
        object.__setattr__(self, '_body', body)

    def to_body(self):
        """Return immutable legacy physics inputs with a nonphysical internal cap."""
        return self._body


@dataclass(frozen=True, eq=False)
class Sea(_Archive):
    """Original carrier inputs, independent of a boat and packet timing.

    Directions propagate CCW from world +x (east). The physical target phase is
    k dot x - omega*(t-dwell_midpoint) + phases_rad. Thus the original phase
    values apply DIRECTLY at the selected dwell midpoint and wave reference
    point; zero specifies a target carrier crest. Fourier sign conventions may
    be converted internally, without shifting that reference time or changing
    these arrays. Finite truncation may depart from the target. Empty arrays
    mean calm water. Solve constructs packets from these four original vectors.
    """

    amplitudes_m: np.ndarray
    angular_frequencies_rad_s: np.ndarray
    directions_rad: np.ndarray
    phases_rad: np.ndarray

    members = ('amplitudes_m', 'angular_frequencies_rad_s', 'directions_rad', 'phases_rad')

    def __post_init__(self):
        waves = WaveComponents(self.amplitudes_m, self.angular_frequencies_rad_s,
                               self.directions_rad, self.phases_rad)
        for name, field in zip(self.members, ('amplitudes', 'angular_frequencies', 'directions', 'phases')):
            object.__setattr__(self, name, getattr(waves, field))

    def to_waves(self):
        return WaveComponents(self.amplitudes_m, self.angular_frequencies_rad_s,
                              self.directions_rad, self.phases_rad)


@dataclass(frozen=True, eq=False)
class Result(_Archive):
    """Self-contained rigid body pose and always-present world water geometry.

    time_s samples absolute world CG position_m and active body-to-world
    orientation_wxyz. Hull vertices are body-frame physical panels (the solved
    mesh may differ from the original boat). World hull points are
    position + (hull_vertices_m-cg_m) @ rotation.T.

    Surface samples have their own time axis with identical start/end times.
    Each frame is a nonempty triangular mesh packed into surface_vertices_m and
    surface_faces. Its faces use LOCAL indices, not offsets into the packed
    vertices. Both offset vectors start at zero and have S+1 entries. Variable
    counts and connectivity are supported; no time interpolation is implicit.
    Export need not cut the surface at the hull: portions inside the posed hull
    are invalid fluid. Waterline intersections are computed separately. A native
    producer may supply its own evolving mesh without an export remeshing step.
    The archive is neither a restart nor a physics-accuracy certificate.
    """

    time_s: np.ndarray
    position_m: np.ndarray
    orientation_wxyz: np.ndarray
    hull_vertices_m: np.ndarray
    hull_faces: np.ndarray
    cg_m: np.ndarray
    surface_time_s: np.ndarray
    surface_vertices_m: np.ndarray
    surface_faces: np.ndarray
    surface_vertex_offsets: np.ndarray
    surface_face_offsets: np.ndarray

    members = ('time_s', 'position_m', 'orientation_wxyz', 'hull_vertices_m', 'hull_faces',
               'cg_m', 'surface_time_s', 'surface_vertices_m', 'surface_faces',
               'surface_vertex_offsets', 'surface_face_offsets')

    def __post_init__(self):
        times = _times(self.time_s, 'time_s')
        surface_times = _times(self.surface_time_s, 'surface_time_s')
        if times[-1] != surface_times[-1]:
            raise ValueError('body and surface time axes must have identical endpoints')
        position = array(self.position_m, (len(times), 3), 'position_m')
        orientation = array(self.orientation_wxyz, (len(times), 4), 'orientation_wxyz')
        if not np.allclose(np.linalg.norm(orientation, axis=1), 1., rtol=0., atol=1e-10):
            raise ValueError('orientation_wxyz must contain unit body-to-world quaternions')
        physical_hull(self.hull_vertices_m, self.hull_faces)
        hull = array(self.hull_vertices_m, (None, 3), 'hull_vertices_m')
        hull_faces = _integer(self.hull_faces, (None, 3), 'hull_faces')
        cg = array(self.cg_m, (3,), 'cg_m')
        vertices = array(self.surface_vertices_m, (None, 3), 'surface_vertices_m')
        faces = _integer(self.surface_faces, (None, 3), 'surface_faces')
        vertex_offsets = _integer(self.surface_vertex_offsets, (len(surface_times)+1,), 'surface_vertex_offsets')
        face_offsets = _integer(self.surface_face_offsets, (len(surface_times)+1,), 'surface_face_offsets')
        for offsets, count, name in ((vertex_offsets, len(vertices), 'surface_vertex_offsets'),
                                      (face_offsets, len(faces), 'surface_face_offsets')):
            if offsets[0] != 0 or offsets[-1] != count or np.any(np.diff(offsets) <= 0):
                raise ValueError(f'{name} must start at zero, increase strictly and span its packed array')
        for i in range(len(surface_times)):
            frame = vertices[vertex_offsets[i]:vertex_offsets[i+1]]
            connectivity = faces[face_offsets[i]:face_offsets[i+1]]
            if len(frame) < 3 or np.any(connectivity < 0) or np.any(connectivity >= len(frame)):
                raise ValueError(f'surface frame {i} has invalid local triangle indices')
            if len(np.unique(np.sort(connectivity, axis=1), axis=0)) != len(connectivity):
                raise ValueError(f'surface frame {i} has duplicate faces')
            if len(np.unique(connectivity)) != len(frame):
                raise ValueError(f'surface frame {i} contains unused vertices')
            triangles = frame[connectivity]
            a, b = triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0]
            magnitude = np.linalg.norm(np.cross(a, b), axis=1)
            if np.any(magnitude <= 64*np.finfo(float).eps*np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1)):
                raise ValueError(f'surface frame {i} has degenerate triangles')
        values = (times, position, orientation, hull, hull_faces, cg, surface_times,
                  vertices, faces, vertex_offsets, face_offsets)
        for name, value in zip(self.members, values):
            object.__setattr__(self, name, value)

    def surface_frame(self, index):
        """Return read-only views of a saved frame; do not interpolate."""
        if not isinstance(index, (int, np.integer)) or isinstance(index, (bool, np.bool_)):
            raise TypeError('surface frame index must be an integer')
        if index < 0 or index >= len(self.surface_time_s):
            raise IndexError('surface frame index out of range')
        return (self.surface_vertices_m[self.surface_vertex_offsets[index]:self.surface_vertex_offsets[index+1]],
                self.surface_faces[self.surface_face_offsets[index]:self.surface_face_offsets[index+1]])
