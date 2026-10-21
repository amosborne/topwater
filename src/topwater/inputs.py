"""Validated physical inputs and their JSON representation.

Arrays are copied and marked read-only. Geometry is the complete, initially dry
and wet hull, expressed in the body frame. These inputs are not a solver state
or a restart file. See class docstrings for units and frame ownership.
"""

from dataclasses import dataclass, field
import json
from numbers import Integral
from types import MappingProxyType
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike
from scipy.sparse import coo_array
from scipy.sparse.csgraph import connected_components

from ._validation import array, quaternion


@dataclass(frozen=True, eq=False)
class HullMesh:
    """One closed, connected triangular shell with outward body normals.

    vertices: (n, 3) body coordinates [m]. faces: (m, 3) integer vertex indices,
    ordered so their right-hand normals point from the body into the water/air.
    Rejects open/non-manifold edges, inconsistent winding, duplicate or unused
    vertices, degenerate panels and non-positive enclosed volume. Does not test
    self-intersections or suitability for free-surface intersection/remeshing.
    """

    vertices: ArrayLike
    faces: ArrayLike

    def __post_init__(self):
        vertices = array(self.vertices, (None, 3), "vertices")
        raw_faces = np.asarray(self.faces)
        if raw_faces.ndim != 2 or raw_faces.shape[1] != 3:
            raise ValueError("faces must have shape (m, 3)")
        if raw_faces.dtype.kind not in "iu":
            raise ValueError("faces must contain integer vertex indices")
        if len(vertices) < 4 or len(raw_faces) < 4:
            raise ValueError("mesh must enclose a volume with at least four vertices/faces")
        if np.any(raw_faces < 0) or np.any(raw_faces >= len(vertices)):
            raise ValueError("faces contain out-of-range vertex indices")
        faces = np.array(raw_faces, dtype=np.int64, copy=True, order="C")
        if len(np.unique(vertices, axis=0)) != len(vertices):
            raise ValueError("vertices must be unique; weld coincident vertices first")
        if len(np.unique(faces)) != len(vertices):
            raise ValueError("mesh contains unused vertices")
        if len(np.unique(np.sort(faces, axis=1), axis=0)) != len(faces):
            raise ValueError("mesh contains duplicate faces")

        triangles = vertices[faces]
        a = triangles[:, 1] - triangles[:, 0]
        b = triangles[:, 2] - triangles[:, 0]
        area = np.linalg.norm(np.cross(a, b), axis=1)
        tolerance = 64 * np.finfo(float).eps
        if np.any(area <= tolerance * np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)):
            raise ValueError("mesh contains degenerate triangles")

        # Pair shared edges with compiled array operations, not per-panel Python.
        directed = np.concatenate([
            faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]],
        ])
        edges = np.sort(directed, axis=1)
        order = np.lexsort((edges[:, 1], edges[:, 0]))
        ordered = edges[order]
        starts = np.r_[0, 1 + np.flatnonzero(np.any(ordered[1:] != ordered[:-1], axis=1))]
        if np.any(np.diff(np.r_[starts, len(edges)]) != 2):
            raise ValueError("mesh must be closed with exactly two faces per edge")
        forward = (directed[:, 0] < directed[:, 1])[order].reshape(-1, 2)
        if np.any(forward[:, 0] == forward[:, 1]):
            raise ValueError("mesh faces must have consistent winding")
        owners = np.tile(np.arange(len(faces)), 3)[order].reshape(-1, 2)
        graph = coo_array(
            (np.ones(len(owners), dtype=bool), (owners[:, 0], owners[:, 1])),
            shape=(len(faces), len(faces)),
        ).tocsr()
        if connected_components(graph, directed=False, return_labels=False) != 1:
            raise ValueError("mesh must contain one connected shell")

        # Translate and scale before summing tetrahedra to reduce cancellation.
        relative = vertices - vertices[0]
        scale = np.max(np.abs(relative))
        normalized = relative / scale
        p, q, r = normalized[faces].transpose(1, 0, 2)
        volume = np.sum(np.einsum("ij,ij->i", p, np.cross(q, r))) / 6
        if not np.isfinite(volume) or volume <= tolerance:
            raise ValueError("mesh must have outward winding and positive enclosed volume")
        faces.setflags(write=False)
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "faces", faces)


@dataclass(frozen=True, eq=False)
class RigidBody:
    """Mass [kg], CG [m] in mesh body coordinates, inertia about CG [kg m^2].

    inertia_cg is a full symmetric tensor in the body axes, not moments about
    the mesh origin. It must be positive definite with physically realizable
    principal moments (each no larger than the sum of the other two).
    """

    mesh: HullMesh
    mass: float
    cg: ArrayLike
    inertia_cg: ArrayLike

    def __post_init__(self):
        if not isinstance(self.mesh, HullMesh):
            raise ValueError("mesh must be a HullMesh")
        mass = array(self.mass, (), "mass").item()
        if mass <= 0:
            raise ValueError("mass must be positive")
        cg = array(self.cg, (3,), "cg")
        inertia = array(self.inertia_cg, (3, 3), "inertia_cg")
        scale = np.max(np.abs(inertia))
        if scale == 0 or not np.allclose(inertia, inertia.T, rtol=0, atol=1e-12 * scale):
            raise ValueError("inertia_cg must be symmetric and positive definite")
        inertia = (inertia + inertia.T) / 2
        moments = np.linalg.eigvalsh(inertia / scale)
        if moments[0] <= 0:
            raise ValueError("inertia_cg must be positive definite")
        if moments[2] > moments[0] + moments[1] + 1e-12:
            raise ValueError("inertia_cg principal moments must satisfy the triangle inequality")
        inertia.setflags(write=False)
        object.__setattr__(self, "mass", mass)
        object.__setattr__(self, "cg", cg)
        object.__setattr__(self, "inertia_cg", inertia)


@dataclass(frozen=True, eq=False)
class WaveComponents:
    """Explicit incident-wave realization, with four equal-length 1-D arrays.

    amplitudes [m] are nonnegative; angular_frequencies [rad/s] are positive;
    directions [rad] are propagation angles counterclockwise from world east;
    phases [rad] use cos(k*(x*cos(direction)+y*sin(direction))-omega*t+phase).
    Empty arrays mean calm water. Directions/phases are retained as supplied.
    A seed records provenance; it does not regenerate or override explicit phases.
    provenance contains string keys/values such as source, station and timestamp.
    """

    amplitudes: ArrayLike
    angular_frequencies: ArrayLike
    directions: ArrayLike
    phases: ArrayLike
    seed: int | None = None
    provenance: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self):
        amplitudes = array(self.amplitudes, (None,), "amplitudes")
        object.__setattr__(self, "amplitudes", amplitudes)
        for name in ("angular_frequencies", "directions", "phases"):
            object.__setattr__(self, name, array(getattr(self, name), amplitudes.shape, name))
        if np.any(amplitudes < 0):
            raise ValueError("amplitudes must be nonnegative")
        if np.any(self.angular_frequencies <= 0):
            raise ValueError("angular_frequencies must be positive")
        if self.seed is not None:
            if isinstance(self.seed, bool) or not isinstance(self.seed, Integral) or self.seed < 0:
                raise ValueError("seed must be a nonnegative integer or None")
            object.__setattr__(self, "seed", int(self.seed))
        if not isinstance(self.provenance, Mapping) or not all(
            isinstance(key, str) and isinstance(value, str)
            for key, value in self.provenance.items()
        ):
            raise ValueError("provenance must map strings to strings")
        object.__setattr__(self, "provenance", MappingProxyType(dict(self.provenance)))


@dataclass(frozen=True, eq=False)
class InitialState:
    """CG position [m] and velocity [m/s] in world axes; angular velocity [rad/s]
    in body axes. orientation is unit (w, x, y, z), actively rotating body to
    world. Defaults are the identity pose at rest. Mesh points transform as
    x_world = position + R @ (x_body - body.cg).
    """

    position: ArrayLike = (0.0, 0.0, 0.0)
    orientation: ArrayLike = (1.0, 0.0, 0.0, 0.0)
    linear_velocity: ArrayLike = (0.0, 0.0, 0.0)
    angular_velocity: ArrayLike = (0.0, 0.0, 0.0)

    def __post_init__(self):
        for name in ("position", "linear_velocity", "angular_velocity"):
            object.__setattr__(self, name, array(getattr(self, name), (3,), name))
        object.__setattr__(self, "orientation", quaternion(self.orientation))


@dataclass(frozen=True, eq=False)
class SimulationInputs:
    """Physical inputs only; serialization version 1 is not a restart format."""

    body: RigidBody
    waves: WaveComponents
    initial_state: InitialState = field(default_factory=InitialState)

    def __post_init__(self):
        for name, kind in (("body", RigidBody), ("waves", WaveComponents),
                           ("initial_state", InitialState)):
            if not isinstance(getattr(self, name), kind):
                raise ValueError(f"{name} must be a {kind.__name__}")

    def to_json(self) -> str:
        """Return portable JSON with SI arrays, actual phases and provenance."""
        return json.dumps({
            "format_version": 1,
            "body": {
                "mesh": {"vertices": self.body.mesh.vertices.tolist(),
                         "faces": self.body.mesh.faces.tolist()},
                "mass": self.body.mass,
                "cg": self.body.cg.tolist(),
                "inertia_cg": self.body.inertia_cg.tolist(),
            },
            "waves": {
                "amplitudes": self.waves.amplitudes.tolist(),
                "angular_frequencies": self.waves.angular_frequencies.tolist(),
                "directions": self.waves.directions.tolist(),
                "phases": self.waves.phases.tolist(),
                "seed": self.waves.seed,
                "provenance": dict(self.waves.provenance),
            },
            "initial_state": {
                name: getattr(self.initial_state, name).tolist()
                for name in ("position", "orientation", "linear_velocity", "angular_velocity")
            },
        }, allow_nan=False, indent=2)

    @classmethod
    def from_json(cls, text: str) -> "SimulationInputs":
        """Read version 1 JSON, applying the same validation as constructors."""
        try:
            data = json.loads(text)
            if not isinstance(data, dict) or set(data) != {
                "format_version", "body", "waves", "initial_state"
            }:
                raise ValueError("expected format_version, body, waves and initial_state")
            if type(data["format_version"]) is not int or data["format_version"] != 1:
                raise ValueError("unsupported input format_version; expected 1")
            body = dict(data["body"])
            body["mesh"] = HullMesh(**body["mesh"])
            return cls(
                body=RigidBody(**body),
                waves=WaveComponents(**data["waves"]),
                initial_state=InitialState(**data["initial_state"]),
            )
        except (KeyError, TypeError, ValueError, OverflowError) as error:
            raise ValueError(f"Invalid simulation inputs: {error}") from error
