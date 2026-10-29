import json

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from topwater import HullMesh
from topwater import InitialState
from topwater import RigidBody
from topwater import SimulationInputs
from topwater import WaveComponents
from topwater.geometry.transforms import body_to_world


def test_opening_metadata_is_owned_and_matches_closure(input_data):
    source = input_data['body']['mesh']
    mesh = HullMesh(**source)
    original = mesh.opening_rim.copy()
    source['opening_rim'][0] = source['opening_rim'][1]
    assert_array_equal(mesh.opening_rim, original)
    assert not mesh.opening_rim.flags.writeable
    assert not mesh.artificial_faces.flags.writeable
    assert not mesh.physical_faces.flags.writeable
    assert len(mesh.physical_faces)+len(mesh.artificial_faces) == len(mesh.faces)
    with pytest.raises(ValueError, match='duplicate'):
        HullMesh(**source)


@pytest.mark.parametrize('change', ['missing_rim', 'missing_cap', 'wrong_loop', 'wrong_cap'])
def test_invalid_opening_identity(input_data, change):
    source = input_data['body']['mesh']
    if change == 'missing_rim':
        source['opening_rim'] = []
    elif change == 'missing_cap':
        source['artificial_faces'] = []
    elif change == 'wrong_loop':
        source['opening_rim'][0], source['opening_rim'][1] = source['opening_rim'][1], source['opening_rim'][0]
    else:
        source['artificial_faces'] = source['artificial_faces'][:-1]
    with pytest.raises(ValueError, match='opening|boundary'):
        HullMesh(**source)


def test_json_requires_explicit_opening_identity(input_data):
    del input_data['body']['mesh']['artificial_faces']
    with pytest.raises(ValueError, match='explicitly identify'):
        SimulationInputs.from_json(json.dumps(input_data))


def test_mesh_accepts_rigid_transforms_and_unit_changes(mesh):
    q = [np.sqrt(0.5), 0, 0, np.sqrt(0.5)]
    for scale in (1e-6, 1, 1e6):
        vertices = body_to_world(mesh.vertices, q) * scale + np.array([4, -5, 6]) * scale
        moved = HullMesh(vertices, mesh.faces)
        assert_array_equal(moved.faces, mesh.faces)
    assert moved.vertices.flags.c_contiguous
    assert moved.faces.flags.c_contiguous


@pytest.mark.parametrize(("change", "message"), [
    ("open", "closed"), ("flip_one", "consistent winding"),
    ("flip_all", "outward winding"), ("out_of_range", "out-of-range"),
    ("negative_index", "out-of-range"), ("float_indices", "integer"),
    ("boolean_indices", "integer"), ("duplicate_face", "duplicate faces"),
    ("degenerate", "degenerate"), ("unused", "unused"),
    ("duplicate_vertex", "unique"), ("nonfinite", "finite"),
    ("disconnected", "connected shell"),
    ("nonmanifold", "exactly two faces"),
])
def test_invalid_meshes(mesh, change, message):
    vertices = mesh.vertices.copy()
    faces = mesh.faces.copy()
    if change == "open":
        faces = faces[:-1]
    elif change == "flip_one":
        faces[0] = faces[0, ::-1]
    elif change == "flip_all":
        faces = faces[:, ::-1]
    elif change == "out_of_range":
        faces[0, 0] = len(vertices)
    elif change == "negative_index":
        faces[0, 0] = -1
    elif change == "float_indices":
        faces = faces.astype(float)
    elif change == "boolean_indices":
        faces = faces.astype(bool)
    elif change == "duplicate_face":
        faces = np.vstack([faces, faces[0]])
    elif change == "degenerate":
        a, b, c = faces[0]
        vertices[c] = (vertices[a] + vertices[b]) / 2
    elif change == "unused":
        vertices = np.vstack([vertices, [2, 2, 2]])
    elif change == "duplicate_vertex":
        vertices[3] = vertices[2]
    elif change == "nonfinite":
        vertices[0, 0] = np.nan
    elif change == "disconnected":
        vertices = np.vstack([vertices, vertices + [10, 0, 0]])
        faces = np.vstack([faces, faces + len(mesh.vertices)])
    elif change == "nonmanifold":
        a, b = faces[0, :2]
        corner = next(i for i in range(len(vertices))
            if i not in (a, b)
            and np.linalg.norm(np.cross(vertices[b]-vertices[a], vertices[i]-vertices[a])) > 0
            and not np.any(np.all(np.sort(faces, axis=1) == np.sort([a, b, i]), axis=1)))
        faces = np.vstack([faces, [a, b, corner]])
    with pytest.raises(ValueError, match=message):
        HullMesh(vertices, faces)


@pytest.mark.parametrize("mass", [0, -1, np.inf, np.nan, True, "6", [6]])
def test_invalid_mass(mesh, mass):
    with pytest.raises(ValueError, match="mass"):
        RigidBody(mesh, mass, [0, 0, 0], np.eye(3))


@pytest.mark.parametrize(("inertia", "message"), [
    (np.diag([1, 1, 0]), "positive definite"),
    (np.diag([1, 1, -1]), "positive definite"),
    (np.diag([1, 1, 3]), "triangle inequality"),
    ([[2, 0.1, 0], [0, 2, 0], [0, 0, 2]], "symmetric"),
    (np.full((3, 3), np.nan), "finite"),
    ([1, 1, 1], "shape"),
])
def test_invalid_inertia(mesh, inertia, message):
    with pytest.raises(ValueError, match=message):
        RigidBody(mesh, 1, [0, 0, 0], inertia)


def test_inertia_accepts_off_diagonal_terms_and_planar_limit(mesh):
    body = RigidBody(mesh, 1, [0, 0, 0], [[5, -1, 0], [-1, 5, 0], [0, 0, 4]])
    assert body.inertia_cg[0, 1] == -1
    RigidBody(mesh, 1, [0, 0, 0], np.diag([1, 1, 2]))


@pytest.mark.parametrize(("field", "value"), [
    ("amplitudes", [-0.1]), ("amplitudes", [np.nan]),
    ("angular_frequencies", [0]), ("angular_frequencies", [-1]),
    ("angular_frequencies", [np.inf]), ("directions", []),
    ("directions", [[0]]), ("phases", [np.inf]), ("phases", ["0"]),
    ("seed", -1), ("seed", True), ("seed", 1.5),
    ("provenance", {"station": 123}), ("provenance", "unstructured"),
])
def test_invalid_wave_inputs(field, value):
    values = dict(amplitudes=[0.1], angular_frequencies=[1], directions=[0], phases=[0])
    values[field] = value
    with pytest.raises(ValueError, match=field):
        WaveComponents(**values)


def test_calm_water_and_zero_amplitude_components():
    assert len(WaveComponents([], [], [], []).amplitudes) == 0
    waves = WaveComponents([0], [0.01], [7 * np.pi], [-9 * np.pi], seed=np.int64(2))
    assert waves.seed == 2
    assert waves.directions[0] == 7 * np.pi
    assert waves.phases[0] == -9 * np.pi


def test_validated_inputs_own_their_data(mesh):
    vertices = mesh.vertices.copy()
    faces = mesh.faces.copy()
    owned = HullMesh(vertices, faces)
    vertices[:] = 0
    faces[:] = 0
    assert_array_equal(owned.vertices, mesh.vertices)
    assert_array_equal(owned.faces, mesh.faces)
    cg = np.zeros(3)
    body = RigidBody(mesh, 1, cg, np.eye(3))
    cg[:] = 20
    assert_array_equal(body.cg, [0, 0, 0])
    amplitudes = np.array([0.1])
    provenance = {"source": "synthetic"}
    waves = WaveComponents(amplitudes, [1], [0], [0], provenance=provenance)
    amplitudes[:] = 4
    provenance["source"] = "changed"
    assert waves.amplitudes[0] == 0.1
    assert waves.provenance["source"] == "synthetic"
    for value in (mesh.vertices, mesh.faces, body.cg, body.inertia_cg, waves.phases,
                  InitialState().orientation):
        with pytest.raises(ValueError):
            value.flat[0] = 4
    with pytest.raises(TypeError):
        waves.provenance["source"] = "changed"


@pytest.mark.parametrize("field", ["position", "linear_velocity", "angular_velocity"])
def test_invalid_initial_vectors(field):
    with pytest.raises(ValueError, match=field):
        InitialState(**{field: [0, np.nan, 0]})


def test_default_state_and_quaternion_roundoff():
    state = InitialState()
    assert_array_equal(state.orientation, [1, 0, 0, 0])
    assert_array_equal(state.linear_velocity, [0, 0, 0])
    assert_array_equal(state.angular_velocity, [0, 0, 0])
    state = InitialState(orientation=[1 + 1e-12, 0, 0, 0])
    assert_array_equal(state.orientation, [1, 0, 0, 0])


def test_input_fixture_round_trip(input_data):
    case = SimulationInputs.from_json(json.dumps(input_data))
    encoded = case.to_json()
    assert json.loads(encoded) == input_data
    decoded = SimulationInputs.from_json(encoded)
    assert_array_equal(decoded.body.mesh.faces, case.body.mesh.faces)
    assert_array_equal(decoded.waves.phases, case.waves.phases)
    assert decoded.waves.seed == 42
    assert dict(decoded.waves.provenance) == input_data["waves"]["provenance"]


def test_calm_water_and_nonidentity_state_round_trip(input_data):
    case = SimulationInputs.from_json(json.dumps(input_data))
    state = InitialState([1, 2, 3], [0.5, 0.5, 0.5, 0.5], [1, 0, 0], [0, 2, 0])
    case = SimulationInputs(case.body, WaveComponents([], [], [], []), state)
    decoded = SimulationInputs.from_json(case.to_json())
    assert_array_equal(decoded.initial_state.orientation, state.orientation)
    assert_array_equal(decoded.initial_state.position, state.position)
    assert_array_equal(decoded.initial_state.linear_velocity, state.linear_velocity)
    assert_array_equal(decoded.initial_state.angular_velocity, state.angular_velocity)
    assert decoded.waves.seed is None
    assert decoded.waves.amplitudes.size == 0


@pytest.mark.parametrize("version", [0, 1, True, 2.0, "2"])
def test_reject_unsupported_input_version(input_data, version):
    input_data["format_version"] = version
    with pytest.raises(ValueError, match="format_version"):
        SimulationInputs.from_json(json.dumps(input_data))


@pytest.mark.parametrize("text", ["{", "null", "[]", "{}"])
def test_reject_malformed_json(text):
    with pytest.raises(ValueError, match="Invalid simulation inputs"):
        SimulationInputs.from_json(text)


def test_json_cannot_bypass_validation(input_data):
    input_data["body"]["mass"] = -1
    with pytest.raises(ValueError, match="mass"):
        SimulationInputs.from_json(json.dumps(input_data))


def test_json_rejects_unknown_fields(input_data):
    input_data["initial_state"]["heading"] = 0
    with pytest.raises(ValueError, match="heading"):
        SimulationInputs.from_json(json.dumps(input_data))
