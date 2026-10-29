"""Standard mass budget, convex physical geometry and CLI publication."""
from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import ConvexHull

from topwater.boat.generation import generate_boat
from topwater.boat.generation import _geometry
from topwater.boat.generation import _structure_moments
from topwater.cli import main
from topwater.cli import parser
from topwater.workflow_data import Boat


LENGTH = 24*.3048
BEAM = 9*.3048
FLAGS = ['--length', str(LENGTH), '--beam', str(BEAM), '--deadrise', '20', '--mass', '2500']


@pytest.fixture(scope='module')
def standard():
    return generate_boat(LENGTH, BEAM, 20, 2500)


def test_generated_geometry_is_the_open_convex_control_envelope(standard):
    boat = standard
    np.testing.assert_allclose(np.ptp(boat.vertices_m, axis=0)[:2], [LENGTH, BEAM], atol=1e-14)
    assert boat.vertices_m[:, 0].min() == boat.vertices_m[:, 2].min() == 0
    mesh = boat.to_body().mesh
    triangles = boat.vertices_m[boat.faces]
    normal = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    normal /= np.linalg.norm(normal, axis=1)[:, None]
    for triangle, n in zip(triangles, normal):
        assert np.max((boat.vertices_m-triangle[0]) @ n) < 1e-12
    stern = np.all(triangles[:, :, 0] == 0, axis=1)
    assert stern.any() and np.all(normal[stern, 0] < -.99)
    edges = np.sort(np.concatenate([boat.faces[:, [0, 1]], boat.faces[:, [1, 2]],
                                    boat.faces[:, [2, 0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    rim = mesh.opening_rim
    assert set(unique[counts == 1].ravel()) == set(rim)
    assert len(mesh.artificial_faces) > 0
    np.testing.assert_array_equal(mesh.faces[mesh.physical_faces], boat.faces)
    points, faces, _, _ = _geometry(BEAM/LENGTH, 20)
    assert ConvexHull(mesh.vertices).volume == pytest.approx(ConvexHull(points).volume*LENGTH**3)


def test_nominal_mass_properties_and_numerical_convergence(standard):
    # Independent planning calculation: 1452.4 kg structure and two 523.8 kg lumps.
    assert .1*(LENGTH/.3048)**3+70 == pytest.approx(1452.4)
    assert standard.mass_kg == 2500
    np.testing.assert_allclose(standard.cg_m, [2.71419579, 0, .76082871], atol=1e-8)
    np.testing.assert_allclose(np.diag(standard.inertia_cg_kg_m2),
                               [1615.52876716, 9523.2623034, 9887.64966325], atol=1e-7)
    assert standard.cg_m[1] == 0
    assert standard.inertia_cg_kg_m2[0, 2] < -1000
    np.testing.assert_array_equal(standard.inertia_cg_kg_m2[[0, 1], [1, 2]], [0, 0])
    points, faces, cap_planes, _ = _geometry(BEAM/LENGTH, 20)
    first, second = _structure_moments(points, faces, cap_planes)
    fine_first, fine_second = _structure_moments(points, faces, cap_planes, stations=401)
    # This is below 1 mm in the assembled nominal CG; the model intentionally
    # retains Marine Gyro's 41-station quadrature rather than exact skin area.
    np.testing.assert_allclose(first, fine_first, atol=2e-4)
    np.testing.assert_allclose(second, fine_second, atol=2e-4)


def test_section_mass_integration_matches_analytic_box():
    # Uniform section-line structure, not a uniform-volume solid: sides+bottom
    # own 70% and a deck at z=.6H owns 30%. x quadrature has a known x^2 error.
    from .test_workflow_data import open_box
    points, physical = open_box()
    points[:, 0] = (points[:, 0]+1)/2
    faces = np.vstack((physical, [[4, 5, 6], [4, 6, 7]]))
    first, second = _structure_moments(points, faces, np.array([[0., 0., 1., -1.]]))
    np.testing.assert_allclose(first, [.5, 0, .7/4+.3*.6], atol=1e-14)
    assert second[1, 1] == pytest.approx(.7*(2/3+2)/4+.3/3)
    assert second[2, 2] == pytest.approx(.7*(2/3)/4+.3*.6**2)
    assert second[0, 2] == pytest.approx(.5*first[2])
    assert abs(second[0, 0]-1/3) < 2e-4


def test_added_total_mass_goes_only_to_the_two_payload_points(standard):
    heavier = generate_boat(LENGTH, BEAM, 20, 2700)
    positions = LENGTH*np.array([[.05, 0, .03], [.50, 0, .20]])
    np.testing.assert_array_equal(heavier.vertices_m, standard.vertices_m)
    np.testing.assert_array_equal(heavier.faces, standard.faces)
    np.testing.assert_allclose(2700*heavier.cg_m-2500*standard.cg_m, 100*positions.sum(axis=0))
    def raw_second(boat):
        I = boat.inertia_cg_kg_m2
        return np.trace(I)/2*np.eye(3)-I+boat.mass_kg*np.outer(boat.cg_m, boat.cg_m)
    np.testing.assert_allclose(raw_second(heavier)-raw_second(standard),
                               100*(positions.T @ positions), atol=1e-10)


def test_public_example_uses_the_documented_generator_inputs():
    project = Path(__file__).resolve().parents[1]
    saved = Boat.load(project/'examples/simple-boat-native/boat.npz')
    generated = generate_boat(7.3152, 2.7432, 20, 2500)
    for name in Boat.members:
        np.testing.assert_array_equal(getattr(saved, name), getattr(generated, name))


@pytest.mark.parametrize('values', [
    (-1, BEAM, 20, 2500), (LENGTH, 0, 20, 2500), (LENGTH, LENGTH, 20, 2500),
    (LENGTH, BEAM, -1, 2500), (LENGTH, BEAM, 60, 2500), (LENGTH, BEAM, 55, 2500),
    (LENGTH, BEAM, 20, 0), (LENGTH, BEAM, 20, 1452.4),
    (np.nan, BEAM, 20, 2500), (LENGTH, np.inf, 20, 2500),
    (LENGTH, BEAM, np.nan, 2500), (LENGTH, BEAM, 20, np.inf),
])
def test_invalid_specifications(values):
    with pytest.raises(ValueError):
        generate_boat(*values)


def test_extrapolation_is_reported():
    with pytest.warns(UserWarning, match='outside the nominal'):
        boat = generate_boat(27*.3048, 9*.3048, 20, 3500)
    assert boat.mass_kg == 3500


def test_cli_defaults_and_archive_publication(tmp_path, monkeypatch):
    assert parser().parse_args(['boat', *FLAGS]).directory == Path('.')
    monkeypatch.chdir(tmp_path)
    assert main(['boat', *FLAGS]) == 0
    first = Boat.load('boat.npz')
    assert main(['boat', *FLAGS]) == 0
    second = Boat.load('boat.npz')
    for name in Boat.members:
        np.testing.assert_array_equal(getattr(first, name), getattr(second, name))
    assert {p.name for p in tmp_path.iterdir()} == {'boat.npz', 'boat.log'}
    log = (tmp_path/'boat.log').read_text()
    assert 'payload_split=50/50' in log and 'status=success' in log
    assert 'inertia_cg_kg_m2=' in log


def test_failed_generation_preserves_previous_archive(tmp_path):
    destination = tmp_path/'boat.npz'
    destination.write_bytes(b'previous boat')
    assert main(['boat', str(tmp_path), *FLAGS[:-1], '100']) == 2
    assert destination.read_bytes() == b'previous boat'
    assert 'output_published=false' in (tmp_path/'boat.log').read_text()


def test_failed_flotation_preserves_previous_archive(tmp_path, monkeypatch):
    from topwater.boat import hydrostatics as boat_hydrostatics
    def failure(body):
        raise ValueError('equilibrium is not stably restoring')
    monkeypatch.setattr(boat_hydrostatics, 'floating_equilibrium', failure)
    (tmp_path/'boat.npz').write_bytes(b'previous boat')
    assert main(['boat', str(tmp_path), *FLAGS]) == 2
    assert (tmp_path/'boat.npz').read_bytes() == b'previous boat'
    assert 'equilibrium is not stably restoring' in (tmp_path/'boat.log').read_text()
