"""Convex geometry comes from input faces, including closure ownership."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from topwater import HullMesh, InitialState, RigidBody, WaveComponents
from topwater.hydrostatics import FloodingError
from .reference_boat import GROUPS, equilibrium, reference_body, simulation_inputs
from .experiments.free_surface.boat_parent_clip_proof import coplanar_groups
from .experiments.free_surface.boat_world_boundary import (
    construct, body_material_samples, interpolate_body_trace,
)
from .experiments.free_surface.boat_geometry_template import evaluate_template


def _box(opened):
    vertices = np.array([[-1,-.5,0], [1,-.5,0], [1,.5,0], [-1,.5,0],
                         [-1,-.5,1], [1,-.5,1], [1,.5,1], [-1,.5,1]], dtype=float)
    faces = np.array([[0,2,1], [0,3,2], [0,1,5], [0,5,4], [1,2,6], [1,6,5],
                      [2,3,7], [2,7,6], [3,0,4], [3,4,7], [4,5,6], [4,6,7]])
    return HullMesh(vertices, faces, [4,5,6,7] if opened else (), [10,11] if opened else ())


def _body(mesh):
    return RigidBody(mesh, 100., [0.,0.,.5], np.diag([40.,50.,60.]))


def _calm():
    return WaveComponents([], [], [], [])


def test_exact_derived_groups_preserve_canonical_order_and_input_twist():
    body = reference_body()
    groups, count = coplanar_groups(body.mesh.vertices, body.mesh.faces, body.mesh.artificial_faces)
    assert groups == GROUPS and count == len(GROUPS)-1
    inputs = simulation_inputs()
    from topwater.frames import rotation_matrix
    assert_allclose(rotation_matrix(inputs.initial_state.orientation)@inputs.initial_state.angular_velocity,
                    [.01,-.015,.02], atol=1e-17, rtol=0)
    assert_allclose(inputs.initial_state.linear_velocity, [.03,-.02,.01], atol=0, rtol=0)


def test_second_boat_fixture_is_explicit_distinct_serializable_physical_input():
    from .reference_convex_boat import box_inputs
    from topwater import SimulationInputs
    from topwater.constants import WATER_DENSITY

    inputs = box_inputs()
    copied = SimulationInputs.from_json(inputs.to_json())
    assert copied.body.mesh.faces.shape == (12,3)
    assert copied.body.mass == 2000.
    extents = np.ptp(copied.body.mesh.vertices, axis=0)
    assert_allclose(extents, [6.5,2.4,1.2], atol=0, rtol=0)
    draft = copied.body.cg[2]-copied.initial_state.position[2]
    assert_allclose(WATER_DENSITY*extents[0]*extents[1]*draft, copied.body.mass, atol=1e-12, rtol=0)


@pytest.mark.parametrize('opened', [False, True])
def test_different_convex_hull_closed_or_capped_has_current_replay_and_material_trace(tmp_path, opened):
    body = _body(_box(opened))
    pose = InitialState([.17,.11,.1])
    mesh, arrays, report = construct(tmp_path/'box.json', body=body, pose=pose, waves=_calm(),
                                    export_template=True)
    assert report['complete'] and 'fixture_sha256' not in report
    assert set(report['input_sha256']) == {'body', 'waves'}
    assert int(arrays['source_physical_group_count']) == (5 if opened else 6)
    assert not np.isin(arrays['body_vertex_source_triangle'], body.mesh.artificial_faces).any()
    replay, _ = evaluate_template(arrays, pose=pose, parent_eta=arrays['parent_nodes'][:,2],
                                  parent_psi=arrays['parent_psi'])
    assert_allclose(replay['points'], mesh.vertices, atol=2e-13, rtol=0)
    points, groups = body_material_samples(replay, np.array([[1/3,1/3,1/3]]))
    phi = .2+replay['points']@np.array([.1,-.3,.4])
    # This pose has identity orientation: material -> world is translation-CG.
    expected = .2+(points+pose.position-body.cg)@np.array([.1,-.3,.4])
    assert_allclose(interpolate_body_trace(replay, phi, points, groups), expected, atol=1e-13, rtol=0)


def test_reordered_and_retessellated_canonical_hull_uses_its_own_source_identity(tmp_path):
    original = reference_body()
    faces = original.mesh.faces.copy()
    # Replace one planar quad's diagonal without moving any physical vertex.
    faces[0], faces[1] = [0,1,5], [1,6,5]
    order = np.roll(np.arange(len(faces)), 5)
    artificial = np.flatnonzero(np.isin(order, original.mesh.artificial_faces))
    mesh = HullMesh(original.mesh.vertices, faces[order], original.mesh.opening_rim, artificial)
    body = RigidBody(mesh, original.mass, original.cg, original.inertia_cg)
    pose = equilibrium()[0]
    _, arrays, report = construct(tmp_path/'retessellated.json', body=body, pose=pose, waves=_calm(),
                                  parent_eta=np.zeros(15), parent_psi=np.zeros(15), export_template=True)
    assert report['complete'] and len(body.mesh.faces) == 18
    groups, count = coplanar_groups(mesh.vertices, mesh.faces, artificial)
    assert count == 9 and len(groups) == 10
    assert sorted(arrays['source_group_faces'].tolist()) == list(range(18))
    assert not np.isin(arrays['body_vertex_source_triangle'], artificial).any()


def test_nonconvex_but_topologically_valid_hull_is_rejected_exactly():
    mesh = _box(False)
    vertices = np.vstack((mesh.vertices, [[0.,0.,.7]]))
    faces = np.vstack((mesh.faces[:10], [[4,5,8], [5,6,8], [6,7,8], [7,4,8]]))
    dented = HullMesh(vertices, faces)
    with pytest.raises(ValueError, match='exactly convex'):
        coplanar_groups(dented.vertices, dented.faces)


def test_multiple_nonplanar_artificial_closures_never_supply_pressure_faces(tmp_path):
    box = _box(False)
    vertices = np.vstack((box.vertices, [[-1.,0.,1.2], [1.,0.,1.2]]))
    caps = np.array([[5,6,9], [4,8,7], [4,5,9], [4,9,8], [8,9,6], [8,6,7]])
    faces = np.vstack((box.faces[:10], caps))
    hull = HullMesh(vertices, faces, [4,5,6,7], np.arange(10,16))
    _, arrays, report = construct(tmp_path/'multiple-caps.json', body=_body(hull),
        pose=InitialState([.17,.11,.1]), waves=_calm(), export_template=True)
    groups, count = coplanar_groups(hull.vertices, hull.faces, hull.artificial_faces)
    assert count == 5 and len(groups)-count == 4
    assert report['complete']
    assert not np.isin(arrays['body_vertex_source_triangle'], hull.artificial_faces).any()


def test_virtual_graph_above_cap_is_rejected_without_claiming_pressure_coverage(tmp_path):
    body = _body(_box(True))
    x, y = np.meshgrid([-2.,0.,2.], [-1.,0.,1.], indexing='ij')
    xy = np.column_stack((x.ravel(), y.ravel()))
    parents = []
    for i in range(2):
        for j in range(2):
            a = i*3+j
            parents.extend(((a,a+3,a+4), (a,a+4,a+1)))
    eta = np.zeros(9); eta[4] = .8
    # Graph at the rim is below .6m, but the arbitrary virtual graph exceeds the
    # artificial cap inside the opening. This implementation rejects that scope.
    with pytest.raises(ValueError, match='every artificial closure dry') as error:
        construct(tmp_path/'unsupported-cap.json', body=body,
            pose=InitialState([.17,.11,.1]), waves=_calm(), parent_xy=xy, parents=parents,
            parent_eta=eta, parent_psi=np.zeros(9))
    assert not isinstance(error.value, FloodingError)
