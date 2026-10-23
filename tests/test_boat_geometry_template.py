"""Actual clipping provenance and safe rejection, without a fluid solve."""
import json
from pathlib import Path
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation

from topwater import InitialState
from .reference_boat import boundary_snapshot
from .experiments.free_surface.boat_geometry_template import evaluate_template, InvalidTemplate


@pytest.mark.parametrize('name', ['boat-state-275ms.json', 'boat-state-337ms.json'])
def test_moving_boat_concave_wet_facet_replays_without_changing_coverage(tmp_path, name):
    case = json.loads(Path(__file__).with_name('data').joinpath(name).read_text())
    state = np.array(case['state'])
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    _, arrays, report = boundary_snapshot(tmp_path/'concave.json', pose=pose,
        parent_eta=state[13:28], parent_psi=state[28:], export_template=True)
    if name == 'boat-state-275ms.json':
        assert report['body_polygon_merge']['additional_nonstar_polygon_pieces'] > 0
    current, _ = evaluate_template(arrays, pose=pose,
        parent_eta=state[13:28], parent_psi=state[28:])
    assert_allclose(current['points'], arrays['points'], rtol=0, atol=1e-12)
    # Independent edge count/winding check on the actual closed fluid boundary.
    faces = arrays['faces']
    edges = np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]))
    _, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    assert np.all(counts == 2)
    assert np.all(np.bincount(inverse, weights=np.where(edges[:, 0] < edges[:, 1], 1, -1)) == 0)


def test_cold_template_replays_current_boat_and_rejects_changed_incidence(tmp_path):
    _, arrays, report = boundary_snapshot(tmp_path/'reference.json', export_template=True)
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    current, diagnostics = evaluate_template(arrays, pose=pose,
        parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
    assert_allclose(current['points'], arrays['points'], rtol=0, atol=1e-12)
    assert_allclose(current['body_material_vertices'], arrays['body_material_vertices'], rtol=0, atol=1e-12)
    assert_allclose(current['free_parent_weights'], arrays['free_parent_weights'], rtol=0, atol=1e-12)
    assert diagnostics['displaced_volume_m3'] > 0
    assert np.array_equal(arrays['node_chart_key'], current['node_chart_key'])
    legacy = {key: value for key, value in arrays.items() if key != 'node_chart_key'}
    replay, _ = evaluate_template(legacy, pose=pose,
        parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
    assert np.array_equal(current['node_chart_key'], replay['node_chart_key'])
    changed = dict(arrays, source_body_cg=arrays['source_body_cg']+1e-8)
    with pytest.raises(InvalidTemplate, match='ownership'):
        evaluate_template(changed, pose=pose, parent_eta=arrays['parent_nodes'][:, 2],
                          parent_psi=arrays['parent_psi'])
    # Cold export covers closed-hull/prism vertices above the current roof too.
    owners = len(arrays['template_owner_triples'])
    active_owners = set(arrays['template_node_to_owner'])-{-1}
    assert owners > len(active_owners)
    assert report['template_export']['dry_and_wet_critical_vertices'] > 0
    orientation = (Rotation.from_rotvec([1e-4, 0, 1e-4])
        *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
    with pytest.raises(InvalidTemplate):
        evaluate_template(arrays, pose=InitialState(pose.position, orientation),
            parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
    altered = dict(arrays, faces=arrays['faces'][:, ::-1])
    with pytest.raises(InvalidTemplate, match='reindexed'):
        evaluate_template(altered, pose=pose, parent_eta=arrays['parent_nodes'][:, 2],
                          parent_psi=arrays['parent_psi'])


def test_constraint_provenance_does_not_promote_an_unrelated_incident_plane():
    from .experiments.free_surface.boat_parent_clip_proof import PlanePoints
    registry = PlanePoints(capture_constraints=True)
    definitions = [([(0, 0, 0), (0, 1, 0), (0, 0, 1)], -1),
                   ([(1, 0, 0), (1, 1, 0), (1, 0, 1)], 1),
                   ([(0, 0, 0), (1, 0, 0), (0, 0, 1)], -1),
                   ([(0, 1, 0), (1, 1, 0), (0, 1, 1)], 1),
                   ([(0, 0, 0), (1, 0, 0), (0, 1, 0)], -1),
                   ([(0, 0, 1), (1, 0, 1), (0, 1, 1)], 1)]
    constraints = []
    for points, interior_sign in definitions:
        plane, _ = registry.plane(points, ('test_box', len(constraints)))
        constraints.append((plane, interior_sign))
    extra, _ = registry.plane([(0, 0, 0), (1, -1, 0), (1, 0, -1)], ('unrelated', 0))
    vertices = registry.polytope(constraints)
    assert len(vertices) == 8
    origin = next(i for i in vertices if all(x == 0 for x in registry.points[i]))
    assert extra in registry.incidence[origin]  # Geometric coincidence is real.
    assert extra not in registry.relevant[origin]  # It was NOT a clipping constraint.
    assert registry.relevant[origin] == {plane for plane, _ in constraints}
