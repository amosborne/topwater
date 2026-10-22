"""Actual clipping provenance and safe rejection, without a fluid solve."""
import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation

from topwater import InitialState
from .experiments.free_surface.boat_world_boundary import construct
from .experiments.free_surface.boat_geometry_template import evaluate_template, InvalidTemplate


def test_cold_template_replays_current_boat_and_rejects_changed_incidence(tmp_path):
    _, arrays, report = construct(tmp_path/'reference.json', export_template=True)
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    current, diagnostics = evaluate_template(arrays, pose=pose,
        parent_eta=arrays['parent_nodes'][:, 2], parent_psi=arrays['parent_psi'])
    assert_allclose(current['points'], arrays['points'], rtol=0, atol=1e-12)
    assert_allclose(current['body_material_vertices'], arrays['body_material_vertices'], rtol=0, atol=1e-12)
    assert_allclose(current['free_parent_weights'], arrays['free_parent_weights'], rtol=0, atol=1e-12)
    assert diagnostics['displaced_volume_m3'] > 0
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
