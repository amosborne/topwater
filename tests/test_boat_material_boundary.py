"""Actual boat geometry and immutable bulk cells under six-axis motion."""
import json
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from .reference_boat import BEAM, GRAVITY
from .reference_hull import surface_integrals
from .experiments.free_surface.boat_world_boundary import construct


def _group_moments(arrays):
    triangles = arrays['points'][arrays['faces'][arrays['body_faces']]]
    normal_area = np.cross(triangles[:, 1]-triangles[:, 0],
                           triangles[:, 2]-triangles[:, 0])/2
    area = np.linalg.norm(normal_area, axis=1)
    rows = np.column_stack((area, area[:, None]*triangles.mean(axis=1), normal_area))
    return {group: rows[arrays['body_source_group'] == group].sum(axis=0)
            for group in np.unique(arrays['body_source_group'])}


def test_actual_nonflat_boat_preserves_boundary_and_bulk_material_triangles(tmp_path):
    _, original, _ = construct(tmp_path/'original.json')
    snapshots = []
    for index, time in enumerate((0., .15)):
        orientation = (Rotation.from_rotvec(time*np.array([.01, -.015, .02]))
            *Rotation.from_quat(original['pose_orientation'], scalar_first=True)).as_quat(scalar_first=True)
        pose = InitialState(original['pose_position']+time*np.array([.03, -.02, .01]), orientation)
        waves = sea_states.realize(beam_m=BEAM, gravity=GRAVITY, seed=42)
        xy = original['parent_nodes'][:, :2]
        k = waves.angular_frequencies**2/GRAVITY
        phase = ((xy[:, :1]*np.cos(waves.directions)+xy[:, 1:]*np.sin(waves.directions))*k
                 +waves.phases-time*waves.angular_frequencies)
        eta = np.cos(phase)@waves.amplitudes
        inputs = dict(pose=pose, parent_eta=eta, parent_psi=original['parent_psi'])
        _, reference, _ = construct(tmp_path/f'reference-{index}.json', **inputs)
        _, material, report = construct(tmp_path/f'material-{index}.json',
                                        material_refinement=1, **inputs)
        assert report['material_overlay']['per_group_area_and_first_moments_exact']
        if time:
            assert report['material_overlay']['disconnected_material_cells'] > 0
        for group, expected in _group_moments(reference).items():
            assert_allclose(_group_moments(material)[group], expected, rtol=3e-14, atol=5e-14)
        for actual, expected in zip(surface_integrals(material['points'], material['faces']),
                                    surface_integrals(reference['points'], reference['faces'])):
            assert_allclose(actual, expected, rtol=3e-14, atol=3e-14)
        indices = material['free_parent_indices']
        support = np.maximum(indices, 0)
        inherited = np.einsum('ni,nij->nj', material['free_parent_weights'],
                              material['parent_nodes'][support])
        assert_allclose(inherited, material['points'][material['free_nodes']], atol=2e-15, rtol=0)
        bulk = {}
        for cell in material['fully_wet_material_cells']:
            rows = np.flatnonzero(material['body_material_cell'] == cell)
            assert len(rows) == 1  # A complete cell keeps its single triangle.
            triangle = material['body_material_vertices'][rows[0]]
            bulk[int(cell)] = np.array(sorted(map(tuple, triangle)))
        snapshots.append(bulk)
    shared = snapshots[0].keys() & snapshots[1].keys()
    assert shared
    for cell in shared:
        assert_allclose(snapshots[0][cell], snapshots[1][cell], rtol=0, atol=1e-15)


def test_accepted_boat_motion_retains_nonstar_waterline_fringe(tmp_path):
    case = json.loads(Path(__file__).with_name('data').joinpath('boat-state-200ms.json').read_text())
    y = np.array(case['state'])
    pose = InitialState(y[:3], y[3:7]/np.linalg.norm(y[3:7]))
    inputs = dict(pose=pose, parent_eta=y[13:28], parent_psi=y[28:])
    _, reference, _ = construct(tmp_path/'reference.json', **inputs)
    _, material, report = construct(tmp_path/'material.json', material_refinement=2, **inputs)
    assert report['material_overlay']['nonstar_fringe_loops'] > 0
    for group, expected in _group_moments(reference).items():
        assert_allclose(_group_moments(material)[group], expected, rtol=3e-14, atol=5e-14)
    for actual, expected in zip(surface_integrals(material['points'], material['faces']),
                                surface_integrals(reference['points'], reference['faces'])):
        assert_allclose(actual, expected, rtol=3e-14, atol=3e-14)
