"""Audit actual native poses against prescribed incident geometry; no dynamics."""
import argparse
import hashlib
import json
import os
import shutil
import sys
from pathlib import Path

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from actual_wave_geometry_core import (AnalyticSea, Budget, BudgetExceeded, area_bounds,
    clip, polygon_area, rim_bounds, wetted_bounds)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def physical_hull(case):
    hull = case['minimum_hull']
    q = hull['equilibrium_quaternion_wxyz']
    rotation = Rotation.from_quat([*q[1:], q[0]]).as_matrix()
    vertices = (np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@rotation.T+case['world_cg_m']
    faces = np.asarray([face for i, face in enumerate(hull['triangular_faces'])
                        if i not in hull['nonphysical_closure_faces']], dtype=int)
    edges = {}
    for face in faces:
        for a, b in zip(face, np.roll(face, -1)):
            key = tuple(sorted((int(a), int(b))))
            edges[key] = edges.get(key, 0)+1
    assert all(count in (1, 2) for count in edges.values())
    rim = np.array([edge for edge, count in edges.items() if count == 1])
    assert len(rim) > 0
    return vertices, faces, rim


def fixtures(case, vertices, faces, heading, global_bank, sea, budget):
    # Exact flat-plane triangle clip: submerged corner has half-edge intercepts,
    # so wet area is one quarter of the complete tilted triangle area.
    triangle = np.array([[0., 0., -.5], [1., 0., .5], [0., 1., .5]])
    expected = polygon_area(triangle)/4
    actual = area_bounds(triangle, triangle, triangle[:, 2], 0.)
    np.testing.assert_allclose(actual[:2], expected, rtol=0, atol=1e-14)
    np.testing.assert_allclose(actual[2:], 0, rtol=0, atol=1e-14)
    mean_area = sum(polygon_area(clip(np.column_stack((tri, tri[:, 2], tri[:, 2])), 3, 0)[:, :3])
                    for tri in vertices[faces])
    native_triangles = np.asarray(case['minimum_wetted_triangles_world_m'])
    expected_mean = sum(polygon_area(tri) for tri in native_triangles)
    np.testing.assert_allclose(mean_area, expected_mean, rtol=2e-13, atol=2e-12)
    # Reflection symmetry of this physical boat's flat-plane material wet area.
    cg = np.asarray(case['world_cg_m'])
    values = []
    for heel in (-2., 2.):
        posed = (vertices-cg)@Rotation.from_euler('x', heel, degrees=True).as_matrix().T+cg
        values.append(sum(polygon_area(clip(np.column_stack((a, b[:, 2], a[:, 2])), 3, 0)[:, :3])
                          for a, b in zip(vertices[faces], posed[faces])))
    np.testing.assert_allclose(values[0], values[1], rtol=0, atol=2e-12)
    points = np.array([[0., 0., 0.], [2., -1., 0.], [-3., 2., 0.]])
    rotation = np.asarray(heading['global_from_boat_rotation'])
    translation = np.asarray(heading['global_from_boat_translation_m'])
    global_sea = AnalyticSea(global_bank['waves'], case['water']['gravity_m_s2'], budget)
    errors = []
    for t in (0., .37, 3.1):
        errors.append(np.max(abs(sea.field(points[:, :2], t)[:, 0]
                                 -global_sea.field((points@rotation.T+translation)[:, :2], t)[:, 0])))
    assert max(errors) < 2e-12
    return dict(exact_flat_plane_clip=True, mean_area_matches_independent_native_surface=True,
        flat_water_reflection_parity=True, analytic_bank_frame_closure_max_m=float(max(errors)),
        mean_wet_area_m2=mean_area)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = read(args.config)
    if config['status'] != 'ready_for_coordinator_execution_allocation':
        raise ValueError('Reviewed source and actual-motion identity binding required')
    for path, digest in config['input_sha256'].items():
        assert sha(path) == digest, path
    assert config['motion_mat'] in config['input_sha256'], 'Actual motion must be frozen'
    case, heading, global_bank = (read(config[key]) for key in ('case', 'heading_file', 'global_bank'))
    assert heading['observation_id'] == 'ndbc-46086-20250105T1140Z' and heading['relative_heading_offset_deg'] == 90
    actual = loadmat(config['motion_mat'])
    time = actual['time'].ravel()
    position, q = actual['position_world'], actual['q']
    cg = np.asarray(case['world_cg_m'])
    assert position.shape == q.shape == (len(time), 6) and len(time) == 3001
    assert np.isfinite(time).all() and np.isfinite(position).all() and np.isfinite(q).all()
    assert np.all(np.diff(time) > 0) and abs(time[0]) < 1e-12 and abs(time[-1]-60) < 1e-12
    np.testing.assert_allclose(time, np.arange(3001)*.02, rtol=0, atol=1e-9)
    if 'native_body' in actual:
        native_time = np.asarray(actual['native_body']['time'].item()).ravel()
        assert np.array_equal(native_time, time), 'Actual raw native time identity differs'
    np.testing.assert_allclose(position[:, :3]-cg, q[:, :3], rtol=0, atol=2e-12)
    np.testing.assert_allclose(position[:, 3:], q[:, 3:], rtol=0, atol=2e-12)
    np.testing.assert_allclose(actual['mean_world_cg'].reshape(-1, 3),
                               np.tile(cg, (actual['mean_world_cg'].size//3, 1)), rtol=0, atol=2e-12)
    args.output.mkdir(parents=True, exist_ok=False)
    snapshot = args.output/'source'
    snapshot.mkdir()
    for i, path in enumerate(dict.fromkeys([str(Path(__file__)), str(args.config), *config['input_sha256']])):
        shutil.copyfile(path, snapshot/f'{i:02d}-{Path(path).name}')
    budget = Budget(300, 2_000_000)
    sea = AnalyticSea(heading['boat_frame_waves'], case['water']['gravity_m_s2'], budget)
    vertices, faces, rim = physical_hull(case)
    rows = []
    report = dict(complete=False, scope='Actual saved-time material-wetting geometry only; no loads or nonlinear trajectory',
        motion_sha256=sha(config['motion_mat']), config_sha256=sha(args.config), command=[sys.executable, *sys.argv],
        input_sha256=config['input_sha256'], saved_samples=len(time), processed_samples=0,
        unprocessed_indices=list(range(len(time))), source_sha256=sha(__file__), cpu_threads=1,
        spatial_bounds='Analytic spectral Hessian bound; clipping lower/upper sets, depth8,2M triangle cap',
        temporal_limit='Saved times only; no continuous-time guarantee between native outputs')

    def save():
        report.update(processed_samples=len(rows), unprocessed_indices=list(range(len(rows), len(time))),
                      triangle_evaluations=budget.triangles, wave_point_evaluations=budget.points)
        (args.output/'geometry-results.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        with threadpool_limits(1):
            report['fixtures'] = fixtures(case, vertices, faces, heading, global_bank, sea, budget)
            mean_area = report['fixtures']['mean_wet_area_m2']
            tolerance = config['relative_area_interval_target']*mean_area
            for index, t in enumerate(time):
                try:
                    budget.check()
                    R = Rotation.from_euler('xyz', position[index, 3:]).as_matrix()
                    posed = (vertices-cg)@R.T+position[index, :3]
                    data = sea.field(position[index, :2][None, :], float(t))[0]
                    # Plane eta=etaCG+(xy-CGxy).gradient at actual CG.
                    plane = np.r_[data[0], position[index, :2], data[1:]]
                    normal = np.r_[-data[1:], 1.]
                    normal /= np.linalg.norm(normal)
                    heel_water = -np.arcsin(normal[1])
                    pitch_water = np.arctan2(normal[0], normal[2])
                    water_R = Rotation.from_euler('xyz', [heel_water, pitch_water, 0.]).as_matrix()
                    relative = Rotation.from_matrix(water_R.T@R).as_euler('xyz')
                    wet = wetted_bounds(vertices[faces], posed[faces], sea, float(t), plane, tolerance, 8)
                    try:
                        rim_result = rim_bounds(posed[rim], sea, float(t), config['rim_interval_target_m'], 8)
                    except BudgetExceeded:
                        rim_result = dict(opening_status='UNKNOWN_BUDGET_BEFORE_RIM_BOUNDS')
                    row = dict(index=index, time_s=float(t), native_CG_m=position[index, :3].tolist(),
                        native_angles_rad=position[index, 3:].tolist(), heave_relative_to_local_eta_m=float(position[index, 2]-data[0]-cg[2]),
                        body_vs_water_relative_euler_rad=relative.tolist(), exact_local_eta_and_gradient=data.tolist(),
                        wetting=wet, rim=rim_result)
                    row['wet_set_turnover_fraction_bounds'] = [v/mean_area for v in wet['wet_set_turnover_bounds_m2']]
                    rows.append(row)
                except BudgetExceeded:
                    report['budget_exhausted'] = True
                    break
                if index % 20 == 0:
                    save()
            report.update(complete=len(rows) == len(time), all_saved_times_covered=len(rows) == len(time),
                all_spatial_targets_met=bool(rows) and all(r['wetting']['area_target_met'] and r['rim'].get('numerical_target_met', False) for r in rows))
        assert all(sha(path) == digest for path, digest in config['input_sha256'].items())
    except BaseException as error:
        report['error'] = repr(error)
    save()
    print(json.dumps({k: report.get(k) for k in ('complete', 'processed_samples', 'all_spatial_targets_met', 'error')}))
    # Valid incomplete bounds are a diagnostic outcome, never a physical PASS.
    return 1 if 'error' in report else 0


if __name__ == '__main__':
    raise SystemExit(main())
