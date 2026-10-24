"""Bounded batched recovery audit; original scalar modules/attempt remain intact."""
import argparse
import hashlib
import json
import os
import shutil
import sys
import time as wall_clock
from pathlib import Path

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[variable] = '1'

import numpy as np
from scipy.io import loadmat
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from actual_wave_geometry_core import Budget, BudgetExceeded, rim_bounds
from actual_wave_geometry import read, sha, physical_hull, fixtures as scalar_fixtures
from actual_wave_geometry_batch_core import BatchSea, batch_oracle_fixtures, wetted_bounds_batch


def prior_interval_overlap(rows, prior):
    counts = dict(wet_area=0, material_XOR=0, rim=0)
    unknown_rim = []
    for current, old in zip(rows, prior):
        assert current['index'] == old['index'] and current['time_s'] == old['time_s']
        pairs = [('wet_area', current['wetting']['wet_area_bounds_m2'], old['wetting']['wet_area_bounds_m2']),
                 ('material_XOR', current['wetting']['wet_set_turnover_bounds_m2'], old['wetting']['wet_set_turnover_bounds_m2'])]
        if 'clearance_min_bounds_m' in current['rim'] and 'clearance_min_bounds_m' in old['rim']:
            pairs.append(('rim', current['rim']['clearance_min_bounds_m'], old['rim']['clearance_min_bounds_m']))
        else:
            unknown_rim.append(current['index'])
        for quantity, new_interval, old_interval in pairs:
            if max(new_interval[0], old_interval[0]) > min(new_interval[1], old_interval[1]):
                raise ValueError(f'Disjoint preserved-prefix {quantity} bounds at index{current["index"]}')
            counts[quantity] += 1
    return dict(matched_samples=min(len(rows), len(prior)), available_old_samples=len(prior),
                compared_interval_counts=counts, rim_without_two_available_intervals=unknown_rim,
                all_compared_intervals_overlap=True, containment_required=False,
                comparison_status='PASS_AVAILABLE_MATCHING_INTERVALS' if counts['wet_area'] else 'UNKNOWN_NO_MATCHED_OUTPUTS',
                scope='Only available matching intervals; missing coverage/rim bounds remain unknown')


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
    assert config['prior_prefix_geometry'] in config['input_sha256'], 'Old prefix must be frozen'
    prior = read(config['prior_prefix_geometry'])
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
    budget = Budget(600, 10_000_000)
    sea = BatchSea(heading['boat_frame_waves'], case['water']['gravity_m_s2'], budget)
    vertices, faces, rim = physical_hull(case)
    rows = []
    report = dict(complete=False, scope='Actual saved-time material-wetting geometry only; no loads or nonlinear trajectory',
        motion_sha256=sha(config['motion_mat']), config_sha256=sha(args.config), command=[sys.executable, *sys.argv],
        input_sha256=config['input_sha256'], saved_samples=len(time), processed_samples=0,
        source_sha256=sha(__file__), cpu_threads=1,
        spatial_bounds='Same analytic Hessian/set math; batched clipping, depth8,10M triangle cap;1%/10mm reporting precision only',
        temporal_limit='Saved times only; no continuous-time guarantee between native outputs')

    stream = (args.output/'geometry-results.jsonl').open('w', encoding='utf-8')

    def save(final=False):
        stream.flush()
        report.update(processed_samples=len(rows),
                      unprocessed_index_range=[len(rows), len(time)-1] if len(rows) < len(time) else [],
                      triangle_evaluations=budget.triangles, wave_point_evaluations=budget.points,
                      elapsed_budget_seconds=wall_clock.perf_counter()-budget.started, stage_timings_seconds=sea.timings,
                      timing_scope='field and batched clipping/adaptation; scalar fixtures/rim overhead outside these counters')
        if final:
            report['unprocessed_indices'] = list(range(len(rows), len(time)))
            (args.output/'geometry-results.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
            report['elapsed_including_final_geometry_reporting_seconds'] = wall_clock.perf_counter()-budget.started
            report['reporting_layout'] = 'Append-only JSONL; compact checkpoints; one final full JSON/list; final elapsed includes geometry serialization'
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')

    save()
    try:
        with threadpool_limits(1):
            report['fixtures'] = scalar_fixtures(case, vertices, faces, heading, global_bank, sea, budget)
            report['batch_oracle_fixtures'] = batch_oracle_fixtures(vertices[faces])
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
                    wet = wetted_bounds_batch(vertices[faces], posed[faces], sea, float(t), plane, tolerance, 8)
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
                    stream.write(json.dumps(row, allow_nan=False)+'\n')
                except BudgetExceeded:
                    report['budget_exhausted'] = True
                    break
                if index % 20 == 0:
                    save()
            report.update(complete=len(rows) == len(time), all_saved_times_covered=len(rows) == len(time),
                all_spatial_targets_met=bool(rows) and all(r['wetting']['area_target_met'] and r['rim'].get('numerical_target_met', False) for r in rows))
            report['preserved_prefix_overlap'] = prior_interval_overlap(rows, prior)
        assert all(sha(path) == digest for path, digest in config['input_sha256'].items())
    except BaseException as error:
        report.update(error=repr(error), complete=False)
    save(final=True)
    stream.close()
    print(json.dumps({k: report.get(k) for k in ('complete', 'processed_samples', 'all_spatial_targets_met', 'error')}))
    # Valid incomplete bounds are a diagnostic outcome, never a physical PASS.
    return 1 if 'error' in report else 0


if __name__ == '__main__':
    raise SystemExit(main())
