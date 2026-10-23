"""Bounded pressure time-derivative audit on the actual nonflat boat snapshot.

Rebuild independent offset meshes and sample identical interior material points.
This tests directional-step convergence, not spatial/ocean/trajectory accuracy.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.reference_boat import BEAM, LENGTH, GRAVITY, WATER_DENSITY, reference_body
from tests.benchmark_rankine import peak_rss_bytes
from .boat_world_boundary import body_material_samples, interpolate_body_trace
from tests.reference_boat import boundary_snapshot
from .boat_world_field import current_field, surface_rates, refine_boundary
from .body_pressure import evaluate_body_pressure, quadrature
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(output, step=1e-4, refinement=0, refinement_region='all', direction_from=None,
        quadrature_order=3):
    from .cuda_operator import prepare_energy, cp
    if not np.isfinite(step) or step <= 0:
        raise ValueError('positive finite directional step required')
    if isinstance(refinement, bool) or not isinstance(refinement, int) or not 0 <= refinement <= 2:
        raise ValueError('bounded audit supports refinement levels zero through two')
    if refinement_region not in ('all', 'body', 'free'):
        raise ValueError('refinement region must be all, body or free')
    if quadrature_order not in (3, 6, 12):
        raise ValueError('bounded audit supports quadrature orders 3, 6 and 12')
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(complete=False, qualified=False, timings=[], derivatives=[], mixed_residuals=[], geometry_checkpoints=[],
        scope='Actual nonflat boat snapshot: pressure directional-step audit in diagnostic finite tank; no spatial/ocean/trajectory qualification',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_pressure_audit.py', 'boat_world_field.py', 'boat_world_boundary.py',
                         'boat_parent_clip_proof.py', 'body_pressure.py', 'canoe.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')},
        gpu_pool_reserved_peak_bytes=0, boundary_refinement=refinement,
        refinement_region=refinement_region,
        quadrature_order=quadrature_order)
    start = perf_counter()
    saved = {}

    def timed(name, function):
        if perf_counter()-start > 120:
            raise RuntimeError('120 s audit budget exhausted between operations')
        before = perf_counter()
        value = function()
        cp.cuda.get_current_stream().synchronize()
        report['gpu_pool_reserved_peak_bytes'] = max(report['gpu_pool_reserved_peak_bytes'],
                                                    cp.get_default_memory_pool().total_bytes())
        report['timings'].append(dict(operation=name, seconds=perf_counter()-before))
        return value

    def checkpoint(name):
        return output.with_name(output.stem+'-'+name+'.json')

    def check_field(field, label):
        residual = field['relative_field_residual']
        report['mixed_residuals'].append(dict(label=label, relative_residual=residual))
        if not np.isfinite(residual) or residual > 1e-9:
            raise ValueError('mixed field residual exceeds 1e-9 algebraic budget')

    def geometry(label, **inputs):
        _, arrays, geometry_report = boundary_snapshot(checkpoint(label), **inputs)
        for level in range(refinement):
            arrays = refine_boundary(arrays, region=refinement_region)
            path = fresh(checkpoint(label+f'-refined{level+1}'))
            np.savez(path.with_suffix('.npz'), **arrays)
            path.write_text(json.dumps(dict(complete=True, qualified=False,
                scope='Conforming subdivision of the checkpointed physical boundary; not an evolved state',
                source_checkpoint=str(checkpoint(label)), boundary_refinement=level+1,
                refinement_region=refinement_region,
                geometry_roundoff_checks_passed=True,
                refinement_moment_errors=arrays['refinement_moment_errors'].tolist(),
                actual_sea_realization=geometry_report['actual_sea_realization'],
                nodes=len(arrays['points']), faces=len(arrays['faces'])), indent=2)+'\n')
        used = checkpoint(label+f'-refined{refinement}') if refinement else checkpoint(label)
        report['geometry_checkpoints'].append(dict(label=label, path=str(used),
            arrays_sha256=sha256(used.with_suffix('.npz').read_bytes()).hexdigest(),
            nodes=len(arrays['points']), faces=len(arrays['faces']),
            face_counts={kind: len(arrays[kind+'_faces']) for kind in ('body', 'free', 'wall', 'bottom')}))
        return arrays

    try:
        body = reference_body()
        arrays = timed('base_geometry', lambda: geometry('base'))
        points, faces = arrays['points'], arrays['faces']
        energy = timed('base_energy', lambda: prepare_energy(faces,
            TrianglePanels(points[faces], 1), order=quadrature_order)).energy
        values = coefficients(sea_states.realize(beam_m=BEAM, gravity=GRAVITY, seed=42))
        twist = np.array([.03, -.02, .01, .01, -.015, .02])
        report['world_twist'] = twist.tolist()
        saved['world_twist'] = twist.copy()
        field = timed('base_field_and_six_responses', lambda: current_field(arrays, energy, values, twist))
        check_field(field, 'base')
        eta_rate, psi_rate, surface = timed('surface_rates', lambda: surface_rates(points,
            faces[arrays['free_faces']], arrays['free_nodes'], field['parent_map'], arrays['parent_psi'], field['flux']))
        saved.update(computed_eta_rate=eta_rate.copy(), computed_psi_rate=psi_rate.copy())
        report['surface_direction'] = dict(kind='self-consistent current discretization')
        if direction_from is not None:
            source_path = Path(direction_from)
            source_report = json.loads(source_path.read_text())
            source_arrays = source_path.with_name(source_path.stem+'-base.npz')
            source_state = dict(np.load(source_arrays, allow_pickle=False))
            source_geometry_path = source_arrays.with_suffix('.json')
            source_geometry = json.loads(source_geometry_path.read_text())
            own_geometry = json.loads(checkpoint('base').read_text())
            source_rates = source_path.with_suffix('.npz')
            if not source_report.get('complete'):
                raise ValueError('prescribed direction requires a complete source audit')
            if not np.array_equal(np.asarray(source_report.get('world_twist')), twist):
                raise ValueError('prescribed direction requires a saved matching world twist; rerun legacy sources')
            if any(source_geometry[key] != own_geometry[key] for key in
                   ('actual_sea_realization', 'body_mass', 'fixture_sha256')):
                raise ValueError('prescribed direction has a different sea or body fixture')
            for key in ('parent_nodes', 'parents', 'parent_psi', 'pose_position', 'pose_orientation',
                        'source_body_vertices', 'source_body_faces'):
                if (arrays[key].shape != source_state[key].shape
                        or not np.allclose(arrays[key], source_state[key], rtol=0, atol=1e-13)):
                    raise ValueError(f'prescribed direction has a different physical state: {key}')
            with np.load(source_rates, allow_pickle=False) as direction:
                imposed_eta, imposed_psi = direction['eta_rate'], direction['psi_rate']
            if (imposed_eta.shape != eta_rate.shape or imposed_psi.shape != psi_rate.shape
                    or not np.isfinite(imposed_eta).all() or not np.isfinite(imposed_psi).all()):
                raise ValueError('prescribed direction does not match finite parent coefficients')
            report['surface_direction'] = dict(kind='prescribed matched-state direction; not a coupled RHS audit',
                source_sha256={str(path): sha256(path.read_bytes()).hexdigest()
                    for path in (source_path, source_arrays, source_rates, source_geometry_path)},
                computed_height_rate_max_difference=float(np.max(abs(eta_rate-imposed_eta))),
                computed_potential_rate_max_difference=float(np.max(abs(psi_rate-imposed_psi))))
            eta_rate, psi_rate = imposed_eta, imposed_psi
        bary, _ = quadrature()
        samples, groups = body_material_samples(arrays, bary)
        rotation = rotation_matrix(arrays['pose_orientation'])
        moving_velocity = twist[:3]+np.cross(twist[3:], points-arrays['pose_position'])
        results = []
        for level, h in enumerate((step, step/2)):
            potentials = []
            for sign in (1, -1):
                dt = sign*h
                orientation = (Rotation.from_rotvec(dt*twist[3:])
                    *Rotation.from_quat(arrays['pose_orientation'], scalar_first=True)).as_quat(scalar_first=True)
                pose = InitialState(arrays['pose_position']+dt*twist[:3], orientation)
                label = f'level{level}-sign{sign}'
                moved = timed(label+'-geometry', lambda: geometry(label, body=body, pose=pose,
                    parent_xy=arrays['parent_nodes'][:, :2], parents=arrays['parents'],
                    parent_eta=arrays['parent_nodes'][:, 2]+dt*eta_rate, parent_psi=arrays['parent_psi']+dt*psi_rate))
                moved_energy = timed(label+'-energy', lambda: prepare_energy(moved['faces'],
                    TrianglePanels(moved['points'][moved['faces']], 1), order=quadrature_order)).energy
                moved_field = timed(label+'-single_field', lambda: current_field(moved,
                    moved_energy, values, twist, acceleration_columns=False, time=dt))
                check_field(moved_field, label)
                potentials.append(timed(label+'-material_trace', lambda: interpolate_body_trace(moved,
                    moved_field['potential'], samples, groups)))
                saved[f'material_trace_{level}_{sign}'] = potentials[-1]
            rate = (potentials[0]-potentials[1])/(2*h)
            result = timed(f'level{level}-pressure', lambda: evaluate_body_pressure(points, faces,
                arrays['body_faces'], arrays['pose_position'], body.mass, rotation@body.inertia_cg@rotation.T,
                twist, field['potential'], rate, moving_velocity, field['acceleration_potentials'],
                density=WATER_DENSITY, gravity=GRAVITY))
            results.append(result)
            report['derivatives'].append(dict(step=h, acceleration_world=result.acceleration.tolist(),
                hydrodynamic_load_world=result.hydrodynamic_load.tolist(),
                known_pressure_load=result.known_pressure_load.tolist(), diagnostics=result.diagnostics))
            saved[f'material_potential_rate_{level}'] = rate
        difference = results[1].acceleration-results[0].acceleration
        load_difference = results[1].hydrodynamic_load-results[0].hydrodynamic_load
        scales = body.mass*GRAVITY*np.r_[np.ones(3), np.full(3, LENGTH)]
        # Reserve 0.1% of weight/weight*length for directional differentiation;
        # this is a numerical audit budget, not a physical accuracy certificate.
        step_error = float(np.max(abs(load_difference)/scales))
        report.update(complete=True, seconds=perf_counter()-start, surface=surface,
            surface_diagnostics_scope='Computed current-field projection and rates; prescribed offsets may use different saved rates',
            acceleration_step_difference_world=difference.tolist(),
            pressure_load_step_difference_world=load_difference.tolist(),
            directional_step_check=dict(normalized_load_max_difference=step_error,
                force_scale_newtons=body.mass*GRAVITY, torque_scale_newton_metres=body.mass*GRAVITY*LENGTH,
                budget=1e-3, passed=step_error <= 1e-3),
            added_mass_reciprocity_relative=report['derivatives'][-1]['diagnostics']['added_mass_reciprocity_relative'])
        saved.update(material_samples=samples, source_groups=groups, eta_rate=eta_rate, psi_rate=psi_rate,
                     potential=field['potential'], acceleration_potentials=field['acceleration_potentials'])
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['process_peak_rss_bytes'] = peak_rss_bytes()
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--step', type=float, default=1e-4)
    parser.add_argument('--refinement', type=int, default=0)
    parser.add_argument('--refinement-region', choices=('all', 'body', 'free'), default='all')
    parser.add_argument('--direction-from', help='Complete audit JSON supplying the same prescribed parent rates')
    parser.add_argument('--quadrature-order', type=int, choices=(3, 6, 12), default=3)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
