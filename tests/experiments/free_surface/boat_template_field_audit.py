"""Compare guarded geometry with a frozen exact snapshot through current BEM.

This is a roundoff-reproduction audit, not a spatial or trajectory validation.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import InitialState, sea_states
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.reference_boat import GRAVITY, WATER_DENSITY, reference_body
from .boat_geometry_template import evaluate_template
from .boat_world_boundary import interpolate_body_trace
from .boat_world_field import current_field, surface_rates
from .body_pressure import evaluate_body_pressure
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def run(input_stem, template, output):
    from .cuda_operator import prepare_energy, cp
    stem, template, output = Path(input_stem), Path(template), fresh(output)
    report = dict(complete=False, qualified=False,
        scope='Current field and pressure sensitivity to guarded geometry roundoff; diagnostic tank only',
        timings=[], input_sha256={}, source_sha256={name: sha256(
            Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('boat_template_field_audit.py', 'boat_geometry_template.py', 'boat_world_field.py',
             'body_pressure.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh')})
    saved = {}

    def load(path):
        report['input_sha256'][str(path)] = sha256(path.read_bytes()).hexdigest()
        return dict(np.load(path, allow_pickle=False))

    def timed(name, function):
        cp.cuda.get_current_stream().synchronize()
        start = perf_counter()
        result = function()
        cp.cuda.get_current_stream().synchronize()
        report['timings'].append(dict(operation=name, seconds=perf_counter()-start))
        return result

    try:
        control = load(stem.with_suffix('.npz'))
        metadata_path = stem.with_suffix('.json')
        metadata = json.loads(metadata_path.read_text())
        report['input_sha256'][str(metadata_path)] = sha256(metadata_path.read_bytes()).hexdigest()
        # Older coarse audits predate the explicit refinement metadata.
        if not metadata['complete'] or metadata.get('boundary_refinement', 0) != 0:
            raise ValueError('complete unrefined frozen pressure audit required')
        base = load(Path(str(stem)+'-base.npz'))
        if len(base['points']) != len(control['potential']):
            raise ValueError('saved base potential does not match coarse checkpoint')
        exact = load(Path(str(stem)+'-level1-sign1.npz'))
        reference = load(template)
        current, diagnostics = timed('guarded_geometry', lambda: evaluate_template(reference,
            pose=InitialState(exact['pose_position'], exact['pose_orientation']),
            parent_eta=exact['parent_nodes'][:, 2], parent_psi=exact['parent_psi']))
        if not np.array_equal(current['faces'], exact['faces']):
            raise ValueError('frozen and guarded charts differ')
        report['geometry'] = dict(diagnostics=diagnostics,
            maximum_coordinate_difference_m=float(np.max(abs(current['points']-exact['points']))))
        sea_path = Path(str(stem)+'-base.json')
        report['input_sha256'][str(sea_path)] = sha256(sea_path.read_bytes()).hexdigest()
        metadata_sea = json.loads(sea_path.read_text())['actual_sea_realization']
        waves = sea_states.loads_realization(json.dumps(metadata_sea))
        report['actual_sea_realization'] = metadata_sea
        values = coefficients(waves)
        twist = np.array([.03, -.02, .01, .01, -.015, .02])
        step = metadata['derivatives'][1]['step']
        traces, fields, rates = [], [], []
        for label, arrays in (('exact', exact), ('guarded', current)):
            energy = timed(label+'-current_energy', lambda: prepare_energy(arrays['faces'],
                TrianglePanels(arrays['points'][arrays['faces']], 1))).energy
            field = timed(label+'-current_field', lambda: current_field(arrays, energy, values,
                twist, acceleration_columns=False, time=step))
            trace = timed(label+'-material_trace', lambda: interpolate_body_trace(arrays,
                field['potential'], control['material_samples'], control['source_groups']))
            try:
                rate = timed(label+'-surface_rates', lambda: surface_rates(arrays['points'],
                    arrays['faces'][arrays['free_faces']], arrays['free_nodes'], field['parent_map'],
                    arrays['parent_psi'], field['flux']))
            except ValueError as error:
                # Offset snapshots are valid mixed-field pressure samples, but
                # their cut-mesh surface mass need not support time evolution.
                # Preserve this independent failure instead of relaxing it.
                report[label+'_surface_rate_failure'] = str(error)
                rate = None
            fields.append(field); traces.append(trace); rates.append(rate)
            saved[label+'_potential'] = field['potential']
            saved[label+'_material_trace'] = trace
            saved[label+'_flux'] = field['flux']
            saved[label+'_parent_map'] = field['parent_map']
            if rate is not None:
                saved[label+'_eta_rate'], saved[label+'_psi_rate'] = rate[:2]
                report[label+'_surface'] = rate[2]
            report[label+'_field_residual'] = field['relative_field_residual']
        report['differences'] = dict(
            potential_max=float(np.max(abs(fields[1]['potential']-fields[0]['potential']))),
            material_trace_max=float(np.max(abs(traces[1]-traces[0]))),
            exact_trace_replay_max=float(np.max(abs(traces[0]-control['material_trace_1_1']))))
        if all(rate is not None for rate in rates):
            report['differences'].update(
                eta_rate_max_m_per_s=float(np.max(abs(rates[1][0]-rates[0][0]))),
                psi_rate_max=float(np.max(abs(rates[1][1]-rates[0][1]))))
        body = reference_body()
        rotation = rotation_matrix(base['pose_orientation'])
        velocity = twist[:3]+np.cross(twist[3:], base['points']-base['pose_position'])
        results = []
        for trace in traces:
            derivative = (trace-control['material_trace_1_-1'])/(2*step)
            results.append(evaluate_body_pressure(base['points'], base['faces'], base['body_faces'],
                base['pose_position'], body.mass, rotation@body.inertia_cg@rotation.T, twist,
                control['potential'], derivative, velocity, control['acceleration_potentials'],
                density=WATER_DENSITY, gravity=GRAVITY))
        report['pressure_load_difference_world'] = (
            results[1].hydrodynamic_load-results[0].hydrodynamic_load).tolist()
        report['acceleration_difference_world'] = (results[1].acceleration-results[0].acceleration).tolist()
        report['exact_load_replay_difference_world'] = (results[0].hydrodynamic_load
            -metadata['derivatives'][1]['hydrodynamic_load_world']).tolist()
        # Numerical reproduction budgets only; physical force convergence is
        # separately unresolved and cannot be inferred from this comparison.
        report['roundoff_reproduction_passed'] = bool(
            report['differences']['potential_max'] < 1e-9
            and report['differences']['exact_trace_replay_max'] < 1e-9
            and max(abs(np.array(report['pressure_load_difference_world']))) < 1e-4
            and max(abs(np.array(report['exact_load_replay_difference_world']))) < 1e-4)
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-stem', required=True)
    parser.add_argument('--template', required=True)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(**vars(parser.parse_args())), indent=2))
