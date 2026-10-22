"""Compare local operator models on actual moving-boat pressure stencils.

Two accepted states train a rank-one geometry map. Optional two-direction models
use the .15 s geometry to choose training directions; .20 s is a later holdout.
The separate past-only fit uses only .05/.10 geometry and pressure stencils.
All current geometry, loads, incident forcing and prescribed surface direction
are shared between models.
Fresh operators/residuals are OFFLINE references, not an online refresh rule.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState, sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import reference_body
from .boat_geometry_template import evaluate_template, InvalidTemplate
from .boat_world_field import current_field, surface_rates, refine_boundary
from .boat_pressure_stencil import evaluate as evaluate_stencil
from .incident_baseline_benchmark import coefficients
from .p4_fluid_ordering import fresh


def _shape_coordinates(points, origin, weight, directions, gram):
    change = points-origin
    alpha = np.linalg.solve(gram, [np.sum(weight[:, None]*g*change) for g in directions])
    residual = change-sum(a*g for a, g in zip(alpha, directions))
    return alpha, float(np.sqrt(np.sum(weight[:, None]*residual*residual)))


def _quadratic_features(alpha, epsilon, *, mixed=False):
    z = np.asarray(alpha)/epsilon
    return np.r_[z, .5*z*z, z[0]*z[1]] if mixed else np.r_[z, .5*z*z]


def _local_models(origin, energy, weight, directions, gram, derivatives, quadratic_coefficients, epsilon):
    def predict(arrays, quadratic):
        alpha, residual = _shape_coordinates(arrays['points'], origin, weight, directions, gram)
        if quadratic:
            result = energy+np.einsum('i,ijk->jk', _quadratic_features(alpha, epsilon,
                mixed=len(quadratic_coefficients) == 5), quadratic_coefficients)
        else:
            result = energy+sum(a*d for a, d in zip(alpha, derivatives))
        return result, dict(coefficients=alpha.tolist(), orthogonal_geometry_rms_m=residual,
            maximum_normalized_coordinate=float(np.max(abs(alpha/epsilon))),
            maximum_feature_magnitude=float(np.max(abs(_quadratic_features(alpha, epsilon, mixed=True)))))
    return lambda arrays: predict(arrays, False), lambda arrays: predict(arrays, True)


def run(pilot, output, operator_from=None, two_directions=False, target_time=.15, tangent_from=None,
        mixed_curvature=False, past_directions=False, past_from=None):
    if target_time not in (.15, .20):
        raise ValueError('bounded target time must be .15 or .20 s')
    if mixed_curvature and tangent_from is None:
        raise ValueError('mixed-curvature extension requires the saved four-sample local model')
    past_enabled = past_directions or past_from is not None
    if past_enabled and (two_directions or tangent_from is not None or mixed_curvature):
        raise ValueError('past-only and geometry-informed models are separate experiments')
    pilot, output = Path(pilot), fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    recorded = json.loads(pilot.read_text())
    if not recorded.get('complete'):
        raise ValueError('requires a complete saved moving pilot')
    with np.load(pilot.with_suffix('.npz'), allow_pickle=False) as data:
        times, states = data['accepted_times'], data['accepted_states']
    selected = []
    requested_times = [.05, .10, .15]+([.20] if target_time == .20 else [])
    for time in requested_times:
        matches = np.flatnonzero(abs(times-time) <= 1e-14)
        if len(matches) != 1:
            raise ValueError('requires unique accepted anchor/training/target states')
        selected.append(int(matches[0]))
    states = states[selected]
    size = (states.shape[1]-13)//2
    report = dict(complete=False, qualified=False,
        scope='One held-out moving-boat geometry/operator/pressure comparison; no ocean, temporal/spatial, cache-policy or throughput qualification',
        anchor_times_seconds=[.05, .10], target_time_seconds=target_time, quadrature_order=12,
        input_sha256={str(p): sha256(p.read_bytes()).hexdigest()
                      for p in (pilot, pilot.with_suffix('.npz'))},
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_operator_prediction.py', 'boat_geometry_template.py',
                         'boat_world_field.py', 'boat_world_boundary.py', 'body_pressure.py',
                         'canoe.py', 'incident_baseline_benchmark.py',
                         'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh', 'boat_pressure_stencil.py')},
        phases=[], models=[])
    if past_enabled:
        for name in ('boat_past_operator.py', 'boat_operator_reuse.py'):
            report['source_sha256'][name] = sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
    previous = None
    if operator_from is not None:
        operator_from = Path(operator_from)
        source = json.loads(operator_from.read_text())
        if (not source.get('complete') or source['input_sha256'] != report['input_sha256']
                or source['target_time_seconds'] != target_time):
            raise ValueError('operator replay requires the same complete moving-state audit')
        for name, expected in source['source_sha256'].items():
            if name != Path(__file__).name and report['source_sha256'][name] != expected:
                raise ValueError(f'operator replay helper changed: {name}')
        with np.load(operator_from.with_suffix('.npz'), allow_pickle=False) as data:
            previous = dict(data)
        report['operator_replay'] = dict(path=str(operator_from),
            sha256=sha256(operator_from.read_bytes()).hexdigest(),
            arrays_sha256=sha256(operator_from.with_suffix('.npz').read_bytes()).hexdigest(),
            scope='Reuse saved reference operators; optional training assemblies are counted separately')
    saved, started = {}, perf_counter()

    def timed(name, function):
        if perf_counter()-started > 120:
            raise RuntimeError('120 s prediction audit budget exhausted between operations')
        before = perf_counter()
        value = function()
        report['phases'].append(dict(operation=name, seconds=perf_counter()-before))
        return value

    def inputs(y):
        return dict(pose=InitialState(y[:3], y[3:7]/np.linalg.norm(y[3:7])),
                    parent_eta=y[13:13+size], parent_psi=y[13+size:])

    try:
        charts = []
        for row in recorded['geometry']:
            path = Path(row['path']).with_suffix('.npz')
            if sha256(path.read_bytes()).hexdigest() != row['sha256']:
                raise ValueError('saved geometry differs from pilot manifest')
            with np.load(path, allow_pickle=False) as data:
                charts.append((path, dict(data)))
        if past_enabled:
            startup = recorded['rhs'][0]
            if startup['time'] != 0. or not startup['complete']:
                raise ValueError('past-only diagnostic requires completed startup chart provenance')
            count = 1+sum(g['mode'] == 'exact chart refresh' for g in startup['geometry'])
            charts = charts[:count]
            report['training_chart_scope'] = 'Only charts available from startup RHS, before both past anchors'
        for chart_path, chart in charts:
            try:
                base = [evaluate_template(chart, **inputs(y))[0] for y in (states[:2] if past_enabled else states)]
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no common guarded geometry chart for the selected states')
        if past_enabled:
            # A future chart failure invalidates the prediction; never select a
            # different training chart using held-out geometry.
            base += [evaluate_template(chart, **inputs(y))[0] for y in states[2:]]
        report['chart'] = dict(path=str(chart_path), sha256=sha256(chart_path.read_bytes()).hexdigest())
        body = reference_body()
        if (not np.array_equal(chart['source_body_vertices'], body.mesh.vertices)
                or not np.array_equal(chart['source_body_faces'], body.mesh.faces)
                or not np.array_equal(chart['source_body_cg'], body.cg)):
            raise ValueError('saved chart differs from current canonical body')
        geometry_metadata = json.loads(chart_path.with_suffix('.json').read_text())
        values = coefficients(sea_states.loads_realization(json.dumps(geometry_metadata['actual_sea_realization'])))

        def refine(arrays):
            for _ in range(2):
                arrays = refine_boundary(arrays, region='body')
            return arrays

        base = [refine(a) for a in base]
        reference_ids = ('faces', 'body_faces', 'free_faces', 'free_nodes', 'body_source_group',
                         'free_source_parent', 'free_parent_indices')

        def same_basis(arrays):
            if any(not np.array_equal(arrays[key], base[0][key]) for key in reference_ids):
                raise ValueError('prediction requires identical chart basis and boundary ownership')

        for arrays in base:
            same_basis(arrays)

        def assemble(name, arrays, key):
            if previous is not None:
                return previous[key]
            from .cuda_operator import prepare_energy
            return timed(name, lambda: prepare_energy(arrays['faces'],
                TrianglePanels(arrays['points'][arrays['faces']], 1), order=12).energy)

        energies = [assemble(f'anchor{i}-energy', arrays, f'anchor{i}_energy')
                    for i, arrays in enumerate(base[:2])]
        target, y = base[-1], states[-1]
        if previous is not None:
            if (not np.array_equal(previous['selected_states'], states)
                    or any(not np.allclose(previous[key], arrays['points'], rtol=0, atol=1e-12)
                           for key, arrays in zip(('anchor0_points', 'anchor1_points', 'target_points'),
                                                  (base[0], base[1], base[-1])))
                    or source['chart'] != report['chart']):
                raise ValueError('operator replay differs from current chart/states/points')
        target_energy = assemble('target-energy', target, 'target_energy')
        twist, pose = y[7:13], inputs(y)['pose']
        actual = timed('target-field', lambda: current_field(target, target_energy, values, twist, time=target_time))
        eta_rate, psi_rate, _ = surface_rates(target['points'], target['faces'][target['free_faces']],
            target['free_nodes'], actual['parent_map'], target['parent_psi'], actual['flux'])
        offsets = {}
        offset_energies = {}
        for h in (.001, .0005):
            for sign in (1, -1):
                dt = sign*h
                orientation = (Rotation.from_rotvec(dt*twist[3:])
                    *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                arrays = timed(f'offset-{dt}-geometry', lambda: refine(evaluate_template(chart,
                    pose=InitialState(pose.position+dt*twist[:3], orientation),
                    parent_eta=y[13:13+size]+dt*eta_rate,
                    parent_psi=y[13+size:]+dt*psi_rate)[0]))
                same_basis(arrays)
                offsets[dt] = arrays
                offset_energies[dt] = assemble(f'offset-{dt}-energy', arrays, f'offset_{dt}_energy')

        # Fixed lumped AREA weights define one explicit geometry projection.
        # Their derivative is zero: never silently change the metric at offsets.
        def nodal_area(arrays):
            panels = TrianglePanels(arrays['points'][arrays['faces']], 1)
            return np.bincount(arrays['faces'].ravel(),
                weights=np.repeat(panels.areas/3, 3), minlength=len(arrays['points']))

        weight = (nodal_area(base[0])+nodal_area(base[1]))/2
        weight /= weight.sum()
        origin, delta = base[0]['points'], base[1]['points']-base[0]['points']

        def inner(a, b):
            return float(np.sum(weight[:, None]*a*b))

        denominator = inner(delta, delta)
        if denominator <= 0:
            raise ValueError('zero anchor geometry displacement')

        def prediction(arrays):
            change = arrays['points']-origin
            alpha = inner(change, delta)/denominator
            return energies[0]+alpha*(energies[1]-energies[0]), dict(alpha=alpha,
                orthogonal_geometry_rms_m=float(np.sqrt(inner(change-alpha*delta, change-alpha*delta))))

        geometry_rate = (offsets[.0005]['points']-offsets[-.0005]['points'])/.001
        alignment = inner(delta, geometry_rate)/np.sqrt(denominator*inner(geometry_rate, geometry_rate))
        report.update(anchor_geometry_rms_m=float(np.sqrt(denominator)),
                      pressure_direction_geometry_cosine=float(alignment),
                      faces=len(target['faces']), nodes=len(target['points']))
        tangent_prediction = quadratic_prediction = None
        extra_local_models = {}
        if two_directions and tangent_from is None:
            training = dict(complete=False, epsilon=.02, anchor_time_seconds=.10,
                scope='Geometry-informed tangent diagnostic: target .15 geometry informs direction one; its operator/force values do not. Not a past-only predictor.',
                assembly_count=0)
            report['two_direction_training'] = training
            training_started = perf_counter()
            try:
                anchor, anchor_state = base[1], states[1]
                anchor_pose = inputs(anchor_state)['pose']
                training_state = states[2]
                training_pose = inputs(training_state)['pose']
                anchor_field = current_field(anchor, energies[1], values, anchor_state[7:13],
                                             time=.10, acceleration_columns=False)
                anchor_eta_rate, _, _ = surface_rates(anchor['points'], anchor['faces'][anchor['free_faces']],
                    anchor['free_nodes'], anchor_field['parent_map'], anchor['parent_psi'], anchor_field['flux'])
                shape1 = np.r_[training_pose.position-anchor_pose.position,
                    (Rotation.from_quat(training_pose.orientation, scalar_first=True)
                     *Rotation.from_quat(anchor_pose.orientation, scalar_first=True).inv()).as_rotvec(),
                    training_state[13:13+size]-anchor_state[13:13+size]]
                shape2 = .05*np.r_[anchor_state[7:13], anchor_eta_rate]

                def shape_geometry(change):
                    orientation = (Rotation.from_rotvec(change[3:6])
                        *Rotation.from_quat(anchor_pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                    arrays = refine(evaluate_template(chart,
                        pose=InitialState(anchor_pose.position+change[:3], orientation),
                        parent_eta=anchor_state[13:13+size]+change[6:],
                        parent_psi=anchor_state[13+size:])[0])
                    same_basis(arrays)
                    return arrays

                epsilon = training['epsilon']

                def geometry_direction(shape):
                    plus, minus = shape_geometry(epsilon*shape), shape_geometry(-epsilon*shape)
                    return (plus['points']-minus['points'])/(2*epsilon), (plus, minus)

                g1, _ = geometry_direction(shape1)
                g2, _ = geometry_direction(shape2)
                coefficient = inner(g1, g2)/inner(g1, g1)
                perpendicular = g2-coefficient*g1
                raw_rms = float(np.sqrt(inner(perpendicular, perpendicular)))
                scale1 = float(np.sqrt(inner(g1, g1)))
                training.update(raw_perpendicular_geometry_rms_m=raw_rms,
                    raw_perpendicular_to_first_ratio=raw_rms/scale1,
                    subtracted_shape_coefficient=coefficient)
                if raw_rms <= 128*np.finfo(float).eps*(np.sqrt(inner(g2, g2))+abs(coefficient)*scale1):
                    raise ValueError('physical training directions are numerically dependent')
                amplification = scale1/raw_rms
                shape2 = amplification*(shape2-coefficient*shape1)
                training['perpendicular_parameter_amplification'] = amplification
                directions, derivatives, training_points, training_energies = [], [], [], []
                for index, shape in enumerate((shape1, shape2)):
                    geometric, pair = geometry_direction(shape)
                    directions.append(geometric)
                    matrices = []
                    for sign, arrays in zip((1, -1), pair):
                        from .cuda_operator import prepare_energy
                        matrices.append(timed(f'tangent{index}-sign{sign}-energy', lambda: prepare_energy(
                            arrays['faces'], TrianglePanels(arrays['points'][arrays['faces']], 1), order=12).energy))
                        training['assembly_count'] += 1
                        training_points.append(arrays['points'])
                        training_energies.append(matrices[-1])
                        saved[f'tangent{index}_sign{sign}_energy'] = matrices[-1]
                    derivatives.append((matrices[0]-matrices[1])/(2*epsilon))
                    saved[f'tangent{index}_geometry'] = geometric
                    saved[f'tangent{index}_energy'] = derivatives[-1]
                    saved[f'tangent{index}_physical_shape'] = shape
                    training.setdefault('physical_perturbations', []).append(dict(
                        maximum_position_step_m=float(epsilon*np.max(abs(shape[:3]))),
                        rotation_step_rad=float(epsilon*np.linalg.norm(shape[3:6])),
                        maximum_parent_height_step_m=float(epsilon*np.max(abs(shape[6:]))),
                        plus_geometry_rms_step_m=float(np.sqrt(inner(pair[0]['points']-anchor['points'],
                                                                    pair[0]['points']-anchor['points'])))))
                gram = np.array([[inner(a, b) for b in directions] for a in directions])
                np.linalg.cholesky(gram)
                condition = float(np.linalg.cond(gram))
                if condition*np.finfo(float).eps > 1e-8:
                    raise ValueError('two-direction metric exceeds the numerical conditioning budget')

                # Fit in the SAME geometry coordinates used at evaluation.
                # Physical ±epsilon paths need not project to exactly ±epsilon
                # on an axis when the chart map itself has curvature.
                feature_rows = []
                training_orthogonal = []
                for points in training_points:
                    alpha, residual = _shape_coordinates(points, anchor['points'], weight, directions, gram)
                    feature_rows.append(_quadratic_features(alpha, epsilon))
                    training_orthogonal.append(residual)
                feature_matrix = np.array(feature_rows)
                feature_condition = float(np.linalg.cond(feature_matrix))
                if feature_condition*np.finfo(float).eps > 1e-8:
                    raise ValueError('quadratic fit exceeds the numerical conditioning budget')
                quadratic_coefficients = np.linalg.solve(feature_matrix,
                    np.stack([e-energies[1] for e in training_energies]).reshape(4, -1)).reshape(4, *energies[1].shape)

                tangent_prediction, quadratic_prediction = _local_models(anchor['points'], energies[1],
                    weight, directions, gram, derivatives, quadratic_coefficients, epsilon)
                saved['quadratic_energy_coefficients'] = quadratic_coefficients
                saved['quadratic_feature_matrix'] = feature_matrix
                saved['local_anchor_energy'] = energies[1]
                training['quadratic_feature_condition'] = feature_condition
                training['training_orthogonal_geometry_rms_m'] = training_orthogonal
                training['quadratic_scope'] = 'Separable quadratic in the same geometry coordinates; mixed alpha1*alpha2 term is unmeasured and omitted.'

                rate_coefficients = np.linalg.solve(gram, [inner(g, geometry_rate) for g in directions])
                rate_residual = geometry_rate-sum(a*g for a, g in zip(rate_coefficients, directions))
                training.update(complete=True, gram=gram.tolist(), gram_condition=condition,
                    pressure_direction_out_of_subspace_fraction=float(np.sqrt(
                        inner(rate_residual, rate_residual)/inner(geometry_rate, geometry_rate))))
                saved['tangent_anchor_points'] = anchor['points']
            except Exception as error:
                tangent_prediction = quadratic_prediction = None
                training['failure'] = dict(type=type(error).__name__, message=str(error))
            finally:
                training['seconds'] = perf_counter()-training_started
        if tangent_from is not None:
            reuse_started = perf_counter()
            tangent_from = Path(tangent_from)
            trained = json.loads(tangent_from.read_text())
            if (not trained.get('complete') or not trained.get('two_direction_training', {}).get('complete')
                    or trained['input_sha256'] != report['input_sha256'] or trained['chart'] != report['chart']):
                raise ValueError('local model requires a complete matching training audit/chart')
            for name, expected in trained['source_sha256'].items():
                if name != Path(__file__).name and report['source_sha256'][name] != expected:
                    raise ValueError(f'local model producer helper changed: {name}')
            with np.load(tangent_from.with_suffix('.npz'), allow_pickle=False) as data:
                training_data = dict(data)
            if (not np.array_equal(training_data['selected_states'][:3], states[:3])
                    or not np.allclose(training_data['tangent_anchor_points'], base[1]['points'], rtol=0, atol=1e-12)
                    or not np.array_equal(training_data['geometry_weight'], weight)):
                raise ValueError('cached local model has different anchor states or metric')
            directions = [training_data[f'tangent{i}_geometry'] for i in range(2)]
            derivatives = [training_data[f'tangent{i}_energy'] for i in range(2)]
            gram = np.array(trained['two_direction_training']['gram'])
            epsilon = trained['two_direction_training']['epsilon']
            anchor_energy = training_data.get('local_anchor_energy', training_data['anchor1_energy'])
            tangent_prediction, quadratic_prediction = _local_models(base[1]['points'], anchor_energy,
                weight, directions, gram, derivatives, training_data['quadratic_energy_coefficients'], epsilon)
            training = dict(trained['two_direction_training'], reused_unchanged=True, new_training_assemblies=0,
                source_path=str(tangent_from), source_sha256=sha256(tangent_from.read_bytes()).hexdigest(),
                source_arrays_sha256=sha256(tangent_from.with_suffix('.npz').read_bytes()).hexdigest())
            if 'seconds' in training:
                training['original_axis_setup_seconds'] = training.pop('seconds')
            training['timing_scope'] = 'Historical axis setup and this-run work are separate, not total six-sample training time.'
            report['two_direction_training'] = training
            for key, value in training_data.items():
                if key.startswith(('tangent', 'quadratic')) or key == 'local_anchor_energy':
                    saved[key] = value
            if mixed_curvature:
                if len(training_data['quadratic_energy_coefficients']) != 4:
                    raise ValueError('mixed extension requires the original separable model')
                anchor_pose = inputs(states[1])['pose']
                shapes = [training_data[f'tangent{i}_physical_shape'] for i in range(2)]
                feature_rows, fit_energies, orthogonal = [], [], []
                for label, shape in [('axis0', shapes[0]), ('axis1', shapes[1]), ('mixed', shapes[0]+shapes[1])]:
                    for sign in (1, -1):
                        change = sign*epsilon*shape
                        orientation = (Rotation.from_rotvec(change[3:6])
                            *Rotation.from_quat(anchor_pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                        arrays = refine(evaluate_template(chart,
                            pose=InitialState(anchor_pose.position+change[:3], orientation),
                            parent_eta=states[1][13:13+size]+change[6:],
                            parent_psi=states[1][13+size:])[0])
                        same_basis(arrays)
                        alpha, residual = _shape_coordinates(arrays['points'], base[1]['points'],
                                                              weight, directions, gram)
                        feature_rows.append(_quadratic_features(alpha, epsilon, mixed=True))
                        orthogonal.append(residual)
                        if label == 'mixed':
                            from .cuda_operator import prepare_energy
                            matrix = timed(f'mixed-sign{sign}-energy', lambda: prepare_energy(arrays['faces'],
                                TrianglePanels(arrays['points'][arrays['faces']], 1), order=12).energy)
                            saved[f'mixed_sign{sign}_energy'] = matrix
                            training['new_training_assemblies'] += 1
                        else:
                            index = int(label[-1])
                            matrix = training_data[f'tangent{index}_sign{sign}_energy']
                        fit_energies.append(matrix)
                feature_matrix = np.array(feature_rows)
                rhs = np.stack([e-anchor_energy for e in fit_energies]).reshape(6, -1)
                # Same six physical samples for both fits: adding an interaction
                # feature must be separated from merely adding oblique samples.
                for label, count in (('separable_refit', 4), ('full', 5)):
                    fit, _, rank, singular = np.linalg.lstsq(feature_matrix[:, :count], rhs, rcond=None)
                    condition = float(singular[0]/singular[-1])
                    if rank != count or condition*np.finfo(float).eps > 1e-8 or not np.isfinite(fit).all():
                        raise ValueError('mixed-curvature fit is rank deficient or numerically uncertain')
                    fit = fit.reshape(count, *anchor_energy.shape)
                    training[label+'_feature_condition'] = condition
                    saved['quadratic_'+label+'_coefficients'] = fit
                saved['quadratic_original_coefficients'] = training_data['quadratic_energy_coefficients']
                saved['quadratic_energy_coefficients'] = saved['quadratic_full_coefficients']
                saved['quadratic_feature_matrix'] = feature_matrix
                _, quadratic_prediction = _local_models(base[1]['points'], anchor_energy, weight, directions,
                    gram, derivatives, saved['quadratic_energy_coefficients'], epsilon)
                training.update(reused_unchanged=False, assembly_count=training['assembly_count']+2,
                    quadratic_scope='Full two-direction quadratic from six physical perturbations; no later target data used.',
                    training_orthogonal_geometry_rms_m=orthogonal)
            for label in ('original', 'separable_refit'):
                key = 'quadratic_'+label+'_coefficients'
                if key in saved:
                    _, predictor = _local_models(base[1]['points'], anchor_energy, weight, directions,
                        gram, derivatives, saved[key], epsilon)
                    extra_local_models['two_direction_quadratic_'+label] = predictor
            rate_alpha = np.linalg.solve(gram, [inner(g, geometry_rate) for g in directions])
            rate_residual = geometry_rate-sum(a*g for a, g in zip(rate_alpha, directions))
            training['pressure_direction_out_of_subspace_fraction'] = float(np.sqrt(
                inner(rate_residual, rate_residual)/inner(geometry_rate, geometry_rate)))
            training['seconds_this_run'] = perf_counter()-reuse_started
        past_models = {}
        if past_enabled:
            from .boat_past_operator import fit, predict, direction_diagnostics
            training_started = perf_counter()
            past_data, past_source = None, None
            if past_from is not None:
                past_from = Path(past_from)
                trained = json.loads(past_from.read_text())
                if (not trained.get('complete') or not trained.get('past_training', {}).get('complete')
                        or trained['input_sha256'] != report['input_sha256'] or trained['chart'] != report['chart']):
                    raise ValueError('past model requires a complete matching training audit/chart')
                for name, expected in trained['source_sha256'].items():
                    if name != Path(__file__).name and report['source_sha256'][name] != expected:
                        raise ValueError(f'past model helper changed: {name}')
                with np.load(past_from.with_suffix('.npz'), allow_pickle=False) as data:
                    past_source = dict(data)
                if (not np.array_equal(past_source['selected_states'][:2], states[:2])
                        or any(not np.allclose(past_source[key], arrays['points'], rtol=0, atol=1e-12)
                            for key, arrays in zip(('anchor0_points', 'anchor1_points'), base[:2]))):
                    raise ValueError('past model anchor states/points differ')
                past_data = {key: value for key, value in past_source.items() if key.startswith('past_')}
                report['past_training'] = dict(trained['past_training'], reused_unchanged=True,
                    new_training_assemblies=0, new_anchor_check_assemblies=0,
                    source_path=str(past_from), source_sha256=sha256(past_from.read_bytes()).hexdigest(),
                    source_arrays_sha256=sha256(past_from.with_suffix('.npz').read_bytes()).hexdigest())

            records, record_energies, current_offsets, current_offset_energies = [], [], {}, {}
            new_training, new_checks = 0, 0
            for index, (arrays, state, time) in enumerate(zip(base[:2], states[:2], (.05, .10))):
                center_field = current_field(arrays, energies[index], values, state[7:13],
                                             acceleration_columns=False, time=time)
                er, pr, _ = surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
                    arrays['free_nodes'], center_field['parent_map'], arrays['parent_psi'], center_field['flux'])
                records.append(arrays); record_energies.append(energies[index])
                pose_at_anchor = inputs(state)['pose']
                for h in ((.001,) if index == 0 else (.001, .0005)):
                    for sign in (1, -1):
                        dt = sign*h
                        orientation = (Rotation.from_rotvec(dt*state[10:13])
                            *Rotation.from_quat(pose_at_anchor.orientation, scalar_first=True)).as_quat(scalar_first=True)
                        moved = refine(evaluate_template(chart,
                            pose=InitialState(pose_at_anchor.position+dt*state[7:10], orientation),
                            parent_eta=state[13:13+size]+dt*er, parent_psi=state[13+size:]+dt*pr)[0])
                        same_basis(moved)
                        key = f'past_anchor{index}_offset_{dt}'
                        if past_source is None:
                            from .cuda_operator import prepare_energy
                            energy = timed(key+'-energy', lambda: prepare_energy(moved['faces'],
                                TrianglePanels(moved['points'][moved['faces']], 1), order=12).energy)
                            if h == .001:
                                new_training += 1
                            else:
                                new_checks += 1
                        else:
                            if not np.allclose(past_source[key+'_points'], moved['points'], rtol=0, atol=1e-12):
                                raise ValueError('past pressure-offset geometry differs from training audit')
                            energy = past_source[key+'_energy']
                        saved[key+'_points'], saved[key+'_energy'] = moved['points'], energy
                        if h == .001:
                            records.append(moved); record_energies.append(energy)
                        if index == 1:
                            current_offsets[dt], current_offset_energies[dt] = moved, energy
            if past_data is None:
                past_data, description = fit(records, record_energies)
                report['past_training'] = dict(description, new_training_assemblies=new_training,
                    new_anchor_check_assemblies=new_checks, reused_unchanged=False,
                    historical_records=6, historical_training_assemblies=6,
                    additional_half_step_reference_assemblies=2)
            saved.update(past_data)
            past_models = dict(past_constrained=lambda arrays: predict(past_data, arrays),
                               past_ordinary=lambda arrays: predict(past_data, arrays, ordinary=True))
            report['past_anchor_check'] = []
            for name in ('fresh', *past_models):
                energy = energies[1] if name == 'fresh' else past_models[name](base[1])[0]
                moved_energies = (current_offset_energies if name == 'fresh' else
                    {dt: past_models[name](arrays)[0] for dt, arrays in current_offsets.items()})
                _, _, _, pressures = evaluate_stencil(body, base[1], states[1], values, .10,
                    energy, current_offsets, moved_energies)
                report['past_anchor_check'].append(dict(name=name, pressure=pressures))
            report['past_training']['pressure_direction'] = direction_diagnostics(past_data,
                offsets[.0005], offsets[-.0005])
            report['past_training']['seconds_this_run'] = perf_counter()-training_started
        def surface_rms(delta_parent):
            faces = target['faces'][target['free_faces']]
            triangles = target['points'][faces]
            area = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])[:, 2]/2
            local = np.full(len(target['points']), -1)
            local[target['free_nodes']] = np.arange(len(target['free_nodes']))
            nodal = (actual['parent_map']@delta_parent)[local[faces]]
            square = (np.sum(nodal*nodal, axis=1)+np.sum(nodal, axis=1)**2)/12
            return float(np.sqrt(area@square/area.sum()))

        names = ['fresh', 'frozen', 'secant']
        if previous is not None:
            names += ['exact_center_secant_slope', 'secant_center_exact_slope', 'exact_center_exact_slope']
        if tangent_prediction is not None:
            local_models = dict(two_direction=tangent_prediction, two_direction_quadratic=quadratic_prediction,
                                **extra_local_models)
            names += list(local_models)
        else:
            local_models = {}
        local_models.update(past_models)
        names += list(past_models)
        predicted_target, predicted_description = prediction(target)
        for name in names:
            row = dict(name=name, complete=False, pressure=[])
            report['models'].append(row)
            before = perf_counter()
            try:
                if name in ('fresh', 'exact_center_secant_slope', 'exact_center_exact_slope'):
                    energy, description = target_energy, {}
                elif name == 'frozen':
                    energy, description = energies[1], {}
                elif name in local_models:
                    energy, description = local_models[name](target)
                else:
                    energy, description = predicted_target, predicted_description
                model_offsets = {}
                for h in (.001, .0005):
                    for dt in (h, -h):
                        arrays = offsets[dt]
                        if name == 'fresh':
                            offset_energy, offset_description = offset_energies[dt], {}
                        elif name == 'frozen':
                            offset_energy, offset_description = energies[1], {}
                        elif name in local_models:
                            offset_energy, offset_description = local_models[name](arrays)
                        elif name == 'exact_center_secant_slope':
                            predicted_offset, offset_description = prediction(arrays)
                            offset_energy = target_energy+predicted_offset-predicted_target
                        elif name in ('secant_center_exact_slope', 'exact_center_exact_slope'):
                            exact_slope = (offset_energies[h]-offset_energies[-h])/(2*h)
                            offset_energy, offset_description = energy+dt*exact_slope, {}
                        else:
                            offset_energy, offset_description = prediction(arrays)
                        model_offsets[dt] = offset_energy
                        if offset_description:
                            row.setdefault('offset_predictions', {})[str(dt)] = offset_description
                field, eta, psi, pressures = evaluate_stencil(body, target, y, values,
                    target_time, energy, offsets, model_offsets)
                defect = target_energy@field['potential']-field['right']
                unknown = np.setdiff1d(np.arange(len(energy)), target['free_nodes'])
                row.update(geometry_prediction=description, pressure=pressures,
                    true_current_interior_relative_defect=float(np.linalg.norm(defect[unknown])
                        /max(np.linalg.norm(field['right'][unknown]), 1e-30)),
                    parent_height_rate_error_surface_rms_m_per_s=surface_rms(eta-eta_rate),
                    parent_potential_rate_error_surface_rms_m2_per_s2=surface_rms(psi-psi_rate))
                row.update(complete=True, seconds=perf_counter()-before)
                saved[name+'_potential'] = field['potential']
            except Exception as error:
                row['failure'] = dict(type=type(error).__name__, message=str(error))
                if name == 'fresh':
                    raise
            print(json.dumps(row), flush=True)
        saved.update(anchor0_energy=energies[0], anchor1_energy=energies[1], target_energy=target_energy,
            anchor0_points=origin, anchor1_points=base[1]['points'], target_points=target['points'],
            selected_states=states, geometry_weight=weight)
        for dt, energy in offset_energies.items():
            saved[f'offset_{dt}_energy'] = energy
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['seconds'] = perf_counter()-started
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--operator-from', help='Complete saved audit for offline value/slope decomposition')
    parser.add_argument('--two-directions', action='store_true')
    parser.add_argument('--target-time', type=float, choices=(.15, .20), default=.15)
    parser.add_argument('--tangent-from', help='Reuse one completed local model without retraining')
    parser.add_argument('--mixed-curvature', action='store_true')
    parser.add_argument('--past-directions', action='store_true')
    parser.add_argument('--past-from', help='Reuse one past-only model without retraining')
    run(**vars(parser.parse_args()))
