"""Unqualified full nonlinear P4 canoe release from the actual saved stage.

Every RHS prepares a fresh current fluid operator and evaluates the complete
canonical action gradient. Only the reference geometry extension is reused.
The rim check is a conservative GLOBAL height separation bound, not an exact
P4 rim intersection/flooding algorithm. No numerical job runs on import.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial, prod
from pathlib import Path
import subprocess
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from .canoe_contact import ContactTopologyError
from .p4_canoe_field import prepare
from .p4_variational_canoe import evaluate
from .p4_whole_geometry import P4WholeCanoeGeometry
from .p4_canoe_stage import verify_master_identity
from .p4_tetra import indices, face_nodes
from .geometry_rk23 import GeometryRK23


@lru_cache(None)
def _bernstein_conversion():
    """Fixed exact rational degree4 triangle conversion, enclosed in float64."""
    alpha = indices(3)
    table = []
    for node in alpha:
        bary = [Fraction(int(value), 4) for value in node]
        row = [Fraction(24, prod(factorial(int(value)) for value in power))
               *prod(bary[j]**int(power[j]) for j in range(3)) for power in alpha]
        table.append(row)
    augmented = [row+[Fraction(int(i == j)) for j in range(15)] for i, row in enumerate(table)]
    for i in range(15):
        pivot = next(j for j in range(i, 15) if augmented[j][i])
        augmented[i], augmented[pivot] = augmented[pivot], augmented[i]
        value = augmented[i][i]
        augmented[i] = [entry/value for entry in augmented[i]]
        for j in range(15):
            if j != i and augmented[j][i]:
                value = augmented[j][i]
                augmented[j] = [a-value*b for a, b in zip(augmented[j], augmented[i])]
    lower = np.empty((15, 15)); upper = np.empty((15, 15))
    for i in range(15):
        for j in range(15):
            exact = augmented[i][15+j]
            value = float(exact)
            rounded = Fraction.from_float(value)
            lower[i, j] = np.nextafter(value, -np.inf) if rounded > exact else value
            upper[i, j] = np.nextafter(value, np.inf) if rounded < exact else value
    return lower, upper


def _rim_bound(owner, points, position, orientation, surface_nodes, outer):
    """min numerical rigid rim z minus an outward-rounded Bernstein upper bound.

    A straight rim segment's minimum z is at an endpoint. The P4 Bernstein
    convex-hull bound applies to every parameter point, including curves inside
    flat hull facets. No value is interpolated/extrapolated into the hull hole.
    """
    lower, upper = _bernstein_conversion()
    height = points[surface_nodes, 2]
    coefficients = np.zeros_like(height)
    for j in range(15):
        value = height[:, j, None]
        multiplier = np.where(value >= 0, upper[:, j], lower[:, j])
        product_upper = np.nextafter(value*multiplier, np.inf)
        coefficients = np.nextafter(coefficients+product_upper, np.inf)
    maximum = float(np.max(coefficients))
    chart = owner.chart.linear
    rim = (chart.rim_vertices-chart.cg)@rotation_matrix(orientation).T+position
    if not len(rim) or not np.all(np.isfinite(rim)) or not np.isfinite(maximum):
        raise ContactTopologyError('P4 global rim separation cannot be certified')
    # Rigid coordinates use the same numerical transform as the actual hull;
    # reserve a stated roundoff margin rather than certify exact real rotations.
    margin = float(64*np.finfo(float).eps*max(1., np.max(abs(rim)), abs(maximum)))
    if np.any(rim[:, :2] <= outer[0]+margin) or np.any(rim[:, :2] >= outer[1]-margin):
        raise ContactTopologyError('physical rim XY is outside the certified outer tank box')
    minimum = float(np.min(rim[:, 2]))
    clearance = float(np.nextafter(minimum-maximum-margin, -np.inf))
    if clearance <= 0:
        raise ContactTopologyError('global P4 height bound cannot certify positive rim separation; '
                                   'this is not a flooding classification')
    return clearance, minimum, maximum, margin


def _revision():
    try:
        return subprocess.check_output(['git', 'rev-parse', 'HEAD'],
            cwd=Path(__file__).resolve().parents[3], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def run(stage_path, *, duration=.001, max_step=.001, rtol=1e-4,
        atol=1e-10, derivative_step=1e-4):
    """Integrate the saved all-six-velocity/nonflat actual-stage initial state."""
    start = perf_counter()
    sources = ['p4_canoe_release.py', 'p4_canoe_stage.py', 'geometry_rk23.py',
        'p4_canoe_field.py', 'p4_variational_canoe.py', 'p4_whole_geometry.py',
        'p4_volume_map.py', 'p4_canoe_contact.py', 'canoe_contact.py', 'p4_geometry.py',
        'p4_surface.py', 'p4_boundary_action.py', 'p4_tetra.py', 'canonical_canoe.py',
        'moving_canoe.py', 'pocket.py', 'p4_whole_canoe.py', 'p4_canoe.py']
    source_hashes = {name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                     for name in sources}
    source_hashes['topwater/frames.py'] = sha256(Path(rotation_matrix.__code__.co_filename).read_bytes()).hexdigest()
    source_hashes['scipy/integrate/_ivp/rk.py'] = sha256(Path(
        GeometryRK23.__mro__[1]._step_impl.__code__.co_filename).read_bytes()).hexdigest()
    if (not all(np.isfinite(value) and value > 0 for value in
                (duration, max_step, rtol, atol, derivative_step)) or rtol >= 1):
        raise ValueError('finite positive run controls and rtol<1 required')
    stage_path = Path(stage_path)
    metadata = json.loads(stage_path.with_suffix('.json').read_text())
    with np.load(stage_path) as saved:
        canonical = np.array(saved['canonical'], dtype=float)
        physical = np.array(saved['physical'], dtype=float)
        saved_rate = np.array(saved['canonical_rate'], dtype=float)
        saved_q = np.array(saved['q'], dtype=float)
    if (canonical.ndim != 1 or len(canonical) % 2 or saved_rate.shape != canonical.shape
            or saved_q.shape != (len(canonical)//2,)
            or not np.array_equal(saved_q, canonical[:len(saved_q)])
            or physical.shape != (13+2*(len(saved_q)-6),)
            or not all(np.all(np.isfinite(value)) for value in (canonical, physical, saved_rate))
            or metadata['case'] != 'actual nonlinear P4 canonical-canoe stage'):
        raise ValueError('finite matching actual-stage NPZ and JSON initial data required')
    owner = P4WholeCanoeGeometry.make()
    master_identity = verify_master_identity(owner, stage_path)
    size = 6+len(owner.free_nodes)
    if canonical.shape != (2*size,) or metadata['eta_modes'] != len(owner.free_nodes):
        raise ValueError('saved stage differs from rebuilt P4 surface master space')
    surface_nodes = face_nodes(owner.mesh, owner.free_faces)
    outer = (owner.mesh.points[:, :2].min(axis=0), owner.mesh.points[:, :2].max(axis=0))
    calls = 0
    initial = {}
    cache = {}

    def rhs(time, value):
        nonlocal calls
        calls += 1
        print(json.dumps(dict(progress=True, stage='fresh nonlinear RHS started',
                              simulated_time=time, rhs_evaluations=calls)), flush=True)
        data = prepare(owner, value[:size])
        clearance, rim_z, surface_z, margin = _rim_bound(owner, data['state'].points,
            value[:3], data['orientation'], surface_nodes, outer)
        current = evaluate(owner, data, value[size:], derivative_step=derivative_step)
        if not initial:
            error = float(np.max(abs(current['physical_state']-physical)))
            if error > 1e-9:
                raise ValueError(f'fresh initial physical state differs from saved stage: {error}')
            initial.update(physical_base_maximum_absolute_error=error,
                rhs_maximum_absolute_difference_from_saved_stage=float(np.max(abs(current['rate']-saved_rate))))
        cache.clear()
        cache.update(time=float(time), canonical=value.copy(), physical=current['physical_state'],
            hamiltonian=float(current['hamiltonian']), volume=float(current['fluid_volume']),
            rim_bound=clearance, rim_z=rim_z, surface_z=surface_z, margin=margin,
            weak_normal_relative=current['full_weak_normal_velocity_relative_residual'],
            fluid_work_residual=current['fluid_work_identity_absolute_residual'])
        print(json.dumps(dict(progress=True, stage='fresh nonlinear RHS complete',
            simulated_time=time, rhs_evaluations=calls, hamiltonian=cache['hamiltonian'],
            rim_clearance_lower_bound=clearance)), flush=True)
        return current['rate']

    # Explicit first_step prevents the extra heuristic trial RHS for a tiny
    # smoke interval. Adaptive RK23 still rejects it if accuracy requires.
    solver = GeometryRK23(rhs, 0., canonical, duration, rtol=rtol, atol=atol,
                  max_step=max_step, first_step=min(duration, max_step))
    times = []; states = []; canonical_states = []; energies = []; volumes = []
    bounds = []; rim_heights = []; surface_heights = []; margins = []
    weak_errors = []; work_errors = []

    def record():
        if cache['time'] != solver.t or not np.array_equal(cache['canonical'], solver.y):
            raise RuntimeError('accepted endpoint lacks exact matching current-state diagnostics')
        times.append(solver.t); states.append(cache['physical'].copy())
        canonical_states.append(solver.y.copy()); energies.append(cache['hamiltonian'])
        volumes.append(cache['volume']); bounds.append(cache['rim_bound'])
        rim_heights.append(cache['rim_z']); surface_heights.append(cache['surface_z'])
        margins.append(cache['margin']); weak_errors.append(cache['weak_normal_relative'])
        work_errors.append(cache['fluid_work_residual'])

    record()
    while solver.status == 'running':
        message = solver.step()
        if solver.status == 'failed':
            raise RuntimeError(message)
        record()
    energies, volumes = np.asarray(energies), np.asarray(volumes)
    physical_states = np.asarray(states)
    energy_change = float(np.max(abs(energies-energies[0])))
    report = dict(qualified=False, case='fully nonlinear P4 canonical canoe release',
        source_stage=str(stage_path), source_stage_metadata=metadata,
        source_npz_sha256=sha256(stage_path.read_bytes()).hexdigest(), source_revision=_revision(),
        runner_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        source_sha256=source_hashes,
        geometry_rk23_sha256=sha256(Path(GeometryRK23._step_impl.__code__.co_filename).read_bytes()).hexdigest(),
        master_identity_assumption=master_identity,
        initial_source_checks=initial, saved_stage_rate_used_for_integration=False,
        reference_points_sha256=sha256(np.asarray(owner.mesh.points, dtype='<f8').tobytes()).hexdigest(),
        elements_sha256=sha256(np.asarray(owner.mesh.elements, dtype='<i8').tobytes()).hexdigest(),
        free_nodes_sha256=sha256(np.asarray(owner.free_nodes, dtype='<i8').tobytes()).hexdigest(),
        reference_linear_node_count=len(owner.domain.reference_points),
        eta_modes=len(owner.free_nodes), geometry_preflight=owner.preflight,
        duration=duration, simulated_seconds=solver.t, rtol=rtol, atol=atol,
        max_step=max_step, first_step=min(duration, max_step), derivative_step=derivative_step,
        rhs_evaluations=calls, accepted_steps=len(times)-1,
        invalid_geometry_trial_rejections=solver.invalid_geometry_rejections,
        initial_excitation_energy_joules=float(energies[0]), maximum_hamiltonian_change_joules=energy_change,
        maximum_hamiltonian_change_over_excitation=energy_change/abs(float(energies[0])) if energies[0] else None,
        maximum_fluid_volume_change_m3=float(np.max(abs(volumes-volumes[0]))),
        maximum_position_change_m=np.max(abs(physical_states[:, :3]-physical_states[0, :3]), axis=0).tolist(),
        maximum_world_velocity_change=np.max(abs(physical_states[:, 7:13]-physical_states[0, 7:13]), axis=0).tolist(),
        world_velocity_units=['m/s']*3+['rad/s']*3,
        minimum_rim_clearance_lower_bound_m=float(min(bounds)),
        rim_bound_method='Minimum transformed straight physical rim vertex z minus the maximum '
            'degree4 Bernstein coefficient over EVERY actual P4 FS patch. Fixed rational conversion '
            'is enclosed in float64 and products/sums round outward. Numerical rigid coordinates '
            'use the hull convention with a recorded64eps geometric margin; rim XY must remain '
            'strictly inside the fixed outer tank box. Nonpositive bound stops as uncertified, '
            'not FloodingError. No P1 interpolation or hull-hole filling.',
        scope='Full current nonlinear six-DOF canonical evolution with fresh current pressure/energy '
            'operator and all physical action terms each RHS. Fixed chart, sampled geometry guards '
            'and conservative positive rim separation only. General P4 flooding/loss-contact event '
            'coverage, remeshing, long-time stability, temporal/spatial/Capytaine and performance '
            'qualification remain unestablished.')
    arrays = dict(times=np.asarray(times), states=physical_states,
        canonical_states=np.asarray(canonical_states), hamiltonian=energies, fluid_volume=volumes,
        rim_clearance_lower_bound=np.asarray(bounds), minimum_rim_z=np.asarray(rim_heights),
        free_surface_height_upper_bound=np.asarray(surface_heights), rim_roundoff_margin=np.asarray(margins),
        weak_normal_relative_residual=np.asarray(weak_errors, dtype=float),
        fluid_work_identity_absolute_residual=np.asarray(work_errors))
    report['diagnostic_wall_seconds'] = perf_counter()-start
    report['diagnostic_wall_scope'] = ('Input loading, one geometry build, current RHS evaluations '
        'and accepted diagnostics; output writing excluded. Concurrency/hardware are uncontrolled. '
        'This is not qualified throughput.')
    return report, arrays


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', type=Path)
    parser.add_argument('--duration', type=float, default=.001)
    parser.add_argument('--max-step', type=float, default=.001)
    parser.add_argument('--rtol', type=float, default=1e-4)
    parser.add_argument('--atol', type=float, default=1e-10)
    parser.add_argument('--derivative-step', type=float, default=1e-4)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-canoe-release.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report, arrays = run(args.stage, duration=args.duration, max_step=args.max_step,
            rtol=args.rtol, atol=args.atol, derivative_step=args.derivative_step)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(args.output.with_suffix('.npz'), **arrays)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
