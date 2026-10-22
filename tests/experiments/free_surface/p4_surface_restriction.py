"""Actual calm-canoe surface restriction falsifier; no new fluid solve.

This initial zero-trace state is a necessary control, not model qualification.
Nonzero waves and trajectory/enrichment tests remain required even if it passes.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
from time import perf_counter

import numpy as np
import scipy
from scipy.linalg import solve
from scipy.sparse.linalg import LinearOperator
from threadpoolctl import threadpool_limits

from .p4_canoe_field import geometry_data
from .p4_canoe_stage import verify_master_identity
from .p4_fluid_ordering import fresh, load_input
from .p4_metric_action import MetricAction
from .p4_surface_basis import rectangle_cosines
from .p4_variational_canoe import action_rate
from .p4_whole_geometry import P4WholeCanoeGeometry


def run(stage, old_matrix, output):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(stage)
    with np.load(stage, allow_pickle=False) as saved:
        q = saved['q'].copy();rate = saved['canonical_rate'].copy()
        physical = saved['physical'].copy()
    _, old, manifest = load_input(old_matrix)
    report = dict(qualified=False, complete=False,
        scope='Initial calm actual canoe, full physical field retained; zero-trace force control is not adequacy proof',
        stage_sha256=sha256(stage.read_bytes()).hexdigest(), old_matrix_npz_sha256=manifest['npz_sha256'],
        runtime=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
        phases={}, cases=[], source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('p4_surface_restriction.py', 'p4_surface_basis.py', 'p4_variational_canoe.py',
                         'p4_metric_action.py', 'p4_canoe_field.py', 'p4_whole_geometry.py')})
    arrays = {}

    def timed(name, function):
        print(json.dumps(dict(phase=name, event='start')), flush=True)
        before = perf_counter()
        try:
            return function()
        finally:
            report['phases'][name] = perf_counter()-before
            print(json.dumps(dict(phase=name, seconds=report['phases'][name])), flush=True)

    try:
        if not np.array_equal(q, old['q']) or np.any(q[6:] != 0) or np.any(old['psi'] != 0):
            raise ValueError('same captured calm zero-height/trace actual canoe required')
        owner = timed('geometry_startup', P4WholeCanoeGeometry.make)
        report['master_identity'] = verify_master_identity(owner, stage)
        size = len(q)
        if (rate.shape != (2*size,) or physical.shape != (13+2*(size-6),)
                or not all(np.all(np.isfinite(v)) for v in (q, rate, physical))):
            raise ValueError('matching finite actual saved field/rate required')
        if (not np.allclose(physical[7:13], old['reference_velocity'], rtol=0, atol=1e-12)
                or not np.array_equal(physical[13+size-6:], old['psi'])):
            raise ValueError('saved physical body velocity/trace differs from the captured field')
        data = timed('current_geometry', lambda: geometry_data(owner, q))
        inside, free = data['inside'], data['free']
        if not np.array_equal(inside, old['inside']) or not np.array_equal(free, old['free']):
            raise ValueError('saved/current field master IDs differ')
        phi = np.zeros(len(owner.mesh.points));phi[inside] = old['reference_potential']
        cache = timed('current_metric_cache', lambda: MetricAction(data['state'].points[owner.mesh.elements],
                                                                  owner.mesh.elements, len(owner.mesh.points)))
        data['energy'] = LinearOperator((len(phi),)*2, dtype=float,
            matvec=lambda value: cache.action(value, chunk_size=256)[0])
        kinetic = timed('current_kinetic', lambda: cache.kinetic(phi))
        kinetic_pair = kinetic['kinetic_energy'], kinetic['kinetic_gradient']
        body_qdot = np.linalg.solve(data['body_map'], physical[7:13])
        flux = (data['energy']@phi-data['load']@physical[7:13])[free]-data['body_transport']@physical[7:13]
        full_eta_rate = data['eta_factor'].solve(flux)
        full_current = dict(potential=phi, physical_state=physical, qdot=np.r_[body_qdot, full_eta_rate])
        baseline = timed('full_reference_action', lambda: action_rate(owner, data, full_current, kinetic=kinetic_pair))
        report['reference_rate_maximum_absolute_difference'] = float(np.max(abs(baseline['rate']-rate)))
        report['reference_weak_normal_relative'] = baseline['full_weak_normal_velocity_relative_residual']
        beam = float(np.ptp(owner.domain.body.mesh.vertices[:, 1]))
        mass = data['mass']
        total_norm = float(full_eta_rate@mass@full_eta_rate)
        arrays.update(full_rate=baseline['rate'], free_nodes=free, reference_xy=owner.mesh.points[free, :2],
                      full_eta_rate=full_eta_rate)
        for enrichment in (1., 1.5, 2.):
            before = perf_counter()
            basis, modes = rectangle_cosines(owner.mesh.points[free, :2], beam, enrichment=enrichment)
            reduced_map = basis.T@(data['eta_map']@basis)
            reduced_rate = solve(reduced_map, basis.T@flux)
            eta_rate = basis@reduced_rate
            residual = flux-data['eta_map']@eta_rate
            current = dict(full_current, qdot=np.r_[body_qdot, eta_rate])
            candidate = action_rate(owner, data, current, kinetic=kinetic_pair)
            error = eta_rate-full_eta_rate
            row = dict(enrichment=enrichment, modes=len(modes), basis_bytes=basis.nbytes,
                reduced_map_condition=float(np.linalg.cond(reduced_map)),
                omitted_surface_rate_mass_relative=float(np.sqrt((error@mass@error)/total_norm)),
                omitted_surface_rate_maximum_absolute=float(np.max(abs(error))),
                projected_kinematic_relative_residual=float(np.linalg.norm(basis.T@residual)/np.linalg.norm(basis.T@flux)),
                full_weak_normal_relative=candidate['full_weak_normal_velocity_relative_residual'],
                fluid_work_residual=candidate['fluid_work_identity_absolute_residual'],
                body_momentum_rate_difference=(candidate['rate'][size:size+6]-baseline['rate'][size:size+6]).tolist(),
                body_momentum_rate_reference=baseline['rate'][size:size+6].tolist(),
                wall_seconds=perf_counter()-before)
            report['cases'].append(row)
            arrays[f'eta_rate_{len(modes)}'] = eta_rate
            arrays[f'momentum_rate_{len(modes)}'] = np.r_[candidate['rate'][size:size+6], basis.T@candidate['rate'][size+6:]]
            print(json.dumps(row), flush=True)
        report.update(complete=True, beam=beam, reference_wavelength_minimum=2*beam)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    from . import cuda_operator
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', required=True)
    parser.add_argument('--old-matrix', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with threadpool_limits(1):
        print(json.dumps(run(args.stage, args.old_matrix, args.output), indent=2))
