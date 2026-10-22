"""Bounded actual-boat probe: spend outer quadrature on slender target rows.

This changes numerical integration only. It does not truncate interactions or
regularize matrices. Production kernels and moving-case defaults stay intact.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from topwater import HullMesh
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.reference_boat import WATER_DENSITY, equilibrium
from .boat_world_field import refine_boundary
from .canoe import body_load
from .p4_fluid_ordering import fresh


def run(output, geometry, high_order=6, pressure_check=False, moving_check=False):
    from . import cuda_operator as gpu
    if high_order not in (6, 12):
        raise ValueError('selected target rows require order 6 or 12')
    output, geometry = fresh(output), Path(geometry)
    output.parent.mkdir(parents=True, exist_ok=True)
    arrays = dict(np.load(geometry, allow_pickle=False))
    for _ in range(2):
        arrays = refine_boundary(arrays, region='body')
    faces, points = arrays['faces'], arrays['points']
    panels = TrianglePanels(points[faces], 1)
    load = body_load(HullMesh(points, faces), panels, arrays['body_faces'], equilibrium()[0].position)
    unknown = np.setdiff1d(np.arange(len(points)), arrays['free_nodes'])
    report = dict(complete=False, qualified=False, faces=len(faces), high_order=high_order, records=[],
        geometry_path=str(geometry), geometry_sha256=sha256(geometry.read_bytes()).hexdigest(),
        scope='Integration allocation on one canonical boat boundary; no spatial or trajectory qualification',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_adaptive_quadrature.py', 'cuda_operator.py', 'cuda_energy.cu',
                         'cuda_pair.cuh', 'boat_world_field.py', 'boat_moving.py', 'canoe.py')})
    saved = {}
    # Same analytic pair integral; select one of two independent outer rules.
    # ALL source panels in a target row use the same rule. Pairwise proximity
    # selection would break the pointwise closed-surface H*1 cancellation.
    source = Path(__file__).with_name('cuda_pair.cuh').read_text()+r'''
extern "C" __global__ void adaptive_pairs(
    const double* tri, const double* normal, const double* outward,
    const double* tangent, const double* length, const double* area,
    const double* center, const double* sign, const long long* faces,
    const double* gradient, const double* points, const double* weights,
    int nf, int nv, int nq, double* v, double* h,
    const double* high_points, const double* high_weights, int high_nq,
    const unsigned char* high) {
    int pair=blockDim.x*blockIdx.x+threadIdx.x;
    if (pair>=nf*nf) return;
    int i=pair/nf, j=pair%nf;
    double vsum,linear[3];
    energy_pair<double,const double*>(tri,normal,outward,tangent,length,area,area,
        center,sign,gradient,high[i]?high_points:points,
        high[i]?high_weights:weights,i,j,high[i]?high_nq:nq,vsum,linear);
    v[pair]=vsum;
    for (int c=0;c<3;++c) atomicAdd(h+i*nv+faces[j*3+c],linear[c]);
}
'''
    kernel = gpu.cp.RawKernel(source, 'adaptive_pairs', options=('--std=c++11',))
    high_points, high_weights = triangle_samples(panels.triangles, high_order)
    high_args = (gpu.cp.asarray(np.ascontiguousarray(high_points)),
                 gpu.cp.asarray(high_weights), np.int32(len(high_weights)))
    aspect = panels.edge_lengths.max(axis=1)**2/(2*panels.areas)
    original = gpu._pairs
    reference = None
    try:
        for name, order, cutoff in [('uniform12', 12, None), ('slender10', 3, 10.),
                                    ('slender50', 3, 50.), ('slender100', 3, 100.),
                                    ('uniform6', 6, None), ('uniform3', 3, None)]:
            row = dict(name=name, order=order, target_aspect_cutoff=cutoff, complete=False)
            report['records'].append(row)
            if cutoff is not None:
                high = np.ascontiguousarray(aspect > cutoff, dtype=np.uint8)
                mask = gpu.cp.asarray(high)
                row.update(high_order_pair_fraction=float(high.mean()),
                           sample_fraction_of_uniform12=float((high.mean()*high_order**2
                                                               +(1-high.mean())*9)/144))
                def dispatch(grid, block, args):
                    kernel(grid, block, (*args, *high_args, mask))
                gpu._pairs = dispatch
            else:
                gpu._pairs = original
            # Warm compilation, then measure three matching full assemblies.
            try:
                assembly = gpu.prepare_energy(faces, panels, order=order)
                seconds = []
                for _ in range(3):
                    start = perf_counter()
                    assembly = gpu.prepare_energy(faces, panels, order=order)
                    seconds.append(perf_counter()-start)
                energy = assembly.energy
                response = cho_solve(cho_factor(energy[np.ix_(unknown, unknown)]), load[unknown])
                mass = WATER_DENSITY*load[unknown].T@response
                if reference is None:
                    reference = energy, mass
                scale = np.r_[np.ones(3), np.full(3, 7.3152)]
                normalized_error = (mass-reference[1])/scale[:, None]/scale[None, :]
                normalized_reference = reference[1]/scale[:, None]/scale[None, :]
                row.update(complete=True, median_assembly_seconds=float(np.median(seconds)),
                    assembly_seconds=seconds, diagnostics=assembly.diagnostics,
                    added_mass_world=mass.tolist(),
                    diagonal_relative_error=(np.diag(mass-reference[1])/np.diag(reference[1])).tolist(),
                    normalized_matrix_relative_error=float(np.linalg.norm(normalized_error)/np.linalg.norm(normalized_reference)),
                    energy_relative_error=float(np.linalg.norm(energy-reference[0])/np.linalg.norm(reference[0])))
                saved[name+'_mass'], saved[name+'_energy'] = mass, energy
                del assembly
            except Exception as error:
                row['failure'] = dict(type=type(error).__name__, message=str(error))
                if name == 'uniform12':
                    raise RuntimeError('uniform order-12 reference failed; comparisons unavailable') from error
            print(json.dumps(row), flush=True)
            np.savez(output.with_suffix('.npz'), **saved)
            output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        if pressure_check or moving_check:
            from .boat_pressure_audit import run as pressure_audit
            original_prepare = gpu.prepare_energy
            audits = {}
            report['pressure_checks'] = {}
            report['moving_checks'] = {}
            try:
                for strategy in ('uniform12', 'slender10'):
                    calls = []

                    def current_prepare(current_faces, current_panels, order=3):
                        if strategy == 'uniform12':
                            return original_prepare(current_faces, current_panels, order=12)
                        samples, weights = triangle_samples(current_panels.triangles, high_order)
                        current_high = (gpu.cp.asarray(np.ascontiguousarray(samples)),
                            gpu.cp.asarray(weights), np.int32(len(weights)))
                        current_aspect = current_panels.edge_lengths.max(axis=1)**2/(2*current_panels.areas)
                        selected = current_aspect > 10
                        current_mask = gpu.cp.asarray(selected.astype(np.uint8))

                        def current_dispatch(grid, block, args):
                            kernel(grid, block, (*args, *current_high, current_mask))

                        gpu._pairs = current_dispatch
                        try:
                            value = original_prepare(current_faces, current_panels, order=3)
                        finally:
                            gpu._pairs = original
                        calls.append(dict(faces=len(current_faces), high_rows=int(selected.sum()),
                            minimum_threshold_distance=float(np.min(abs(current_aspect-10))),
                            h_constant_norm=value.diagnostics['h_constant_norm']))
                        return value

                    gpu.prepare_energy = current_prepare
                    if pressure_check:
                        path = output.with_name(output.stem+'-pressure-'+strategy+'.json')
                        audits[strategy] = pressure_audit(path, step=.001, refinement=2, refinement_region='body')
                        report['pressure_checks'][strategy] = dict(path=str(path),
                            sha256=sha256(path.read_bytes()).hexdigest(), integration_calls=list(calls),
                            derivatives=audits[strategy]['derivatives'],
                            directional_step_check=audits[strategy]['directional_step_check'])
                    if moving_check:
                        from .boat_moving import run as moving_pilot
                        calls.clear()
                        path = output.with_name(output.stem+'-moving-'+strategy+'.json')
                        pilot = moving_pilot(path, duration=.2, refinement=2, max_rhs=25,
                                             budget_seconds=90., max_step=.1)
                        if not pilot['complete'] or pilot['accepted_simulated_seconds'] != .2:
                            raise RuntimeError('moving comparison requires both complete 0.2 s endpoints')
                        with np.load(path.with_suffix('.npz'), allow_pickle=False) as data:
                            if data['accepted_times'][-1] != .2:
                                raise RuntimeError('moving state checkpoint lacks the target endpoint')
                            saved[strategy+'_moving_times'] = data['accepted_times']
                            saved[strategy+'_moving_states'] = data['accepted_states']
                        report['moving_checks'][strategy] = dict(path=str(path),
                            sha256=sha256(path.read_bytes()).hexdigest(), integration_calls=list(calls),
                            **{key: pilot[key] for key in ('seconds', 'accepted_simulated_seconds',
                                'accepted_steps', 'cold_geometry_count', 'controls', 'source_sha256')},
                            rhs_seconds=[row['seconds'] for row in pilot['rhs']],
                            energy_seconds=[sum(field['energy_seconds'] for field in row['fields'])
                                            for row in pilot['rhs']])
                    output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
                if pressure_check:
                    report['pressure_load_difference_world'] = (
                        np.array(audits['slender10']['derivatives'][-1]['hydrodynamic_load_world'])
                        -np.array(audits['uniform12']['derivatives'][-1]['hydrodynamic_load_world'])).tolist()
                    report['pressure_acceleration_difference_world'] = (
                        np.array(audits['slender10']['derivatives'][-1]['acceleration_world'])
                        -np.array(audits['uniform12']['derivatives'][-1]['acceleration_world'])).tolist()
                if moving_check:
                    from scipy.spatial.transform import Rotation
                    reference_state = saved['uniform12_moving_states'][-1]
                    candidate_state = saved['slender10_moving_states'][-1]
                    report['moving_endpoint_state_difference'] = (
                        candidate_state-reference_state).tolist()
                    report['moving_endpoint_rotation_difference_rad'] = float((
                        Rotation.from_quat(candidate_state[3:7], scalar_first=True)
                        *Rotation.from_quat(reference_state[3:7], scalar_first=True).inv()).magnitude())
            finally:
                gpu.prepare_energy = original_prepare
        report['complete'] = True
    finally:
        gpu._pairs = original
        np.savez(output.with_suffix('.npz'), **saved)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--geometry', required=True)
    parser.add_argument('--high-order', type=int, choices=(6, 12), default=6)
    parser.add_argument('--pressure-check', action='store_true')
    parser.add_argument('--moving-check', action='store_true')
    run(**vars(parser.parse_args()))
