"""Exact current BEM factor sharing on the actual boat; no trajectory timing.

Prescribed six-axis body velocity and observed-sea potential supply a nonzero
mixed field on the flat-water snapshot. The derivative check holds its nodal
potential fixed, as required by the shape kernel. It is NOT a pressure audit of
the proposed moving-world basis or a consistent nonlinear wave initialization.
"""
import argparse
from dataclasses import replace
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import HullMesh
from topwater.rankine import TrianglePanels
from topwater.sea_states import realize
from tests.laplace_reference import triangle_samples
from tests.reference_boat import BEAM, GRAVITY, WATER_DENSITY
from .boat_contact_preflight import make_chart, pose, target_heights
from .canoe import body_load, solve_mixed
from .incident_baseline_benchmark import coefficients, fused_cpu
from .p4_fluid_ordering import fresh


def run(output):
    from .cuda_operator import prepare_energy, cp
    from .cuda_shape_gradient import vertex_gradient

    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    report = dict(complete=False, qualified=False,
        scope='Current flat-water boat field/shape kernel only; no nonlinear trajectory or world-basis force qualification',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_current_factor.py', 'cuda_operator.py', 'cuda_shape_gradient.py',
                         'cuda_shape.cu', 'cuda_dual.cuh', 'cuda_pair.cuh', 'cuda_energy.cu', 'canoe.py')},
        density=WATER_DENSITY, gravity=GRAVITY, timings=[])
    arrays = {}

    def budget():
        if perf_counter()-start > 120:
            raise RuntimeError('120 s soft budget exhausted between operations')

    def timed(name, function):
        budget()
        before = perf_counter()
        result = function()
        cp.cuda.get_current_stream().synchronize()
        report['timings'].append(dict(name=name, seconds=perf_counter()-before))
        return result

    try:
        body, rest, chart = make_chart()
        points, faces = chart.reference_points, chart.faces
        mesh = HullMesh(points, faces)
        panels = TrianglePanels(points[faces], 1)
        report.update(body_mass=body.mass, vertices=len(points), triangles=len(faces))
        current = timed('current_energy_and_factor', lambda: prepare_energy(faces, panels))
        report['retained_factor_and_trace_bytes'] = current.factor.nbytes+current.mixed_trace.nbytes
        waves = realize(beam_m=BEAM, gravity=GRAVITY, seed=42)
        values = coefficients(waves)
        incident = fused_cpu(points, values, .37)
        twist = np.array([.1, .1, .1, .01, .01, .01])
        rhs = body_load(mesh, panels, chart.body_faces, rest.position)@twist
        # Prescribed incoming normal velocity at the finite walls/bed, not an
        # open-ocean boundary claim. Free-surface trace fixes the potential gauge.
        walls = chart.wall_faces
        samples, weights = triangle_samples(points[faces[walls]], 3)
        bary, _ = triangle_samples(np.eye(3)[None], 3)
        velocity = fused_cpu(samples.reshape(-1, 3), values, .37)[:, 6:].reshape(samples.shape)
        flux = np.einsum('fqi,fi->fq', velocity, panels.normals[walls])
        for corner in range(3):
            np.add.at(rhs, faces[walls, corner], panels.areas[walls]
                      *np.einsum('q,fq->f', weights*bary[0, :, corner], flux))
        potential = timed('current_mixed_field', lambda: solve_mixed(
            current.energy, chart.free_nodes, incident[chart.free_nodes, 4], rhs))
        unknown = np.setdiff1d(np.arange(len(points)), chart.free_nodes)
        residual = current.energy@potential-rhs
        report['field_relative_residual'] = float(np.linalg.norm(residual[unknown])/np.linalg.norm(rhs[unknown]))
        if report['field_relative_residual'] > 1e-9:
            raise ValueError('current mixed field residual failed')
        gradients = []
        for repetition in range(2):
            raw = timed('fresh_force_'+str(repetition), lambda: vertex_gradient(
                points, faces, potential, density=WATER_DENSITY))
            reused = timed('retained_force_'+str(repetition), lambda: vertex_gradient(
                points, faces, potential, assembly=current, density=WATER_DENSITY))
            gradients.append(reused)
            difference = float(np.max(abs(raw-reused)))
            report.setdefault('force_comparisons', []).append(dict(
                maximum_absolute_difference=difference,
                relative_l2_difference=float(np.linalg.norm(raw-reused)/np.linalg.norm(raw))))
            np.testing.assert_allclose(reused, raw, rtol=2e-10, atol=2e-9)
        legacy_density = timed('retained_force_density_1000', lambda: vertex_gradient(
            points, faces, potential, assembly=current, density=1000.))
        report['density_scaling_relative_error'] = float(
            np.linalg.norm(gradients[-1]-WATER_DENSITY/1000*legacy_density)/np.linalg.norm(gradients[-1]))
        if report['density_scaling_relative_error'] > 2e-12:
            raise ValueError('shape force density scaling failed')
        eta_direction = target_heights(points[chart.free_nodes], waves)
        xy = points[chart.free_nodes, :2]
        half_width = np.max(abs(xy), axis=0)
        eta_direction *= np.prod(np.maximum(0., 1-(xy/half_width)**2), axis=1)
        eta_direction *= .1/np.max(abs(eta_direction))
        base = chart.evaluate(rest.position, rest.orientation, np.zeros(len(chart.free_nodes)))
        direction = base.velocity(eta_direction, twist)
        if np.max(abs(direction[chart.wall_nodes])) > 1e-13:
            raise ValueError('shape derivative control moved its outer boundary')
        expected = float(np.sum(gradients[-1]*direction))
        if abs(expected) < 1e-6:
            raise ValueError('shape derivative control has negligible contraction')
        report['energy_directional_checks'] = []
        for h in (2e-5, 1e-5):
            energy = []
            for sign in (-1, 1):
                position, orientation = pose(rest, sign*h*twist[:3], sign*h*twist[3:])
                trial = chart.evaluate(position, orientation, sign*h*eta_direction).points
                assembled = timed('directional_energy_'+str(sign*h), lambda: prepare_energy(
                    faces, TrianglePanels(trial[faces], 1)))
                energy.append(float(.5*WATER_DENSITY*potential@assembled.energy@potential))
            actual = (energy[1]-energy[0])/(2*h)
            report['energy_directional_checks'].append(dict(step=h, finite_difference=actual,
                force_contraction=expected, relative_error=abs(actual-expected)/max(abs(actual), abs(expected), 1e-30)))
        if report['energy_directional_checks'][-1]['relative_error'] > 2e-4:
            raise ValueError('current energy shape derivative failed')
        report['invalidation_checks'] = []
        changed = points.copy()
        changed[faces[0, 0], 0] += 1e-8
        for name, assembly, coordinates, numbering in (
            ('moved_geometry', current, changed, faces),
            ('renumbered_faces', current, points, faces[::-1]),
            ('changed_order', replace(current, quadrature_order=4), points, faces),
        ):
            try:
                vertex_gradient(coordinates, numbering, potential, assembly=assembly, density=WATER_DENSITY)
            except ValueError:
                report['invalidation_checks'].append(name)
            else:
                raise ValueError('stale setup incorrectly accepted: '+name)
        arrays.update(points=points, faces=faces, potential=potential,
                      gradient=gradients[-1], shape_direction=direction, body_twist=twist)
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-start
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), indent=2))
