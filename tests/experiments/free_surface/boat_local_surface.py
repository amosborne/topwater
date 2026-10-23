"""Actual nonflat hull: local versus uniform fluid-trace refinement.

Manufactured affine and five harmonic-quadratic acceleration potentials prescribe
exact fine-node Dirichlet data and weak Neumann data. The physical P1 graph is
unchanged. This deliberately bypasses coarse parent-state interpolation; it is
not a moving simulation, open-water boundary, or alternate baseline path.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from topwater import SimulationInputs, HullMesh
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from topwater.surface import body_load
from tests.benchmark_rankine import peak_rss_bytes
from tests.laplace_reference import triangle_samples
from .boat_world_boundary import construct
from .boat_world_field import refine_boundary


def field(points, center, beam):
    x, y, z = np.moveaxis((points-center)/beam, -1, 0)
    value = np.stack((.2*x+.3*y+.4*z, x*y, x*z, y*z, x*x-y*y, x*x-z*z), axis=-1)
    gradient = np.empty((*x.shape, 6, 3))
    gradient[..., 0, :] = np.array([.2, .3, .4])/beam
    zero = np.zeros_like(x)
    for column, components in enumerate(((y, x, zero), (z, zero, x), (zero, z, y),
                                         (2*x, -2*y, zero), (2*x, zero, -2*z)), start=1):
        gradient[..., column, :] = np.stack(components, axis=-1)/beam
    return value, gradient


def exact_data(points, faces, free_faces, body_faces, center, beam, order=3):
    panels = TrianglePanels(points[faces], 1)
    samples, weights = triangle_samples(panels.triangles, order)
    basis, _ = triangle_samples(np.eye(3)[None], order)
    value, gradient = field(samples, center, beam)
    flux = np.einsum('fqmi,fi->fqm', gradient, panels.normals)
    load = np.zeros((len(points), value.shape[-1]))
    nonfree = np.setdiff1d(np.arange(len(faces)), free_faces)
    for corner in range(3):
        contribution = panels.areas[:, None]*np.einsum('q,q,fqm->fm',
            weights, basis[0, :, corner], flux)
        np.add.at(load, faces[nonfree, corner], contribution[nonfree])
    normals = np.broadcast_to(panels.normals[:, None], samples.shape)
    modes = np.concatenate((normals, np.cross(samples-center, normals)), axis=-1)
    moments = -WATER_DENSITY*np.einsum('f,q,fqk,fqm->km', panels.areas[body_faces],
        weights, modes[body_faces], value[body_faces])
    energy = np.einsum('f,q,fqm,fqm->m', panels.areas, weights, value, flux)
    return load, moments, energy


def solve_case(arrays, center, beam, kind, level, order):
    from .cuda_operator import prepare_energy, cp
    start = perf_counter()
    points, faces = arrays['points'], arrays['faces']
    panels = TrianglePanels(points[faces], 1)
    exact, _ = field(points, center, beam)
    load, reference, exact_energy = exact_data(points, faces, arrays['free_faces'],
        arrays['body_faces'], center, beam)
    _, reference4, energy4 = exact_data(points, faces, arrays['free_faces'],
        arrays['body_faces'], center, beam, order=4)
    np.testing.assert_allclose(reference, reference4, rtol=2e-12, atol=2e-10)
    np.testing.assert_allclose(exact_energy, energy4, rtol=2e-12, atol=2e-12)
    before = perf_counter()
    current = prepare_energy(faces, panels, order=order)
    cp.cuda.get_current_stream().synchronize()
    assembly_seconds = perf_counter()-before
    fixed = arrays['free_nodes']
    unknown = np.setdiff1d(np.arange(len(points)), fixed)
    potential = np.zeros_like(exact)
    potential[fixed] = exact[fixed]
    right = (load-current.energy@potential)[unknown]
    before = perf_counter()
    factor = cho_factor(current.energy[np.ix_(unknown, unknown)])
    potential[unknown] = cho_solve(factor, right)
    solve_seconds = perf_counter()-before
    pressure_map = -WATER_DENSITY*body_load(HullMesh(points, faces), panels,
        arrays['body_faces'], center).T
    moment = pressure_map@potential
    interpolated_exact_moment = pressure_map@exact
    body = arrays['body_faces']
    samples, weights = triangle_samples(panels.triangles[body], 3)
    basis, _ = triangle_samples(np.eye(3)[None], 3)
    analytic, _ = field(samples, center, beam)
    interpolated = np.einsum('qc,fcm->fqm', basis[0], potential[faces[body]])
    norm2 = lambda values: np.einsum('f,q,fqm->m', panels.areas[body], weights, values*values)
    residual = (current.energy@potential-load)[unknown]
    return dict(kind=kind, refinement_level=level, kernel_quadrature_order=order,
        nodes=len(points), panels=len(faces), body_panels=len(body), free_panels=len(arrays['free_faces']),
        geometry_sha256=sha256(points.tobytes()+faces.tobytes()).hexdigest(),
        body_potential_relative_l2=np.sqrt(norm2(interpolated-analytic)/norm2(analytic)).tolist(),
        pressure_moments=moment.tolist(), exact_pressure_moments=reference.tolist(),
        pressure_moment_error=(moment-reference).tolist(),
        pressure_solver_error=(moment-interpolated_exact_moment).tolist(),
        pressure_body_interpolation_error=(interpolated_exact_moment-reference).tolist(),
        energy_relative_error=(np.einsum('im,im->m', potential, current.energy@potential)/exact_energy-1).tolist(),
        algebraic_relative_residual=(np.linalg.norm(residual, axis=0)/np.linalg.norm(right, axis=0)).tolist(),
        assembly_seconds=assembly_seconds, mixed_factor_solve_seconds=solve_seconds,
        total_case_seconds=perf_counter()-start, operator_diagnostics=current.diagnostics)


def run(inputs, output, body_refinement):
    started = perf_counter()
    geometry = output.with_name(output.stem+'-geometry.json')
    _, base, metadata = construct(geometry, body=inputs.body, pose=inputs.initial_state, waves=inputs.waves)
    beam = float(np.ptp(inputs.body.mesh.vertices[:, 1]))
    rows = []
    # Shared body enrichment; conforming free-edge splits can subsequently
    # enrich adjacent body/wall triangles too. Counts are reported explicitly.
    for _ in range(body_refinement):
        base = refine_boundary(base, region='body')
    cases = [('base', 0, base)]
    for kind in ['local_free', 'uniform_free']:
        current = base
        for level in [1, 2]:
            selected = None
            if kind == 'local_free':
                contact = np.intersect1d(current['free_nodes'],
                    np.unique(current['faces'][current['body_faces']]))
                free = current['free_faces']
                selected = free[np.any(np.isin(current['faces'][free], contact), axis=1)]
            current = refine_boundary(current, region='free', marked_faces=selected)
            cases.append((kind, level, current))
    for kind, level, arrays in cases:
        for order in [3, 6]:
            row = solve_case(arrays, inputs.initial_state.position, beam, kind, level, order)
            rows.append(row)
            print(json.dumps(row), flush=True)
    return dict(complete=True, qualified=False,
        scope=__doc__, input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        geometry_report_sha256=sha256(geometry.read_bytes()).hexdigest(),
        geometry_arrays_sha256=sha256(geometry.with_suffix('.npz').read_bytes()).hexdigest(),
        initial_geometry_seconds=metadata['seconds'],
        common_body_refinement=body_refinement,
        graph_height_range_m=np.ptp(base['parent_nodes'][:, 2]).item(),
        fields=['affine', 'xy', 'xz', 'yz', 'x^2-y^2', 'x^2-z^2'],
        force_rows=['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'],
        moment_units=['N']*3+['Nm']*3,
        field_interpretation='Unit-scaled acceleration potentials; loads equal minus rho integral of phi times fluid-normal rigid mode',
        refinement='Exact same physical P1 graph and hull; fine-node analytic Dirichlet traces, not prolonged 15-parent state',
        limits=['Finite diagnostic tank, not ocean radiation or a moving sixDOF case',
                'Conforming free-edge refinement also subdivides neighboring body/wall panels; counts expose this coupling',
                'No new production tolerance, integrator or runtime refinement switch'],
        rows=rows, peak_process_rss_bytes=peak_rss_bytes(), total_seconds=perf_counter()-started)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--body-refinement', type=int, choices=[1, 2, 3], default=1)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('new output required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path(__file__),
        *[Path('tests')/name for name in ['laplace_reference.py', 'benchmark_rankine.py', 'reference_hull.py']],
        *[Path(__file__).with_name(name) for name in ['boat_world_boundary.py', 'boat_parent_clip_proof.py',
          'boat_integer_intersections.py', 'boat_world_field.py', 'incident_baseline_benchmark.py',
          'p4_fluid_ordering.py', 'cuda_operator.py', 'cuda_energy.cu', 'cuda_pair.cuh']]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    result = run(SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8')),
                 args.output, args.body_refinement)
    result['producer_sha256'] = hashes
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
