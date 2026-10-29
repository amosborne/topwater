"""Compare an opt-in potential-only GPU evaluator with unmodified Capytaine.

uv run --locked --with pyopencl python scripts/prototype_gpu_surface.py
Use --points 5911 --frequencies 1920 for the full field workload. No public
example is overwritten. Requires an existing fp64 OpenCL GPU driver.
"""
import argparse
import json
from pathlib import Path
import tempfile
from time import perf_counter

import numpy as np
from threadpoolctl import threadpool_limits

from topwater.solve.surface import exterior_trace_points
from topwater.solve.fourier import surface_radiation_frequencies
from topwater.solve.wecsim.exchange import build_hydrodynamics
from topwater.solve.preparation import prepare_case
from topwater.inputs import WaveComponents
from topwater.solve.gpu import OpenCLPotential
from topwater.workflow_data import Boat, Sea


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--points', type=int, default=128)
    parser.add_argument('--frequencies', type=int, default=128)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--device', default='gfx1201')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--reconstruct', action='store_true',
                        help='Also reconstruct the entire saved case with the GPU adapter and compare the published CPU result')
    parser.add_argument('--reconstruct-cpu', action='store_true',
                        help='Also reconstruct with stock CPU radiation evaluation using --threads')
    parser.add_argument('--reconstruct-only', action='store_true',
                        help='Skip the paired matrix benchmark and run the selected reconstruction only')
    parser.add_argument('--potential-only-cpu', action='store_true',
                        help='Compatibility flag; CPU reconstruction always uses potential-only kernels')
    parser.add_argument('--gpu-fields',action='store_true',
                        help='Compatibility flag; GPU reconstruction always uses all GPU field kernels')
    args = parser.parse_args()
    if min(args.points, args.frequencies, args.threads) < 1:
        parser.error('counts must be positive')
    if args.reconstruct_only and not (args.reconstruct or args.reconstruct_cpu):
        parser.error('--reconstruct-only requires --reconstruct or --reconstruct-cpu')
    if args.potential_only_cpu and not args.reconstruct_cpu:
        parser.error('--potential-only-cpu requires --reconstruct-cpu')
    if args.gpu_fields and not args.reconstruct:
        parser.error('--gpu-fields requires --reconstruct')
    root = Path(__file__).resolve().parents[1]
    case = root/'examples/simple-boat-native'
    sea = Sea.load(case/'sea.npz')
    waves = WaveComponents(sea.amplitudes_m, sea.angular_frequencies_rad_s,
                           sea.directions_rad, sea.phases_rad)
    setup = prepare_case(Boat.load(case/'boat.npz').to_body(), waves, 100., rest_s=10., ramp_s=5.)
    output = args.output or Path(tempfile.mkdtemp(prefix='topwater-gpu-prototype-'))
    output.mkdir(parents=True, exist_ok=True)
    print(f'PROTOTYPE output={output}', flush=True)
    with tempfile.TemporaryDirectory(prefix='topwater-gpu-hydro-') as scratch:
        h = setup.hydrostatic
        last = max(34, int(np.ceil(setup.train.omega[-1]*10)))
        hydro = build_hydrodynamics(setup.wet_triangles_m, setup.mean_cg_m, setup.mass_matrix,
            Path(scratch)/'hydro', hydrostatic_stiffness=h.stiffness, cb=h.center_of_buoyancy,
            volume=h.displaced_volume, excitation_omega=np.arange(1,last+1)/10,
            radiation_omega=np.arange(1,max(120,last)+1)/10,
            field_omega=np.unique(np.r_[setup.train.omega[1:],surface_radiation_frequencies(setup.time_s)[1:]]),
            green_cache_dir=Path.home()/'.cache/topwater/capytaine-3.0.0')
        mesh = hydro.body.mesh_including_lid.merged()
        with np.load(case/'result.npz') as data:
            count = int(data['surface_vertex_offsets'][1])
            xy = (data['surface_vertices_m'][:count]@setup.world_from_native)[:, :2]
        mapped, _ = exterior_trace_points(hydro, xy, 0., allow_boundary=True)
        points = np.unique(mapped, axis=0)
        random = np.random.default_rng(42)
        if len(points) > args.points:
            points = points[random.choice(len(points), args.points, replace=False)]
        # Include surface panel centers, waterline vertices and edges explicitly.
        surface = mesh.vertices[abs(mesh.vertices[:, 2]) < 1e-12, :2]
        centers = mesh.faces_centers[abs(mesh.faces_centers[:, 2]) < 1e-12, :2]
        special = np.unique(np.vstack((surface, centers)), axis=0)
        points = np.c_[np.vstack((points, special)), np.zeros(len(points)+len(special))]
        omega = np.linspace(.025, 48., args.frequencies)
        ks = omega**2/hydro.gravity
        green = hydro.solver.engine.green_function
        if args.reconstruct_only and args.reconstruct_cpu:
            report = dict(cpu_threads=args.threads,
                reconstruction_backend='cpu-potential-only')
            reconstruct(hydro, setup, None, case, output, report,
                        cpu_threads=args.threads,potential_only_cpu=args.potential_only_cpu)
            return
        start = perf_counter()
        gpu = OpenCLPotential(green, mesh, device_name=args.device)
        setup_seconds = perf_counter()-start
        if args.reconstruct_only:
            report = dict(device=gpu.device.name, cpu_threads=args.threads,
                          reconstruction_backend='cpu' if args.reconstruct_cpu else 'gpu')
            reconstruct(hydro, setup, gpu, case, output, report,
                        cpu_threads=args.threads if args.reconstruct_cpu else None,
                        potential_only_cpu=args.potential_only_cpu,gpu_fields=args.gpu_fields)
            return
        errors = dict(max_abs=0., peak_reference=0., squared_error=0., squared_reference=0., values=0)
        cpu_seconds = gpu_seconds = gpu_kernel_seconds = 0.
        # Warm the driver/compiler independently of the measured steady-state work.
        gpu.evaluate(points[:2], ks[:2])
        with threadpool_limits(args.threads):
            for start in range(0, len(points), 128):
                chunk = points[start:start+128]
                began = perf_counter()
                actual = gpu.evaluate(chunk, ks)
                gpu_seconds += perf_counter()-began
                gpu_kernel_seconds += gpu.last_kernel_seconds
                began = perf_counter()
                expected = np.array([green.evaluate(chunk, mesh, wavenumber=float(k))[0] for k in ks])
                cpu_seconds += perf_counter()-began
                delta = abs(actual-expected)
                errors['max_abs'] = max(errors['max_abs'], float(delta.max()))
                errors['peak_reference'] = max(errors['peak_reference'], float(abs(expected).max()))
                errors['squared_error'] += float(np.sum(delta**2))
                errors['squared_reference'] += float(np.sum(abs(expected)**2))
                errors['values'] += actual.size
                print(f'PROTOTYPE points={min(start+128,len(points))}/{len(points)} '
                      f'max_abs={errors["max_abs"]:.6g}', flush=True)
        report = dict(device=gpu.device.name, driver=gpu.device.driver_version,
            capytaine_version='3.0.0', precision='float64', source_panels=mesh.nb_faces,
            quadrature_points=gpu.nq, frequencies=len(ks), field_points=len(points),
            cpu_threads=args.threads, setup_seconds=setup_seconds, build_seconds=gpu.build_seconds,
            cpu_seconds=cpu_seconds, gpu_seconds=gpu_seconds, gpu_kernel_seconds=gpu_kernel_seconds,
            speedup=cpu_seconds/gpu_seconds, **errors)
        report['relative_rms_error'] = np.sqrt(errors['squared_error']/errors['squared_reference'])
        report['relative_max_error'] = errors['max_abs']/errors['peak_reference']
        (output/'comparison.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(report, indent=2), flush=True)
        if report['relative_max_error'] > 1e-8 or report['relative_rms_error'] > 1e-9:
            raise RuntimeError('GPU potential did not meet the prototype equivalence limits')
        if args.reconstruct or args.reconstruct_cpu:
            reconstruct(hydro, setup, gpu, case, output, report,
                        cpu_threads=args.threads if args.reconstruct_cpu else None,
                        potential_only_cpu=args.potential_only_cpu,gpu_fields=args.gpu_fields)


def reconstruct(hydro, setup, gpu, case, output, report, *, cpu_threads=None,
                potential_only_cpu=False,gpu_fields=False):
    """Benchmark the production reconstruction pipeline without substitutions."""
    from dataclasses import replace
    from topwater.solve import surface as producer
    from topwater.solve.native.dynamics import integrate_motion
    from topwater.solve.field_backend import select_field_backend
    from topwater.solve.sampling import surface_sample_times
    motion = integrate_motion(setup, hydro)
    offset = setup.world_vertical_offset_m
    computational = replace(motion, position_m=motion.position_m-np.array([0.,0.,offset]))
    if gpu is not None:
        gpu.close()
    backend = select_field_backend(hydro, 'gpu' if cpu_threads is None else 'cpu')
    if cpu_threads is not None:
        backend.evaluator.threads = cpu_threads
    started = perf_counter()
    try:
        surface_time = surface_sample_times(motion.time_s)
        vertices,faces,vo,fo = producer.reconstruct_total_surface(
            hydro, setup, computational, surface_time, field_backend=backend)
    finally:
        backend.close()
    elapsed = perf_counter()-started
    vertices[:,2] += offset
    with np.load(case/'result.npz') as expected:
        np.testing.assert_array_equal(faces, expected['surface_faces'])
        np.testing.assert_array_equal(vo, expected['surface_vertex_offsets'])
        np.testing.assert_array_equal(fo, expected['surface_face_offsets'])
        delta = abs(vertices-expected['surface_vertices_m'])
        position_error = np.linalg.norm(motion.position_m-expected['position_m'],axis=1).max()
    report.update(reconstruction_seconds=elapsed,
                  gpu_incident_and_diffraction=cpu_threads is None,
                  reconstruction_backend='gpu' if cpu_threads is None else 'cpu-potential-only',
                  surface_vertex_max_abs_difference_m=float(delta.max()),
                  motion_max_position_difference_m=float(position_error),
                  surface_topology_identical=True)
    (output/'comparison.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('RECONSTRUCTION '+json.dumps(report,indent=2),flush=True)
    if delta.max() > 1e-8 or position_error > 1e-8:
        raise RuntimeError('GPU reconstruction differs from the frozen CPU field')


if __name__ == '__main__':
    main()
