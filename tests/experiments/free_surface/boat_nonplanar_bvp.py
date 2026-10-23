"""Bounded BVP control on owned curved-graph geometry reports.

One fixed internal source fixture supplies independent acceleration potentials
and exact six pressure moments. Harmonic exterior-total-trace fill and the known
full auxiliary trace are verification inputs to the SAME solve. No moving state,
time pressure, or qualified physical hull load is inferred from this probe.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull

from topwater import SimulationInputs
from topwater.constants import WATER_DENSITY
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from .body_internal_field import internal_field
from .boat_graph_bvp import GraphBody
from .boat_graph_images import image_height
from .boat_hos_continuation import tilted_geometry
from .boat_spatial_work import subdivide
from .cuda_graph import ExteriorTrace, cp
from .graph_geometry_reference import tensor_fourier, fourier_points


def run(parent_path, args, report):
    parent = json.loads(parent_path.read_text(encoding='utf-8'))
    if not parent['complete']:
        raise ValueError('complete curved-geometry parent report required')
    for name, digest in parent['producer_sha256'].items():
        if sha256(Path(name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'parent producer changed: {name}')
    inputs = SimulationInputs.from_json((parent_path.with_suffix('.source')/'inputs.json').read_text(encoding='utf-8'))
    if sha256(inputs.to_json().encode()).hexdigest() != parent['input_sha256']:
        raise ValueError('parent physical input identity changed')
    graph_path = parent_path.with_name(parent_path.stem+'-graph.npz')
    graph = np.load(graph_path)
    original = graph['eta']
    origin, length = graph['origin'], float(graph['length'])
    expected_shape = (parent['graph_grid_size'],)*2
    if (original.shape != expected_shape or original.dtype != np.float64
            or length != parent['graph_box_length_m']
            or not np.array_equal(origin, inputs.initial_state.position[:2]-length/2)):
        raise ValueError('parent Fourier graph shape, dtype, origin or length changed')
    if sha256(original.tobytes()).hexdigest() != parent['graph_sha256']:
        raise ValueError('parent Fourier graph changed')
    size = (original.shape[0]-1)*args.surface_refinement+1
    x = origin[0]+np.arange(size)*length/size
    y = origin[1]+np.arange(size)*length/size
    eta = tensor_fourier(original, x, y, length, origin)[0]
    xx, yy = np.meshgrid(x, y)
    points = np.column_stack((xx.ravel(), yy.ravel(), eta.ravel()))
    body, pose = inputs.body, inputs.initial_state
    cg = pose.position
    world = (body.mesh.vertices-body.cg)@rotation_matrix(pose.orientation).T+cg
    equations = ConvexHull(world).equations
    maximum = np.max(points@equations[:, :3].T+equations[:, 3], axis=1)
    roundoff = 64*np.finfo(float).eps*max(1., float(abs(world).max()), float(abs(points).max()))
    if np.any(abs(maximum) <= roundoff):
        raise ValueError('a graph node is numerically on the hull; exterior ownership is ambiguous')
    extension = ExteriorTrace((maximum < 0).reshape(size, size))
    # Freeze the manufactured poles BEFORE either curved geometry or its panels.
    fixture = internal_field(tilted_geometry(inputs)[0])
    if np.any(fixture.points@equations[:, :3].T+equations[:, 3] >= 0):
        raise ValueError('manufactured poles must remain strictly inside the physical solid')
    pole_gap = fourier_points(original, fixture.points[:, :2], length, origin)[0]-fixture.points[:, 2]
    if np.any(pole_gap <= 0):
        raise ValueError('this manufactured fixture must stay below the auxiliary graph')
    cases = []
    for row in parent['rows']:
        if args.roof_refinement is not None and row['roof_refinement'] not in args.roof_refinement:
            continue
        path = Path(row['geometry_report']).with_suffix('.npz')
        if sha256(path.read_bytes()).hexdigest() != row['geometry_npz_sha256']:
            raise ValueError('parent wet geometry archive changed')
        arrays = np.load(path)
        frame_origin = np.asarray(row['geometry_frame_origin_world_m'])
        if (not np.array_equal(frame_origin, np.r_[cg[:2], 0.])
                or not np.array_equal(arrays['geometry_frame_origin_world_m'], frame_origin)
                or row['body_triangles_coordinate_frame'] != 'world'):
            raise ValueError('parent geometry frame identity changed')
        local_triangles = arrays['points'][arrays['faces'][arrays['body_faces']]]
        if sha256(local_triangles.tobytes()).hexdigest() != row['local_body_triangles_sha256']:
            raise ValueError('parent local body triangles changed')
        fluid_triangles = local_triangles+frame_origin
        if sha256(fluid_triangles.tobytes()).hexdigest() != row['body_triangles_sha256']:
            raise ValueError('parent body triangle identity changed')
        # GraphBody owns BODY-outward triangles and supplies normal_sign=-1.
        triangles = fluid_triangles[:, [0, 2, 1]]
        for _ in range(args.body_refinement):
            triangles = subdivide(triangles)
        cases.append((row, TrianglePanels(triangles, -1)))
    if (not cases or args.roof_refinement is not None
            and set(args.roof_refinement) != {row['roof_refinement'] for row, _ in cases}):
        raise ValueError('requested roof inputs are absent from the owned parent report')
    required_image = max(image_height(p.triangles, eta, [length]*2, origin)[0] for _, p in cases)
    image = required_image if args.image_height is None else args.image_height
    if not np.isfinite(image) or image < required_image:
        raise ValueError('prescribed image plane does not meet the current graph bound')
    height = max(float(p.triangles[..., 2].max()) for _, p in cases)
    true_trace = fixture.field(points)[0]
    physical_trace = cp.asnumpy(extension.extend(true_trace[extension.exterior]))
    report.update(input_sha256=parent['input_sha256'], parent_graph_sha256=parent['graph_sha256'],
        parent_graph_npz_sha256=sha256(graph_path.read_bytes()).hexdigest(),
        surface_refinement=args.surface_refinement, body_refinement=args.body_refinement,
        roof_refinements=[row['roof_refinement'] for row, _ in cases],
        graph_grid_size=size, graph_box_length_m=length, hos_order=args.order,
        image_height_m=image, reference_height_m=height,
        minimum_required_image_height_m=required_image,
        exterior_nodes=len(extension.exterior), hidden_nodes=len(extension.hidden),
        source_fixture=dict(points=fixture.points.tolist(), strengths=fixture.strengths.tolist(),
            minimum_pole_depth_below_auxiliary_graph_m=float(pole_gap.min()),
            interpretation='Acceleration potentials; one unit m^3/s^2 per source/sink strength'))
    for geometry, panels in cases:
        case = dict(roof_refinement=geometry['roof_refinement'], body_panels=len(panels), rows=[],
            panel_sha256=sha256(panels.triangles.tobytes()).hexdigest())
        report['cases'].append(case)
        depth = fourier_points(original, panels.centroids[:, :2], length, origin)[0]-panels.centroids[:, 2]
        case['minimum_centroid_depth_below_fourier_graph_m'] = float(depth.min())
        if np.any(depth <= 0):
            raise ValueError('a body collocation point lies above the Fourier graph')
        context = GraphBody(panels, eta, length, origin, height, image, args.order)
        case['setup'] = dict(total_seconds=context.setup_seconds,
            source_operator_seconds=context.source_operator_seconds,
            harmonic_map_seconds=context.lift_map_seconds, source_trace_seconds=context.source_trace_seconds)
        modes = np.column_stack((panels.normals, np.cross(panels.centroids-cg, panels.normals)))
        weighted = panels.areas[:, None]*modes
        exact, gradient = fixture.field(panels.centroids)
        q = np.einsum('nci,ni->nc', gradient, panels.normals)
        integrated_load = fixture.pressure_load(panels, cg, WATER_DENSITY)
        centroid_load = -WATER_DENSITY*weighted.T@exact
        solutions = {}
        for col, name in enumerate(['vertical internal dipole', 'transverse internal dipole']):
            for label, data in [('physical exterior total fill', physical_trace[:, col]),
                                ('known analytic auxiliary trace', true_trace[:, col])]:
                solution = context.solve(data, q[:, col])
                load = -WATER_DENSITY*weighted.T@solution.potential
                metric = dict(field=name, auxiliary_trace=label, **solution.diagnostics,
                    body_potential_relative_l2=float(np.sqrt(np.sum(panels.areas*(solution.potential-exact[:, col])**2)
                        /np.sum(panels.areas*exact[:, col]**2))),
                    pressure_centroid_functional=load.tolist(),
                    exact_integrated_pressure_functional=integrated_load[:, col].tolist(),
                    solution_minus_exact_centroid_functional=(load-centroid_load[:, col]).tolist(),
                    centroid_quadrature_error=(centroid_load[:, col]-integrated_load[:, col]).tolist(),
                    total_integrated_pressure_functional_error=(load-integrated_load[:, col]).tolist())
                case['rows'].append(metric)
                solutions[col, label] = solution.potential
        case['extension_sensitivity'] = [dict(field=name,
            pressure_centroid_functional_change=(-WATER_DENSITY*weighted.T@(
                solutions[col, 'physical exterior total fill']-solutions[col, 'known analytic auxiliary trace'])).tolist())
            for col, name in enumerate(['vertical internal dipole', 'transverse internal dipole'])]
        for column, name in [(2, 'unit heave'), (3, 'unit world-x rotation')]:
            solution = context.solve(np.zeros(size*size), modes[:, column])
            case['rows'].append(dict(field=name, **solution.diagnostics,
                instantaneous_added_mass_column=(WATER_DENSITY*weighted.T@solution.potential).tolist()))
        print(json.dumps({k:v for k,v in case.items() if k != 'rows'}), flush=True)
        del context, solution
    report['between_roof_reference_changes'] = []
    for before, after in zip(report['cases'][:-1], report['cases'][1:]):
        fields = []
        for index in [0, 2]:
            a, b = before['rows'][index], after['rows'][index]
            change = np.asarray(b['exact_integrated_pressure_functional'])-a['exact_integrated_pressure_functional']
            fields.append(dict(field=a['field'], exact_integrated_pressure_functional_change=change.tolist()))
        report['between_roof_reference_changes'].append(dict(
            from_roof=before['roof_refinement'], to_roof=after['roof_refinement'], fields=fields))
    report.update(complete=True, cupy_pool_used_bytes=cp.get_default_memory_pool().used_bytes(),
        cupy_pool_reserved_bytes=cp.get_default_memory_pool().total_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry-report', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--body-refinement', type=int, choices=[0, 1, 2], default=1)
    parser.add_argument('--roof-refinement', type=int, choices=[1, 2], nargs='+', help='select owned input roof cases')
    parser.add_argument('--surface-refinement', type=int, choices=[1, 2, 4], default=1)
    parser.add_argument('--order', type=int, choices=range(2, 9), default=4)
    parser.add_argument('--image-height', type=float, help='fixed valid representation plane for controlled comparisons')
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output and archive required')
    parent = json.loads(args.geometry_report.read_text(encoding='utf-8'))
    archive = args.output.with_suffix('.source')
    shutil.copytree(args.geometry_report.with_suffix('.source'), archive)
    hashes = dict(parent['producer_sha256'])
    source = Path(__file__).resolve().relative_to(Path.cwd())
    shutil.copyfile(source, archive/source)
    hashes[str(source)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.geometry_report, archive/'parent-report.json')
    shutil.copyfile(args.geometry_report.with_name(args.geometry_report.stem+'-graph.npz'), archive/'graph.npz')
    for row in parent['rows']:
        shutil.copyfile(Path(row['geometry_report']).with_suffix('.npz'),
                        archive/f"roof{row['roof_refinement']}-geometry.npz")
    report = dict(complete=False, qualified=False, scope=__doc__, producer_sha256=hashes,
        parent_report_sha256=sha256(args.geometry_report.read_bytes()).hexdigest(), cases=[], limits=[
            'Static curved-graph BVP; no evolving surface, total time pressure or sixDOF trajectory',
            'P1 clipping differs from the fixed Fourier graph; parent geometry mismatch remains explicit',
            'Known auxiliary trace is an independent verification input, not data available to a simulation',
            'Source fixture is fixed across meshes; thin-hull internal data can be underresolved by P0 panels',
            'Source positions/fields are inside the physical solid, never singularities in physical fluid',
            'Aperiodic sources/periodic HOS, finite-frequency radiation and body-work convergence remain unqualified',
            'Ending CUDA pool is not total/peak VRAM; component timings are not nonlinear throughput'])
    started = perf_counter()
    try:
        run(args.geometry_report, args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        report['peak_process_rss_bytes'] = peak_rss_bytes()
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
