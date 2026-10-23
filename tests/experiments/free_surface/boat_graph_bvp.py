"""Matrix-free GPU graph/body BVP on one actual nonflat hull snapshot.

No full effective matrix or six-mode response basis is assembled. GMRES applies
the current graph lift and uses a flat source operator only as a preconditioner.
Two rigid-mode columns are solved solely as verification right-hand sides.
This is not a six-DOF trajectory or a baseline promotion.
"""
import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, gmres

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels, panel_influence
from tests.benchmark_rankine import peak_rss_bytes
from .body_internal_field import internal_field
from .boat_graph_images import image_height
from .boat_hos_continuation import tilted_geometry, tangent_graph
from .boat_image_plane import elevated_potential
from .boat_radiation_time import harmonic_maps
from .boat_spatial_work import source_matrix, subdivide
from .boat_parent_clip_proof import coplanar_groups
from .cuda_graph import GraphDirichlet, body_maps, cp, _resize
from .boat_source_projection import project_trace


@dataclass(frozen=True)
class GraphSolution:
    potential: np.ndarray
    density: np.ndarray
    spectrum: object
    diagnostics: dict


class GraphBody:
    def __init__(self, panels, eta, length, origin, reference_height, image_plane_height, order):
        started = perf_counter()
        self.panels, self.shape = panels, eta.shape
        self.length = length
        self.origin = np.array(origin, dtype=float, copy=True)
        self.origin.setflags(write=False)
        self.reference_height, self.image_height = reference_height, image_plane_height
        shifted = TrianglePanels(panels.triangles-[0., 0., image_plane_height], -1)
        self.matrix = cp.asarray(source_matrix(shifted))
        self.single = cp.asarray(assemble_dirichlet_halfspace(shifted).single_layer)
        preconditioner_panels = TrianglePanels(panels.triangles-[0., 0., reference_height], -1)
        self.preconditioner = lu_factor(source_matrix(preconditioner_panels))
        self.source_operator_seconds = perf_counter()-started
        before = perf_counter()
        self.graph = GraphDirichlet(eta-reference_height, [length, length], order=order)
        target = panels.centroids.copy()
        target[:, 2] -= reference_height
        padded = self.graph.padded_shape[0]
        self.potential, self.normal = body_maps(target, panels.normals, padded, length/padded, origin)
        cp.cuda.get_current_stream().synchronize()
        self.lift_map_seconds = perf_counter()-before
        before = perf_counter()
        size = self.shape[0]
        x, y = np.meshgrid(origin[0]+np.arange(size)*length/size,
                           origin[1]+np.arange(size)*length/size)
        self.surface_points = np.column_stack((x.ravel(), y.ravel(), eta.ravel()))
        self.surface_points.setflags(write=False)
        self.source_size = 2*size-1
        source_eta = cp.asnumpy(_resize(cp.asarray(eta), (self.source_size, self.source_size)))
        sx, sy = np.meshgrid(origin[0]+np.arange(self.source_size)*length/self.source_size,
                             origin[1]+np.arange(self.source_size)*length/self.source_size)
        self.source_points = np.column_stack((sx.ravel(), sy.ravel(), source_eta.ravel()))
        self.source_points.setflags(write=False)
        trace = elevated_potential(self.source_points, panels, image_plane_height)
        self.trace = self.project_source_trace(trace)
        cp.cuda.get_current_stream().synchronize()
        self.source_trace_seconds = perf_counter()-before
        self.evaluations = 0
        count = len(panels)
        self.operator = LinearOperator((count, count), matvec=self.apply, dtype=np.float64)
        self.precondition = LinearOperator((count, count),
            matvec=lambda value: lu_solve(self.preconditioner, value), dtype=np.float64)
        self.setup_seconds = perf_counter()-started

    def project_source_trace(self, values):
        return project_trace(values, self.source_size, self.shape[0])

    def spectrum(self, trace):
        spectra = self.graph.coefficients(trace.reshape(self.shape))
        return sum(spectra).ravel()

    def apply(self, value):
        self.evaluations += 1
        density = cp.asarray(value)
        spectrum = self.spectrum(self.trace@density)
        answer = cp.asnumpy(self.matrix@density-(self.normal@spectrum).real)
        if not np.isfinite(answer).all():
            raise ValueError('nonfinite current body-operator result')
        return answer

    def right_hand_side(self, trace, normal_load):
        lift = self.spectrum(cp.asarray(trace).reshape(-1))
        return np.asarray(normal_load)-cp.asnumpy((self.normal@lift).real)

    def solve_density(self, right):
        started, initial_evaluations = perf_counter(), self.evaluations
        residual_history = []
        density, info = gmres(self.operator, right, M=self.precondition,
            rtol=1e-9, atol=0., restart=min(40, len(right)), maxiter=5,
            callback=lambda value: residual_history.append(float(value)), callback_type='pr_norm')
        defect = float(np.linalg.norm(self.apply(density)-right))
        norm = float(np.linalg.norm(right))
        residual = defect/norm if norm else defect
        if info != 0 or not np.isfinite(residual) or residual > 2e-9:
            raise RuntimeError(f'current body GMRES failed: info={info}, true_relative_residual={residual}')
        return density, dict(gmres_iterations=len(residual_history),
            operator_applications=self.evaluations-initial_evaluations,
            true_relative_residual=residual, preconditioned_residual_history=residual_history,
            solve_seconds=perf_counter()-started)


    def solve(self, trace, normal_load):
        started = perf_counter()
        density, diagnostics = self.solve_density(self.right_hand_side(trace, normal_load))
        spectrum = self.spectrum(cp.asarray(trace).reshape(-1)-self.trace@cp.asarray(density))
        potential = self.body_potential(density, spectrum)
        cp.cuda.get_current_stream().synchronize()
        diagnostics['solve_seconds'] = perf_counter()-started
        return GraphSolution(potential, density, spectrum, diagnostics)


    def body_potential(self, density, spectrum):
        return cp.asnumpy((self.potential@spectrum).real+self.single@cp.asarray(density))


def tangent_plane_reference(panels, cg, level, slope, normal_load):
    """Independent aperiodic flat-plane BVP in a rotated coordinate frame.

    Same P0 body discretization, but no HOS, graph sampling or raised images.
    The plane continues to infinity; the tapered graph differs far from hull.
    This checks coupling plus far-graph differences, not continuum accuracy.
    """
    normal = np.r_[-slope, 1.]
    normal /= np.linalg.norm(normal)
    tangent = np.array([1., 0., slope[0]])
    tangent /= np.linalg.norm(tangent)
    transform = np.vstack((tangent, np.cross(normal, tangent), normal))
    origin = np.r_[cg[:2], level]
    triangles = (panels.triangles-origin)@transform.T
    # Undo only positive rotation roundoff on the known planar contact boundary.
    # No clearance is introduced; a resolved above-plane point is an error.
    roundoff = 64*np.finfo(float).eps*max(1., float(np.max(abs(panels.triangles-origin))))
    maximum_height = float(triangles[..., 2].max())
    if maximum_height > roundoff:
        raise ValueError('tangent-plane reference has resolved above-plane vertices')
    triangles[..., 2] = np.minimum(triangles[..., 2], 0.)
    rotated = TrianglePanels(triangles, -1)
    density = np.linalg.solve(source_matrix(rotated), normal_load)
    return (assemble_dirichlet_halfspace(rotated).single_layer@density,
            dict(maximum_vertex_height_before_roundtrip_repair_m=maximum_height,
                 roundtrip_roundoff_bound_m=roundoff))


def run(inputs, args, report):
    started = perf_counter()
    mesh = inputs.body.mesh
    coplanar_groups(mesh.vertices, mesh.faces, mesh.artificial_faces)
    coarse, level, slope = tilted_geometry(inputs)
    fixture = internal_field(coarse)  # Fixed before changing body discretization.
    triangles = coarse
    for _ in range(args.body_refinement-1):
        triangles = subdivide(triangles)
    panels = TrianglePanels(triangles, -1)
    cg = inputs.initial_state.position
    modes = np.column_stack((panels.normals, np.cross(panels.centroids-cg, panels.normals)))
    weighted = panels.areas[:, None]*modes
    x, y = np.meshgrid(cg[0]-args.length/2+np.arange(args.grid)*args.length/args.grid,
                       cg[1]-args.length/2+np.arange(args.grid)*args.length/args.grid)
    eta, _, outer = tangent_graph(inputs, x, y, level, slope)
    if np.any(outer >= args.length/2):
        raise ValueError('graph support must fit strictly inside the box')
    origin = cg[:2]-args.length/2
    height, bound = image_height(triangles, eta, [args.length, args.length], origin)
    reference_height = float(triangles[..., 2].max())
    true_plane = level+(triangles[..., :2]-cg[:2])@slope
    height = max(height, .5*(reference_height+float(true_plane.max())))
    plane_points = panels.centroids.copy()
    plane_points[:, 2] = 0
    eta_map, _, _ = harmonic_maps(plane_points, np.zeros_like(plane_points),
                                 args.grid, args.length/args.grid, origin)
    represented_eta = (eta_map@np.fft.fft2(eta).ravel()).real
    minimum_depth = float(np.min(represented_eta-panels.centroids[:, 2]))
    if minimum_depth <= 0:
        raise ValueError('body centroids must be below the represented Fourier graph')
    image_gap = float(np.min(2*height-triangles[..., 2]-true_plane))
    if image_gap < 0:
        raise ValueError('image panels intersect the true tangent-plane fluid')
    report.update(input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        panel_sha256=sha256(triangles.tobytes()).hexdigest(), body_panels=len(panels),
        grid_size=args.grid, box_length_m=args.length, hos_order=args.order,
        body_refinement=args.body_refinement, reference_height_m=reference_height,
        image_height_m=height, local_image_bound=bound,
        minimum_image_gap_above_true_plane_m=image_gap,
        minimum_centroid_depth_below_represented_graph_m=minimum_depth,
        represented_minus_true_graph_at_centroids_max_m=float(np.max(np.abs(
            represented_eta-(level+(panels.centroids[:, :2]-cg[:2])@slope)))),
        sources=dict(points=fixture.points.tolist(), strengths=fixture.strengths.tolist(),
                     center=fixture.center.tolist(), inscribed_radius=fixture.inscribed_radius))
    preparation_seconds = perf_counter()-started
    context = GraphBody(panels, eta, args.length, origin, reference_height, height, args.order)
    report['setup'] = dict(geometry_and_fields_seconds=preparation_seconds,
        source_operators_seconds=context.source_operator_seconds,
        lift_map_seconds=context.lift_map_seconds, source_trace_seconds=context.source_trace_seconds,
        context_seconds=context.setup_seconds)
    points = np.column_stack((x.ravel(), y.ravel(), eta.ravel()))
    trace = fixture.field(points)[0]
    exact, gradient = fixture.field(panels.centroids)
    right = np.einsum('fci,fi->fc', gradient, panels.normals)
    integrated_load = fixture.pressure_load(panels, cg, WATER_DENSITY)
    exact_centroid_load = -WATER_DENSITY*weighted.T@exact
    exact_flux = panel_influence(fixture.points, panels).double_layer.sum(axis=1)@fixture.strengths
    before = perf_counter()
    plane_potential, plane_geometry = tangent_plane_reference(panels, cg, level, slope, modes[:, [2, 3]])
    report['tangent_plane_reference_geometry'] = plane_geometry
    report['setup']['verification_tangent_plane_seconds'] = perf_counter()-before
    names = ['internal vertical dipole', 'internal transverse dipole', 'unit heave', 'unit world-x rotation']
    for col, name in enumerate(names):
        manufactured = col < 2
        psi = trace[:, col] if manufactured else np.zeros(args.grid**2)
        q = right[:, col] if manufactured else modes[:, [2, 3][col-2]]
        solution = context.solve(psi, q)
        potential, diagnostics = solution.potential, solution.diagnostics
        load = -WATER_DENSITY*weighted.T@potential
        row = dict(name=name, pressure_centroid_functional=load.tolist(), **diagnostics)
        if manufactured:
            row.update(body_potential_relative_l2=float(np.sqrt(
                np.sum(panels.areas*(potential-exact[:, col])**2)/np.sum(panels.areas*exact[:, col]**2))),
                exact_integrated_pressure_load=integrated_load[:, col].tolist(),
                centroid_functional_error=(load-exact_centroid_load[:, col]).tolist(),
                centroid_quadrature_error=(exact_centroid_load[:, col]-integrated_load[:, col]).tolist(),
                total_integrated_load_error=(load-integrated_load[:, col]).tolist(),
                prescribed_body_flux_centroid=float(panels.areas@q),
                exact_body_flux=float(exact_flux[col]))
        else:
            row['instantaneous_added_mass_column'] = (-load).tolist()
            reference_load = -WATER_DENSITY*weighted.T@plane_potential[:, col-2]
            row['tangent_plane_added_mass_column'] = (-reference_load).tolist()
            row['difference_from_tangent_plane_added_mass'] = (reference_load-load).tolist()
        report['rows'].append(row)
        print(json.dumps(row), flush=True)
    report.update(complete=True, device=cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)['name'].decode(),
        cupy_version=cp.__version__, cupy_pool_used_bytes=cp.get_default_memory_pool().used_bytes(),
        cupy_pool_reserved_bytes=cp.get_default_memory_pool().total_bytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--grid', type=int, default=65)
    parser.add_argument('--length', type=float, default=64.)
    parser.add_argument('--order', type=int, choices=range(2, 9), default=4)
    parser.add_argument('--body-refinement', type=int, choices=[1, 2, 3], default=1)
    args = parser.parse_args()
    if args.output.exists() or args.grid < 17 or args.grid % 2 != 1 or not np.isfinite(args.length) or args.length <= 0:
        parser.error('new output, odd grid>=17 and positive finite length required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path(__file__),
        *[Path('tests')/name for name in ['laplace_reference.py', 'benchmark_rankine.py', 'reference_hull.py']],
        *[Path(__file__).with_name(name) for name in ['body_internal_field.py', 'boat_hos_continuation.py',
          'boat_graph_images.py', 'boat_image_plane.py', 'boat_radiation_time.py', 'boat_spatial_work.py',
          'hos_graph.py', 'cuda_graph.py', 'cuda_operator.py', 'cuda_pair.cuh', 'cuda_energy.cu',
          'incident_baseline_benchmark.py', 'panel_gradient_reference.py',
          'boat_parent_clip_proof.py', 'p4_fluid_ordering.py']]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    report = dict(complete=False, qualified=False, scope=__doc__, rows=[], producer_sha256=hashes,
        limits=['Static current-graph BVP only; no moving mesh, pressure time derivative or sixDOF advancement',
                'Full manufactured auxiliary trace, including inside hull; no evolving hidden water states',
                'Preconditioner is a flat reference, not the current physical operator',
                'Aperiodic source fields and periodic HOS require independent domain control',
                'Point-source body data can be underresolved; centroid versus exact load/flux diagnostics expose this',
                'Tangent-plane control shares P0 hull discretization and differs from the tapered graph far from the hull',
                'Centroid depth checks do not certify whole panels below the Fourier graph; wet cut is to the true tangent plane',
                'Bulk layer evaluator still rejects exact edge/vertex target coincidence; no grid shift or physical gap is introduced'])
    started = perf_counter()
    try:
        run(SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8')), args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        report['peak_process_rss_bytes'] = peak_rss_bytes()
        report['memory_scope'] = 'Process peak RSS; ending CUDA pool is not total or peak VRAM.'
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
