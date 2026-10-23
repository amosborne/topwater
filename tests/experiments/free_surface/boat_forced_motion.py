"""Free six-DOF body motion under a FIXED, FORCED graph/Dirichlet boundary.

This intermediate integration gate does not evolve a sea or enforce p=0 on it.
Fresh current fields and exact discrete pressure actions drive every body RHS.
Representation heights and exterior labels stay fixed; leaving them fails closed.
"""
import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.integrate import RK23
from scipy.spatial import ConvexHull

from topwater import InitialState, SimulationInputs, sea_states
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.frames import rotation_matrix
from tests.benchmark_rankine import peak_rss_bytes
from tests.forced_graph_case import boundary_data
from .boat_parent_clip_proof import coplanar_groups
from .boat_graph_geometry import constant_twist_pose, wet_snapshot, fourier_dry_certificate
from .boat_graph_images import image_height
from .graph_geometry_reference import fourier_points
from .boat_graph_bvp import GraphBody
from .boat_graph_fields import body_velocity
from .boat_frozen_field import frozen_field
from .boat_coupled_pressure import coupled_pressure, body_dynamics
from .boat_motion_metrics import attitude_degrees
from .cuda_graph import ExteriorTrace, cp


def state_from_vector(vector):
    vector = np.asarray(vector)
    if vector.shape != (13,) or not np.isfinite(vector).all() or np.linalg.norm(vector[3:7]) == 0:
        raise ValueError('finite13-vector with nonzero quaternion required')
    quaternion = vector[3:7]/np.linalg.norm(vector[3:7])
    rotation = rotation_matrix(quaternion)
    return InitialState(vector[:3], quaternion, vector[7:10], rotation.T@vector[10:13])


@dataclass
class Snapshot:
    state: InitialState
    panels: object
    mask: np.ndarray
    owners: np.ndarray
    signs: np.ndarray
    diagnostics: dict


@dataclass
class Evaluation:
    rate: np.ndarray
    pressure: np.ndarray
    acceleration: np.ndarray
    diagnostics: dict


class ForcedBody:
    """One guarded body RHS; boundary data and input hull are explicit owners."""
    def __init__(self, inputs, boundary, derivative_step=.005, refinement=1, order=4):
        self.body, self.boundary = inputs.body, boundary
        self.step, self.refinement, self.order = derivative_step, refinement, order
        coplanar_groups(self.body.mesh.vertices, self.body.mesh.faces, self.body.mesh.artificial_faces)
        beam = float(np.ptp(self.body.mesh.vertices[:, 1]))
        sea_states.validate_envelope(inputs.waves, beam_m=beam, gravity=GRAVITY)
        initial = self.stencil(inputs.initial_state, self.step)
        self.mask = initial[0].mask.copy()
        self.extension = ExteriorTrace(self.mask.reshape(boundary['eta'].shape))
        # Fixed representation over this bounded pilot. These margins are setup
        # padding, not physical operating limits or a refresh policy.
        margin = .01*beam
        self.height = max(s.panels.triangles[..., 2].max() for s in initial)+margin
        self.image = max(self.required_image(s.panels) for s in initial)+margin
        self.setup = dict(reference_height_m=float(self.height), image_height_m=float(self.image),
            representation_padding_m=margin, derivative_step_seconds=self.step,
            exterior_mask_sha256=sha256(self.mask.tobytes()).hexdigest(),
            policy='Fixed image/expansion heights and exterior mask; fail on invalidation, no automatic refresh',
            physical_scope='Freely moving body, fixed forced graph and analytic exterior trace; no evolving free surface')

    def required_image(self, panels):
        b = self.boundary
        graph_bound = image_height(panels.triangles, b['eta'], [b['length']]*2, b['origin'])[0]
        normal = b['frame'][2]
        slope = -normal[:2]/normal[2]
        plane_height = b['plane_origin'][2]+(panels.triangles[..., :2]-b['plane_origin'][:2])@slope
        plane_bound = float(np.max(panels.triangles[..., 2]+plane_height))/2
        return max(graph_bound, float(np.nextafter(plane_bound, np.inf)))

    def snapshot(self, state):
        started = perf_counter()
        b = self.boundary
        world = (self.body.mesh.vertices-self.body.cg)@rotation_matrix(state.orientation).T+state.position
        if np.any(np.max(abs(world[:, :2]-b['plane_origin'][:2]), axis=0) >= b['inner']):
            raise ValueError('body left the forced graph tangent-plane plateau')
        panels, mask, owners, signs = wet_snapshot(self.body, state, b['plane_origin'], b['frame'],
                                                  b['grid_xy'], b['length'], self.refinement)
        if hasattr(self, 'mask') and not np.array_equal(mask, self.mask):
            raise ValueError('exterior-mask event: pilot stops before harmonic-fill state changes')
        depth = fourier_points(b['eta'], panels.centroids[:, :2], b['length'], b['origin'])[0] \
            -panels.centroids[:, 2]
        if depth.min() <= 0:
            raise ValueError('body collocation point crossed represented Fourier graph')
        # Actual Fourier samples and the clipping plane must classify nodes alike.
        equations = ConvexHull(world).equations
        graph_points = np.column_stack((b['grid_xy'], b['eta'].ravel()))
        actual_mask = np.max(graph_points@equations[:, :3].T+equations[:, 3], axis=1) <= 0
        if not np.array_equal(mask, actual_mask):
            raise ValueError('Fourier graph and planar clipping disagree on exterior-node ownership')
        dry = fourier_dry_certificate(self.body, world, b['eta'], b['length'], b['origin'])
        flat = panels.triangles.reshape(-1, 3)
        contact = flat[abs((flat-b['plane_origin'])@b['frame'][2]) < 64*np.finfo(float).eps*b['length']]
        contact_height = fourier_points(b['eta'], contact[:, :2], b['length'], b['origin'])[0]
        if hasattr(self, 'height'):
            if panels.triangles[..., 2].max() > self.height or self.required_image(panels) > self.image:
                raise ValueError('body invalidated fixed expansion/image representation')
        return Snapshot(state, panels, mask, owners, signs,
            dict(total_seconds=perf_counter()-started, panels=len(panels),
                triangle_sha256=sha256(panels.triangles.tobytes()).hexdigest(),
                minimum_fourier_centroid_depth_m=float(depth.min()),
                waterline_vertex_fourier_mismatch_max_m=float(np.max(abs(contact[:, 2]-contact_height))),
                dry_fourier_certificate=dry))

    def stencil(self, state, step):
        snapshots = [self.snapshot(pose) for pose in
            [state, constant_twist_pose(state, step), constant_twist_pose(state, -step)]]
        center = snapshots[0]
        for offset in snapshots[1:]:
            if (len(offset.panels) != len(center.panels) or not np.array_equal(offset.owners, center.owners)
                    or not np.array_equal(offset.signs, center.signs) or not np.array_equal(offset.mask, center.mask)):
                raise ValueError('derivative stencil changed clipping recipes/exterior labels')
        return snapshots

    def trace(self, time):
        b = self.boundary
        exterior = self.extension.exterior
        return self.extension.extend(b['amplitude']*np.cos(
            b['grid_xy'][exterior]@b['wavevector']-b['frequency']*time))

    def context(self, panels):
        b = self.boundary
        return GraphBody(panels, b['eta'], b['length'], b['origin'], self.height, self.image, self.order)

    @staticmethod
    def normal_speed(snapshot):
        state, panels = snapshot.state, snapshot.panels
        omega = rotation_matrix(state.orientation)@state.angular_velocity
        speed = state.linear_velocity+np.cross(omega, panels.centroids-state.position)
        return np.sum(speed*panels.normals, axis=1)

    def evaluate(self, time, vector, step=None):
        started = perf_counter()
        step = self.step if step is None else step
        state = state_from_vector(vector)
        geometry_started = perf_counter()
        snapshots = self.stencil(state, step)
        geometry_seconds = perf_counter()-geometry_started
        center, plus, minus = snapshots
        context = self.context(center.panels)
        solution = context.solve(self.trace(time), self.normal_speed(center))
        velocity_started = perf_counter()
        velocity = body_velocity(context, solution.density, solution.spectrum, self.boundary['length'])
        velocity_seconds = perf_counter()-velocity_started
        actions = []
        for offset, offset_time in [(plus, time+step), (minus, time-step)]:
            action = frozen_field(offset.panels, self.boundary['eta'], self.boundary['length'],
                self.boundary['origin'], self.height, self.image, self.order, self.trace(offset_time), solution.density)
            action['residual'] = self.normal_speed(offset)-action['normal']
            actions.append(action)
        chart_velocity = (plus.panels.centroids-minus.panels.centroids)/(2*step)
        normal_rhs = (actions[0]['residual']-actions[1]['residual'])/(2*step)
        known_phi_t = (actions[0]['potential']-actions[1]['potential'])/(2*step) \
            -np.sum(chart_velocity*velocity, axis=1)
        coupled = coupled_pressure(context, self.body, state, normal_rhs, known_phi_t, velocity)
        # State angular velocity is WORLD; InitialState above owns BODY components.
        raw_q, omega = vector[3:7], vector[10:13]
        qdot = .5*np.r_[-omega@raw_q[1:], raw_q[0]*omega+np.cross(omega, raw_q[1:])]
        rate = np.r_[vector[7:10], qdot, coupled['acceleration']]
        return Evaluation(rate, coupled['pressure'], coupled['acceleration'],
            dict(time_seconds=float(time), total_seconds=perf_counter()-started,
                geometry_seconds=geometry_seconds, geometry=[s.diagnostics for s in snapshots],
                context_seconds=context.setup_seconds, current_solve=solution.diagnostics,
                context_breakdown=dict(source_operator_seconds=context.source_operator_seconds,
                    harmonic_maps_seconds=context.lift_map_seconds, source_trace_seconds=context.source_trace_seconds),
                body_velocity_seconds=velocity_seconds, offsets=[a['diagnostics'] for a in actions],
                coupled_solve=coupled['diagnostics'], acceleration_world=coupled['acceleration'].tolist(),
                pressure_functionals=coupled['pressure_functionals'].tolist(),
                derivative_step_seconds=step))

    def audit(self, time, vector, evaluation):
        """Fresh full current/offset contexts and seven modal solves are the oracle."""
        started = perf_counter()
        state, step = state_from_vector(vector), self.step
        center, plus, minus = self.stencil(state, step)
        context = self.context(center.panels)
        solution = context.solve(self.trace(time), self.normal_speed(center))
        velocity = body_velocity(context, solution.density, solution.spectrum, self.boundary['length'])
        offsets = []
        for snapshot, offset_time in [(plus, time+step), (minus, time-step)]:
            full = self.context(snapshot.panels)
            trace = self.trace(offset_time)
            right = full.right_hand_side(trace, self.normal_speed(snapshot))
            spectrum = full.spectrum(trace-full.trace@cp.asarray(solution.density))
            offsets.append((right-full.apply(solution.density), full.body_potential(solution.density, spectrum)))
        rhs = (offsets[0][0]-offsets[1][0])/(2*step)
        density_dot, _ = context.solve_density(rhs)
        spectrum_dot = -context.spectrum(context.trace@cp.asarray(density_dot))
        chart_velocity = (plus.panels.centroids-minus.panels.centroids)/(2*step)
        phi_t = (offsets[0][1]-offsets[1][1])/(2*step)+context.body_potential(density_dot, spectrum_dot) \
            -np.sum(chart_velocity*velocity, axis=1)
        known_pressure = -WATER_DENSITY*(phi_t+.5*np.sum(velocity**2, axis=1)+GRAVITY*center.panels.centroids[:, 2])
        relative = center.panels.centroids-state.position
        modes = np.column_stack((center.panels.normals, np.cross(relative, center.panels.normals)))
        weighted = center.panels.areas[:, None]*modes
        response = [context.solve(np.zeros(self.boundary['eta'].size), modes[:, axis]) for axis in range(6)]
        potential = np.column_stack([s.potential for s in response])
        fluid_mass = WATER_DENSITY*weighted.T@potential
        mass, external = body_dynamics(self.body, state)
        acceleration = np.linalg.solve(mass+fluid_mass, weighted.T@known_pressure+external)
        pressure = known_pressure-WATER_DENSITY*potential@acceleration
        return dict(total_seconds=perf_counter()-started,
            acceleration_difference=(evaluation.acceleration-acceleration).tolist(),
            pressure_max_difference_pa=float(np.max(abs(evaluation.pressure-pressure))),
            six_functional_difference=(weighted.T@(evaluation.pressure-pressure)).tolist(),
            fluid_mass_matrix_world=fluid_mass.tolist(),
            scope='Fresh full-context discrete pressure oracle at evolved pose; not independent physical validation')


def run(inputs, args, report):
    setup_started = perf_counter()
    boundary = boundary_data(inputs, args.grid, args.length)
    boundary_path = args.output.with_name(args.output.stem+'-boundary.npz')
    np.savez(boundary_path, **boundary)
    model = ForcedBody(inputs, boundary, args.derivative_step, args.body_refinement, args.order)
    report.update(setup=model.setup, setup_seconds=perf_counter()-setup_started,
                  graph_sha256=sha256(boundary['eta'].tobytes()).hexdigest(),
                  boundary_sha256=sha256(boundary_path.read_bytes()).hexdigest())
    initial = inputs.initial_state
    vector = np.r_[initial.position, initial.orientation, initial.linear_velocity,
                    rotation_matrix(initial.orientation)@initial.angular_velocity]
    times, states, dense = [0.], [vector.copy()], []
    rhs_started = perf_counter()

    def rhs(time, values):
        attempted = dict(time_seconds=float(time), state=values.tolist())
        report['rhs_attempts'] = report.get('rhs_attempts', 0)+1
        try:
            if len(report['rhs']) >= args.max_rhs or perf_counter()-rhs_started > args.budget_seconds:
                raise RuntimeError('bounded forced-motion RHS/time budget exhausted')
            result = model.evaluate(time, values)
        except BaseException:
            report['failed_attempted_rhs'] = attempted
            raise
        report['rhs'].append(result.diagnostics)
        print(json.dumps(dict(event='rhs', index=len(report['rhs']), time=float(time),
                              seconds=result.diagnostics['total_seconds'])), flush=True)
        return result.rate

    solver = None
    try:
        solver = RK23(rhs, 0., vector, args.duration, rtol=1e-5,
            atol=np.r_[np.full(3, 1e-7), np.full(4, 1e-8), np.full(3, 1e-7), np.full(3, 1e-8)],
            max_step=args.max_step, first_step=min(args.max_step, args.duration))
        while solver.status == 'running':
            message = solver.step()
            if solver.status == 'failed':
                raise RuntimeError(message)
            times.append(float(solver.t)); states.append(solver.y.copy())
            dense.append(solver.dense_output().Q.copy())
        report.update(integration_complete=True,
            final_attitude_degrees=attitude_degrees(states[-1][3:7]).tolist(),
            maximum_accepted_quaternion_norm_error=float(np.max(abs(np.linalg.norm(np.asarray(states)[:, 3:7], axis=1)-1))))
    finally:
        report.update(integration_seconds=perf_counter()-rhs_started,
            simulated_seconds=float(times[-1]), accepted_steps=len(times)-1,
            completed_rhs_evaluations=len(report['rhs']),
            rhs_evaluations=solver.nfev if solver is not None else report.get('rhs_attempts', 0))
        state_path = args.output.with_name(args.output.stem+'-states.npz')
        np.savez(state_path, times=np.asarray(times), states=np.asarray(states), dense_q=np.asarray(dense))
        report['states_sha256'] = sha256(state_path.read_bytes()).hexdigest()
    # Held-out dense-output state, not a starting fixture or one of its stencils.
    audit_time = args.duration/2
    index = min(np.searchsorted(times, audit_time, side='right')-1, len(times)-2)
    width = times[index+1]-times[index]
    powers = ((audit_time-times[index])/width)**np.arange(1, 4)
    audit_state = states[index]+width*dense[index]@powers
    audit_started = perf_counter()
    result = model.evaluate(audit_time, audit_state)
    oracle = model.audit(audit_time, audit_state, result)
    half = model.evaluate(audit_time, audit_state, step=args.derivative_step/2)
    report['heldout'] = dict(time_seconds=audit_time, state=audit_state.tolist(), full_context_oracle=oracle,
        derivative_half_step_acceleration_difference=(half.acceleration-result.acceleration).tolist(),
        derivative_half_step_pressure_max_difference_pa=float(np.max(abs(half.pressure-result.pressure))),
        total_seconds=perf_counter()-audit_started)
    report['complete'] = True


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--duration', type=float, default=.1)
    parser.add_argument('--max-step', type=float, default=.05)
    parser.add_argument('--derivative-step', type=float, default=.005)
    parser.add_argument('--grid', type=int, default=129)
    parser.add_argument('--length', type=float, default=64.)
    parser.add_argument('--body-refinement', type=int, choices=[1, 2], default=1)
    parser.add_argument('--order', type=int, choices=range(2, 9), default=4)
    parser.add_argument('--max-rhs', type=int, default=40)
    parser.add_argument('--budget-seconds', type=float, default=90.)
    parser.add_argument('--resource-context', default='Uncharacterized shared-machine load',
                        help='declared competing workload/idle context for timing interpretation')
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or args.grid < 17 or args.grid % 2 != 1 or args.max_rhs < 4
            or not all(np.isfinite(x) and x > 0 for x in [args.duration, args.max_step,
                args.derivative_step, args.length, args.budget_seconds])):
        parser.error('fresh outputs, odd grid>=17, >=4 RHS and finite positive controls required')
    archive_started = perf_counter()
    archive = args.output.with_suffix('.source'); archive.mkdir(parents=True)
    sources = [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'),
               *Path(__file__).parent.glob('*.py'), *Path(__file__).parent.glob('*.cu'),
               *Path(__file__).parent.glob('*.cuh')]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target); hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    physical = args.inputs.read_text(encoding='utf-8')
    inputs = SimulationInputs.from_json(physical)
    report = dict(complete=False, integration_complete=False, qualified=False, scope=__doc__,
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(), producer_sha256=hashes,
        resource_context=args.resource_context,
        controls={k:str(v) if isinstance(v, Path) else v for k,v in vars(args).items()}, rhs=[],
        limitations=['Fixed forced graph; no evolving free-wave dynamics or free-wave energy-conservation claim',
                     'Planar clipping differs from the Fourier interpolant near contact; mismatch/depth are reported',
                     'Stops on exterior-mask events; general contact-state transfer is not implemented',
                     'No open-water, spatial/physical, full-minute accuracy or <=20x speed qualification'])
    report['archive_seconds'] = perf_counter()-archive_started
    started = perf_counter()
    try:
        run(inputs, args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['numerical_run_seconds'] = perf_counter()-started
        report['total_seconds'] = report['archive_seconds']+report['numerical_run_seconds']
        report['timing_scope'] = 'Archive plus numerical setup/integration/audit; excludes Python imports and uv environment startup.'
        report['peak_process_rss_bytes'] = peak_rss_bytes()
        report['ending_cuda_pool_reserved_bytes'] = cp.get_default_memory_pool().total_bytes()
        report['memory_scope'] = 'Process peak RSS and ending CUDA pool; ending pool is not peak VRAM.'
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
