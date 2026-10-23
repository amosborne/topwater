"""Fixed-hull LINEAR time-domain radiation checkpoint, not a nonlinear solver.

Physical exterior eta/psi are the only fluid states, evolved on the GPU.
Geometry and fixed boundary maps are prepared on the CPU and uploaded once.
A fixed harmonic interior
fill supplies the fictitious plane trace for Fourier continuation; the current
body correction solves the prescribed total normal velocity. Odd grids avoid
ambiguous Nyquist continuation. The periodic lift/nonperiodic source field and
finite sponge require domain refinement before comparison is qualified.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.spatial import ConvexHull

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import WATER_DENSITY, GRAVITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from .boat_spatial_work import source_matrix, subdivide
from .panel_gradient_reference import checked_influence


def drive(time, omega, ramp_seconds):
    s = np.clip(time/ramp_seconds, 0., 1.)
    ramp = s**3*(10-15*s+6*s*s)
    rate = 30*s*s*(1-s)**2/ramp_seconds if time < ramp_seconds else 0.
    return ramp*np.cos(omega*time), rate*np.cos(omega*time)-omega*ramp*np.sin(omega*time)


def harmonic_maps(panels, size, spacing, origin):
    kx, ky = np.meshgrid(2*np.pi*np.fft.fftfreq(size, d=spacing),
                         2*np.pi*np.fft.fftfreq(size, d=spacing))
    wave = np.hypot(kx, ky)
    phase = ((panels.centroids[:, 0]-origin[0])[:, None]*kx.ravel()
             +(panels.centroids[:, 1]-origin[1])[:, None]*ky.ravel())
    potential = np.exp(1j*phase+panels.centroids[:, 2, None]*wave.ravel())/(size*size)
    normal = potential*(1j*panels.normals[:, 0, None]*kx.ravel()
        +1j*panels.normals[:, 1, None]*ky.ravel()+panels.normals[:, 2, None]*wave.ravel())
    return potential, normal, wave


class Radiation:
    """One fixed discrete plane/body BVP and its exact time derivative."""
    def __init__(self, inputs, size, length, refinement, sponge_strength):
        start = perf_counter()
        if size % 2 != 1:
            raise ValueError('odd grid required for unambiguous Fourier continuation')
        wet = hydrostatics(inputs.body, inputs.initial_state, density=WATER_DENSITY, gravity=GRAVITY)
        triangles = wet.wetted_triangles
        for _ in range(refinement):
            triangles = subdivide(triangles)
        self.panels = panels = TrianglePanels(triangles, -1)
        self.modes = np.column_stack((panels.normals,
            np.cross(panels.centroids-inputs.initial_state.position, panels.normals)))
        self.body_velocity = np.zeros((6, 2))
        self.body_velocity[2, 0], self.body_velocity[3, 1] = .01, .005
        self.velocity = self.modes@self.body_velocity
        self.amplitudes = np.array([.01, .005])
        self.size, self.length, self.spacing = size, length, length/size
        self.origin = origin = inputs.initial_state.position[:2]-length/2
        x, y = np.meshgrid(origin[0]+np.arange(size)*self.spacing,
                            origin[1]+np.arange(size)*self.spacing)
        xy = np.column_stack((x.ravel(), y.ravel()))
        footprint = ConvexHull(wet.waterline_segments[:, :, :2].reshape(-1, 2))
        if np.any(np.abs(footprint.points-inputs.initial_state.position[:2]) >= .5*length):
            raise ValueError('periodic box must strictly contain the actual waterplane')
        if not np.isclose(footprint.volume, wet.waterplane_area, rtol=1e-12, atol=1e-12):
            raise ValueError('convex physical waterplane required')
        distance = xy@footprint.equations[:, :2].T+footprint.equations[:, 2]
        # Only roundoff allowance: no diagnostic 5 mm exclusion in the fluid domain.
        inside = np.all(distance <= 64*np.finfo(float).eps*length, axis=1)
        self.exterior, self.hidden = np.flatnonzero(~inside), np.flatnonzero(inside)
        if not len(self.hidden):
            raise ValueError('grid must resolve the physical waterplane footprint')
        lookup = {int(index): row for row, index in enumerate(self.hidden)}
        boundary = sorted({neighbor for index in self.hidden for neighbor in self.neighbors(index)
                           if neighbor not in lookup})
        boundary_lookup = {index: row for row, index in enumerate(boundary)}
        matrix, rhs = np.eye(len(lookup))*4, np.zeros((len(lookup), len(boundary)))
        for index, row in lookup.items():
            for neighbor in self.neighbors(index):
                if neighbor in lookup:
                    matrix[row, lookup[neighbor]] -= 1
                else:
                    rhs[row, boundary_lookup[neighbor]] += 1
        self.fill_weights = np.linalg.solve(matrix, rhs)
        self.fill_boundary = np.array(boundary)
        potential, self.normal_map, self.wave = harmonic_maps(panels, size, self.spacing, origin)
        self.single = assemble_dirichlet_halfspace(panels).single_layer
        self.factor = lu_factor(source_matrix(panels))
        weighted_modes = panels.areas[:, None]*self.modes
        adjoint = lu_solve(self.factor, self.single.T@weighted_modes, trans=1)
        self.impulse_trace = WATER_DENSITY*(weighted_modes.T@potential-adjoint.T@self.normal_map)
        self.impulse_velocity = WATER_DENSITY*weighted_modes.T@self.single@lu_solve(self.factor, self.modes)
        del potential
        surface = np.column_stack((xy[self.exterior], np.zeros(len(self.exterior))))
        self.flux_map = np.empty((len(surface), len(panels)))
        for first in range(0, len(surface), 128):
            last = min(first+128, len(surface))
            self.flux_map[first:last] = 2*checked_influence(surface[first:last], panels).gradient_single_layer[:, :, 2]
        radius = np.max(np.abs(xy-inputs.initial_state.position[:2]), axis=1)/(length/2)
        position = np.clip((radius-.6)/.4, 0, 1)
        self.sponge = sponge_strength*position[self.exterior]**2
        self.cpu_static_setup_seconds = perf_counter()-start
        before = perf_counter()
        from .cuda_operator import cp
        from cupyx.scipy.linalg import lu_solve as device_lu_solve
        self.cp, self.device_lu_solve = cp, device_lu_solve
        self.factor = (cp.array(self.factor[0], order='F'), cp.asarray(self.factor[1]))
        for name in ('fill_weights', 'normal_map', 'velocity',
                     'impulse_trace', 'impulse_velocity', 'flux_map', 'wave', 'sponge', 'body_velocity'):
            setattr(self, name, cp.asarray(getattr(self, name)))
        self.device_exterior = cp.asarray(self.exterior)
        self.device_hidden = cp.asarray(self.hidden)
        self.device_fill_boundary = cp.asarray(self.fill_boundary)
        cp.cuda.get_current_stream().synchronize()
        self.device_setup_seconds = perf_counter()-before
        self.device_name = cp.cuda.runtime.getDeviceProperties(cp.cuda.Device().id)['name'].decode()
        self.evaluations = 0
        self.setup_seconds = perf_counter()-start

    def neighbors(self, index):
        row, col = divmod(int(index), self.size)
        return [row*self.size+(col-1) % self.size, row*self.size+(col+1) % self.size,
                ((row-1) % self.size)*self.size+col, ((row+1) % self.size)*self.size+col]

    def extend(self, exterior):
        cp = self.cp
        exterior = cp.asarray(exterior)
        full = cp.zeros((self.size*self.size, exterior.shape[1]), dtype=exterior.dtype)
        full[self.device_exterior] = exterior
        full[self.device_hidden] = self.fill_weights@full[self.device_fill_boundary]
        return full.reshape(self.size, self.size, -1)

    def rhs(self, time, state, omega, ramp):
        cp = self.cp
        state = cp.asarray(state)
        eta, psi = state
        spectrum = cp.fft.fft2(self.extend(psi), axes=(0, 1)).reshape(self.size*self.size, 2)
        speed, acceleration = drive(time, omega, ramp)
        density = self.device_lu_solve(self.factor, self.velocity*speed-(self.normal_map@spectrum).real)
        flat_flux = cp.fft.ifft2((spectrum.reshape(self.size, self.size, 2)
                                  *self.wave[:, :, None]), axes=(0, 1)).real.reshape(-1, 2)[self.device_exterior]
        flux = flat_flux+self.flux_map@density
        eta_rate = flux-self.sponge[:, None]*eta
        psi_rate = -GRAVITY*eta-self.sponge[:, None]*psi
        rate_spectrum = cp.fft.fft2(self.extend(psi_rate), axes=(0, 1)).reshape(-1, 2)
        force = -(self.impulse_trace@rate_spectrum).real-self.impulse_velocity@self.body_velocity*acceleration
        impulse = (self.impulse_trace@spectrum).real+self.impulse_velocity@self.body_velocity*speed
        energy = .5*WATER_DENSITY*self.spacing**2*cp.sum(GRAVITY*eta*eta+psi*flux, axis=0)
        energy += .5*cp.sum(impulse*self.body_velocity*speed, axis=0)
        nominal_sponge_power = -WATER_DENSITY*self.spacing**2*cp.sum(
            self.sponge[:, None]*(GRAVITY*eta*eta+psi*flux), axis=0)
        self.evaluations += 1
        return cp.stack((eta_rate, psi_rate)), force, energy, nominal_sponge_power


def fit(times, force, omega, start, amplitudes):
    use = times >= start
    t = times[use]
    design = np.column_stack((np.cos(omega*t), np.sin(omega*t), np.ones(len(t)), t-np.mean(t)))
    coefficients = np.linalg.lstsq(design, force[use].reshape(len(t), -1), rcond=None)[0]
    residual = force[use].reshape(len(t), -1)-design@coefficients
    cosine, sine = coefficients[:2].reshape(2, 6, 2)
    return dict(start_seconds=start, end_seconds=float(times[-1]),
                added_mass=(sine/(omega*amplitudes)).tolist(), damping=(-cosine/amplitudes).tolist(),
                force_fit_rms=np.sqrt(np.mean(residual*residual, axis=0)).reshape(6, 2).tolist())


def run(inputs, output, size, length, refinement, duration, step, ramp, sponge_strength, reference):
    started = perf_counter()
    system = Radiation(inputs, size, length, refinement, sponge_strength)
    cp = system.cp
    omega = 2*np.pi*.1
    if step*np.sqrt(GRAVITY*float(cp.max(system.wave).item())) > .5:
        raise ValueError('RK4 step too large for this resolved wave band')
    panel_hash = sha256(system.panels.triangles.tobytes()).hexdigest()
    references = [row for row in reference['rows'] if row['panel_sha256'] == panel_hash
                  and row['reference']['angular_frequency'] != 'infinite'
                  and np.isclose(row['reference']['angular_frequency'], omega, rtol=0, atol=1e-14)
                  and row['reference']['density'] == WATER_DENSITY
                  and row['reference']['gravity'] == GRAVITY
                  and row['reference']['water_depth'] == 'infinite'
                  and row['reference']['radiating_modes'] == ['Heave', 'Roll']
                  and row['reference']['influenced_modes'] == ['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw']]
    if len(references) != 1:
        raise ValueError('one exact-panel finite-frequency reference required')
    ref = references[0]['reference']
    count = round(duration/step)
    if not np.isclose(count*step, duration) or duration < ramp+30:
        raise ValueError('duration must be a multiple of step and leave >=30s after ramp')
    state = cp.zeros((2, len(system.exterior), 2))
    times = np.arange(count+1)*step
    forces = cp.empty((count+1, 6, 2))
    energies, sponge_power, peaks = [cp.empty((count+1, 2)) for _ in range(3)]
    cp.cuda.get_current_stream().synchronize()
    before = perf_counter()
    for index in range(count+1):
        time = index*step
        k1, force, energy, power = system.rhs(time, state, omega, ramp)
        forces[index], energies[index], sponge_power[index] = force, energy, power
        peaks[index] = cp.max(cp.abs(state[0]), axis=0)
        if index == count:
            break
        k2 = system.rhs(time+step/2, state+step*k1/2, omega, ramp)[0]
        k3 = system.rhs(time+step/2, state+step*k2/2, omega, ramp)[0]
        k4 = system.rhs(time+step, state+step*k3, omega, ramp)[0]
        state += step*(k1+2*k2+2*k3+k4)/6
        if not bool(cp.isfinite(state).all().item()):
            raise FloatingPointError('nonfinite radiation state')
        if index % max(1, count//10) == 0:
            print(json.dumps(dict(time=time, elapsed=perf_counter()-started,
                                  peak_eta=cp.asnumpy(peaks[index]).tolist())), flush=True)
    cp.cuda.get_current_stream().synchronize()
    evolution_seconds = perf_counter()-before
    before = perf_counter()
    forces, energies, sponge_power, peaks, state = [cp.asnumpy(value) for value in
                                                   (forces, energies, sponge_power, peaks, state)]
    output_transfer_seconds = perf_counter()-before
    before = perf_counter()
    fits = [fit(times, forces, omega, start, system.amplitudes) for start in [ramp+10, duration-20]]
    for result in fits:
        result['added_mass_error'] = (np.array(result['added_mass'])-ref['added_mass']).tolist()
        result['damping_error'] = (np.array(result['damping'])-ref['radiation_damping']).tolist()
    np.savez(output.with_suffix('.npz'), times=times, forces=forces, boundary_energy_proxy=energies,
             nominal_sponge_power=sponge_power, peak_eta=peaks, final_state=state,
             exterior=system.exterior, hidden=system.hidden)
    fit_and_output_seconds = perf_counter()-before
    return dict(complete=True, qualified=False,
        scope='LINEAR fixed-hull time-domain radiation; not nonlinear moving-boat throughput or open-water qualification',
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(), panel_sha256=panel_hash,
        grid_size=size, box_length_m=length, exterior_nodes=len(system.exterior), hidden_nodes=len(system.hidden),
        extension='Fixed 2D harmonic Dirichlet fill; no hidden eta/psi states', integrator='fixed-step RK4',
        time_step_seconds=step, duration_seconds=duration, ramp_seconds=ramp, frequency_hz=.1,
        sponge_strength_per_second=sponge_strength, sponge_start_fraction_of_half_box=.6,
        driven_modes=['Heave', 'Roll'], velocity_amplitudes=system.amplitudes.tolist(),
        force_rows=['Surge', 'Sway', 'Heave', 'Roll', 'Pitch', 'Yaw'],
        execution='GPU float64/complex128 evolution, CPU fixed geometry/operator preparation',
        device=system.device_name, cupy_version=cp.__version__,
        cpu_static_setup_seconds=system.cpu_static_setup_seconds,
        device_setup_seconds=system.device_setup_seconds,
        setup_seconds=system.setup_seconds, evolution_seconds=evolution_seconds,
        output_transfer_seconds=output_transfer_seconds, fit_and_output_seconds=fit_and_output_seconds,
        timing_scope='Run total includes static setup, lazy CUDA import/context/transfers, first-use compilation/plans, synchronized evolution, output download, fits and NPZ; source archiving and final JSON are outside this total.',
        total_seconds=perf_counter()-started, rhs_evaluations=system.evaluations,
        minimum_boundary_energy_proxy=energies.min(axis=0).tolist(), fits=fits,
        limits=['Periodic lift plus nonperiodic radiation; domain and sponge controls outstanding',
                'Sampled surface hole and raw body collocation require spatial refinement',
                'Nominal sponge power is a diagnostic; variable damping is not guaranteed passive',
                'Boundary energy is a mixed-domain collocation/quadrature proxy, not certified total physical energy',
                'No hydrostatic restoring, quadratic pressure, moving geometry or nonlinear graph in radiation comparison'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--grid', type=int, default=129)
    parser.add_argument('--length', type=float, default=128.)
    parser.add_argument('--refinement', type=int, default=1)
    parser.add_argument('--duration', type=float, default=60.)
    parser.add_argument('--step', type=float, default=.05)
    parser.add_argument('--ramp', type=float, default=20.)
    parser.add_argument('--sponge', type=float, default=.5)
    args = parser.parse_args()
    if (args.output.exists() or args.grid < 17 or args.grid % 2 != 1 or not 0 <= args.refinement <= 2
            or not np.isfinite([args.length, args.duration, args.step, args.ramp, args.sponge]).all()
            or min(args.length, args.duration, args.step, args.ramp) <= 0 or args.sponge < 0):
        parser.error('new output, odd grid>=17, refinement0..2 and finite positive length/time controls required')
    inputs = SimulationInputs.from_json(args.inputs.read_text(encoding='utf-8'))
    reference = json.loads(args.reference.read_text(encoding='utf-8'))
    if reference['input_sha256'] != sha256(inputs.to_json().encode()).hexdigest():
        parser.error('reference and physical inputs differ')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    sources = [*Path('src/topwater').glob('*.py'), Path('tests/laplace_reference.py'), Path(__file__),
               *[Path(__file__).with_name(name) for name in
                 ('boat_spatial_work.py', 'panel_gradient_reference.py', 'cuda_operator.py',
                  'cuda_energy.cu', 'cuda_pair.cuh')]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs, archive/'inputs.json')
    shutil.copyfile(args.reference, archive/'reference.json')
    result = run(inputs, args.output, args.grid, args.length, args.refinement, args.duration,
                 args.step, args.ramp, args.sponge, reference)
    result['producer_sha256'] = hashes
    result['reference_sha256'] = sha256(args.reference.read_bytes()).hexdigest()
    result['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
