"""WEC-Sim adapter: native body integration and total-surface output.

Imports are lazy so the core package and renderer do not require MATLAB or
Capytaine. This adapter only attaches to the existing CAPYTAINE_REFERENCE
shared engine. Temporary BEMIO, Simulink and source-field products never become
public outputs. The caller owns atomic publication of the returned Result.
"""
from dataclasses import dataclass, replace
from contextlib import contextmanager
from importlib import import_module, metadata, resources, util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from time import perf_counter
import uuid

import numpy as np
from scipy.io import loadmat, savemat
from scipy.spatial.transform import Rotation

from ._wecsim_setup import prepare_case, TIME_STEP_S
from .constants import GRAVITY, WATER_DENSITY
from .hydrostatics import FloodingError, LossOfContactError, hydrostatics
from .inputs import InitialState


ENGINE_NAME = 'CAPYTAINE_REFERENCE'
SUPPORTED_FORK_COMMIT = '6b774228e1f6b31dee473e626204b0ef65e5ec38'
APPLICABILITY_WARNING_FRACTION = .10


class BackendUnavailableError(RuntimeError):
    """A required dependency or reviewed physical capability is unavailable."""


class EngineCleanupError(RuntimeError):
    """Shared-engine state is unverified; scratch must remain for recovery."""


class ScratchCleanupError(RuntimeError):
    """A complete run cannot publish while its scratch cleanup has failed."""


@dataclass(frozen=True)
class NativeMotion:
    time_s: np.ndarray
    q: np.ndarray
    velocity_native: np.ndarray
    acceleration_native: np.ndarray
    position_m: np.ndarray
    orientation_wxyz: np.ndarray
    excitation_native: np.ndarray


def _find_fork():
    git = shutil.which('git')
    if git is None:
        raise BackendUnavailableError('WEC-Sim: Git is required to verify the maintained fork')
    explicit = os.environ.get('TOPWATER_WECSIM_PATH')
    checkout = (Path(explicit).expanduser() if explicit else
                Path.home() / '.cache/topwater/wec-sim')
    if not (checkout / 'source/objects/bodyClass.m').is_file():
        raise BackendUnavailableError('WEC-Sim: maintained fork not found; set TOPWATER_WECSIM_PATH to its checkout')
    try:
        head = subprocess.run([git, '-C', str(checkout), 'rev-parse', 'HEAD'],
                              check=True, capture_output=True, text=True).stdout.strip()
        changed = subprocess.run([git, '-C', str(checkout), 'status', '--porcelain',
                                  '--untracked-files=all', '--', 'source'],
                                 check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as error:
        raise BackendUnavailableError('WEC-Sim: cannot establish the maintained fork revision') from error
    if head != SUPPORTED_FORK_COMMIT:
        raise BackendUnavailableError('WEC-Sim: fork revision differs from the supported reviewed version')
    if changed:
        raise BackendUnavailableError('WEC-Sim: maintained fork source differs from its reviewed revision')
    return checkout.resolve()


def preflight():
    """Check static capabilities before equilibrium or expensive producers.

    Discovering the shared engine and verifying licenses happens when executing
    the already scheduled solve, not as an import-time side effect.
    """
    for distribution, expected in (('capytaine', '3.0.0'), ('matlabengine', '25.2.2')):
        try:
            installed = metadata.version(distribution)
        except metadata.PackageNotFoundError as error:
            raise BackendUnavailableError(
                f'WEC-Sim: {distribution} is not installed; install the wecsim dependency group') from error
        if installed != expected:
            raise BackendUnavailableError(f'WEC-Sim: expected {distribution} {expected}, found {installed}')
    fork = _find_fork()
    if SUPPORTED_FORK_COMMIT == '9952634f6f2e2746a1626484e27c3003d1e1d5a9':
        raise BackendUnavailableError(
            'WEC-Sim: the held fork lacks corrected FullDir phases and directional weights; '
            'a separately validated fork revision is required')
    if util.find_spec('topwater._capytaine_surface') is None:
        raise BackendUnavailableError('WEC-Sim: required total-surface reconstruction is not available')
    producer = import_module('topwater._capytaine_surface')
    if not all(callable(getattr(producer, name, None)) for name in
               ('build_hydrodynamics', 'require_surface_readiness')):
        raise BackendUnavailableError('WEC-Sim: required total-surface reconstruction is not implemented')
    producer.require_surface_readiness()
    if not callable(getattr(producer, 'reconstruct_total_surface', None)):
        raise BackendUnavailableError('WEC-Sim: total-surface export is not implemented')
    engine_module = import_module('matlab.engine')
    if ENGINE_NAME not in engine_module.find_matlab():
        raise BackendUnavailableError(
            'WEC-Sim: existing CAPYTAINE_REFERENCE shared MATLAB engine is not visible; no new engine was started')
    engine = engine_module.connect_matlab(ENGINE_NAME)
    try:
        _check_engine(engine)
    finally:
        del engine
    return fork, producer


def solve(boat, sea, dwell_s, *, rest_s=10., ramp_s=5.):
    """Return complete standard motion AND total water-air geometry.

    No incident-only substitute is permitted. No public file is written here.
    """
    from .workflow_data import Result

    fork, producer = preflight()
    body = boat.to_body()
    waves = _waves(sea)
    started = _stage_started('equilibrium and incoming-train preparation')
    setup = prepare_case(body, waves, dwell_s, rest_s=rest_s, ramp_s=ramp_s)
    _stage_complete('equilibrium and incoming-train preparation', started)
    started = _stage_started('incident field startup assessment')
    setup.train.certify_incident_startup(setup)
    _stage_complete('incident field startup assessment', started)
    with _scratch_directory() as workspace:
        h = setup.hydrostatic
        # Extend genuine coefficient coverage for the complete selected packet.
        last_excitation_index = max(34, int(np.ceil(setup.train.omega[-1]*10)))
        excitation_omega = np.arange(1, last_excitation_index+1)/10
        radiation_omega = np.arange(1, max(120, last_excitation_index)+1)/10
        hydro = producer.build_hydrodynamics(
            setup.wet_triangles_m, setup.mean_cg_m, setup.mass_matrix, workspace / 'hydrodynamics',
            density=WATER_DENSITY, gravity=GRAVITY, depth=np.inf,
            hydrostatic_stiffness=h.stiffness,
            cb=h.center_of_buoyancy, volume=h.displaced_volume,
            excitation_omega=excitation_omega, radiation_omega=radiation_omega,
            green_cache_dir=Path.home()/'.cache/topwater/capytaine-3.0.0')
        motion = integrate_motion(setup, hydro, workspace / 'body', fork)
        started = _stage_started('calm-restoring diagnostics')
        warn_restoring_departures(setup, motion)
        _stage_complete('calm-restoring diagnostics', started)
        surface_time = np.linspace(0., motion.time_s[-1],
                                   max(2, int(round(motion.time_s[-1] / .08)) + 1))
        offset = float(setup.world_vertical_offset_m)
        computational_motion = replace(motion,
            position_m=motion.position_m-np.array([0., 0., offset]))
        offset_marker = workspace / 'surface-world-offset.json'
        offset_marker.write_text(json.dumps(dict(
            world_vertical_offset_m=offset, stored_surface_frame='computational_mean_zero',
            recovery='Translate recovered hull/body/surface world z by this offset once.'))+'\n',
            encoding='utf-8')
        try:
            vertices, faces, vertex_offsets, face_offsets = producer.reconstruct_total_surface(
                hydro, setup, computational_motion, surface_time)
        except BaseException as error:
            recovery = getattr(error, 'recovery_directory', None)
            if recovery is not None:
                try:
                    shutil.copy2(offset_marker, Path(recovery)/offset_marker.name)
                except OSError as marker_error:
                    error.add_note(f'World-offset recovery marker remains at {offset_marker}: {marker_error}')
            raise
        vertices = np.asarray(vertices).copy()
        vertices[:, 2] += offset
        return Result(
            time_s=motion.time_s, position_m=motion.position_m,
            orientation_wxyz=motion.orientation_wxyz,
            hull_vertices_m=np.asarray(boat.vertices_m), hull_faces=np.asarray(boat.faces),
            cg_m=np.asarray(boat.cg_m), surface_time_s=surface_time,
            surface_vertices_m=vertices, surface_faces=faces,
            surface_vertex_offsets=vertex_offsets, surface_face_offsets=face_offsets)


def _stage_started(name):
    started = perf_counter()
    print(f'STAGE {name} started', file=sys.stderr)
    return started


def _stage_complete(name, started):
    print(f'STAGE {name} complete elapsed_s={perf_counter()-started:.1f}', file=sys.stderr)


def _waves(sea):
    from .inputs import WaveComponents
    return WaveComponents(sea.amplitudes_m, sea.angular_frequencies_rad_s,
                          sea.directions_rad, sea.phases_rad)


def _validate_excitation_support(setup, hydro):
    """Validate coefficient/whole-packet coverage without synthesizing forces."""
    omega = np.asarray(hydro.excitation_omega, dtype=float)
    headings = np.asarray(hydro.excitation_headings, dtype=float)
    transfer = np.asarray(hydro.excitation)
    packet = np.asarray(setup.train.omega, dtype=float)
    if (omega.ndim != 1 or len(omega) < 2 or not np.isfinite(omega).all()
            or np.any(np.diff(omega) <= 0) or abs(omega[0]-.1) > 1e-12
            or packet.ndim != 1 or not len(packet) or not np.isfinite(packet).all()
            or np.any(packet < 0) or packet[-1] > omega[-1]):
        raise ValueError('excitation: retained whole packet is outside the supported transfer table')
    if (headings.ndim != 1 or len(headings) < 2 or not np.isfinite(headings).all()
            or np.any(np.diff(headings) <= 0) or abs(headings[0]) > 1e-12
            or abs(headings[-1]-2*np.pi) > 1e-12
            or transfer.shape != (len(omega), len(headings), 6)
            or not np.isfinite(transfer).all()
            or not np.allclose(transfer[:, 0], transfer[:, -1], rtol=0, atol=1e-10)):
        raise ValueError('excitation: coefficient grid or periodic heading seam is invalid')


def _prepare_native_inputs(setup, hydro, workspace, *, force_only=False):
    """Write physical spectral components and ordinary public hydroStruct coefficients."""
    from ._wecsim_spectral import prepare_spectral_inputs

    _validate_excitation_support(setup, hydro)
    spectral = prepare_spectral_inputs(setup, hydro)
    workspace = Path(workspace)
    workspace.mkdir(parents=True, exist_ok=False)
    library = resources.files('topwater').joinpath('data/wecsim')
    for name in ('wecSimInputFile.m', 'run_workflow_case.m', 'build_workflow_model.m',
                 'WorkflowEngineGuard.m', 'validate_workflow_hydrodynamics.m',
                 'prepare_workflow_hydrodynamics.m', 'native_excitation_startup.m',
                 'native_spectral_excitation.m'):
        workspace.joinpath(name).write_text(library.joinpath(name).read_text(encoding='utf-8'), encoding='utf-8')
    mean_vertices_cg = setup.mean_hull_vertices_m - setup.mean_cg_m
    _write_stl(workspace / 'boat_cg.stl', mean_vertices_cg, setup.physical_faces)
    spectral_file = workspace / 'spectral-wave.mat'
    savemat(spectral_file, spectral.wave, oned_as='column')
    # Physical component diagnostics only, never an applied force history.
    savemat(workspace / 'spectral-diagnostics.mat',
            dict(coefficients=spectral.coefficients,
                 original_omega_rad_s=spectral.original_omega_rad_s,
                 native_omega_rad_s=spectral.native_omega_rad_s,
                 omega_roundtrip_error_rad_s=spectral.omega_roundtrip_error_rad_s,
                 domega_rad_s=spectral.domega_rad_s, dtheta_rad=spectral.dtheta_rad))
    offset = float(getattr(setup, 'world_vertical_offset_m', 0.))
    native = dict(M=setup.mass_matrix, K=hydro.K, Ainf=hydro.Ainf,
                  cg=setup.mean_cg_m, cb=hydro.cb, volume=hydro.volume,
                  rho=WATER_DENSITY, g=GRAVITY, water_level_offset_m=offset,
                  initial_q=np.array([0., 0., -offset, 0., 0., 0.]),
                  time_step_s=TIME_STEP_S)
    carriers = setup.train.carriers
    savemat(workspace / 'prepared-case.mat',
            dict(n=native, time=setup.time_s[:, None],
                 spectral_file=str(spectral_file.resolve()), spectral_hydro=spectral.hydro,
                 incident_origin_elevation_m=setup.train.elevation([[0., 0.]], setup.time_s)[:, 0]-offset,
                 carrier_omega_rad_s=carriers.angular_frequencies,
                 carrier_amplitudes_m=carriers.amplitudes,
                 carrier_directions_rad=setup.train.theta,
                 boat_length_m=float(np.ptp(setup.body.mesh.vertices, axis=0).max()),
                 native_mat_path=str(Path(hydro.native_mat_path).resolve()),
                 netcdf_path=str(Path(hydro.netcdf_path).resolve()), force_only=bool(force_only),
                 model_name='topwater_' + uuid.uuid4().hex[:12]))
    return workspace


def integrate_motion(setup, hydro, workspace, fork):
    """Execute the unmodified maintained native WEC-Sim dynamics.

    WEC-Sim computes excitation from the exact physical spectral components. Only
    disposable model loggers are configured. Mean equilibrium
    orientation is already in the geometry and full mass tensor; it is applied
    once again only when converting native perturbations back to body pose.
    """
    started = _stage_started('native MATLAB preparation, integration and audits')
    workspace = _prepare_native_inputs(setup, hydro, workspace)
    _run_matlab(workspace, Path(fork))
    path = workspace / 'motion.mat'
    if not path.is_file():
        raise RuntimeError('WEC-Sim: native integration did not produce the required complete motion')
    raw = loadmat(path)
    stamp = np.asarray(raw['time']).reshape(-1)
    q = np.asarray(raw['q'])
    velocity, acceleration = np.asarray(raw['velocity']), np.asarray(raw['acceleration'])
    excitation = np.asarray(raw['excitation'])
    if (stamp.shape != setup.time_s.shape or not np.allclose(stamp, setup.time_s, rtol=0, atol=1e-9)
            or any(value.shape != (len(stamp), 6) or not np.isfinite(value).all()
                   for value in (q, velocity, acceleration, excitation))):
        raise RuntimeError('WEC-Sim: native output is incomplete or nonfinite')
    perturbation = Rotation.from_euler('xyz', q[:, 3:]).as_matrix()
    equilibrium = Rotation.from_quat(setup.equilibrium.orientation, scalar_first=True).as_matrix()
    orientation = Rotation.from_matrix(setup.world_from_native @ perturbation @ equilibrium).as_quat(scalar_first=True)
    position = q[:, :3] @ setup.world_from_native.T + setup.mean_cg_m
    position[:, 2] += float(getattr(setup, 'world_vertical_offset_m', 0.))
    _stage_complete('native MATLAB preparation, integration and audits', started)
    return NativeMotion(stamp, q, velocity, acceleration, position, orientation, excitation)


def _run_matlab(workspace, fork):
    engine_module = import_module('matlab.engine')
    if ENGINE_NAME not in engine_module.find_matlab():
        raise BackendUnavailableError(
            'WEC-Sim: existing CAPYTAINE_REFERENCE shared MATLAB engine is not visible; no new engine was started')
    engine = engine_module.connect_matlab(ENGINE_NAME)
    out, error = io.StringIO(), io.StringIO()
    original_directory = engine.pwd(nargout=1)
    original_path = engine.path(nargout=1)
    primary = None
    try:
        _check_engine(engine)
        # Register the unique owned model and snapshot session settings BEFORE
        # initializeWecSim or either template/model construction can fail.
        engine.eval(f"topwaterCaseDir={_quote(workspace)}; topwaterFork={_quote(fork)}; "
                    "addpath(topwaterCaseDir); "
                    "topwaterGuardInput=load(fullfile(topwaterCaseDir,'prepared-case.mat'),'model_name'); "
                    f"topwaterEngineGuard=WorkflowEngineGuard(topwaterGuardInput.model_name,{_quote(original_path)}); "
                    "topwaterEngineGuard.prepare(); "
                    f"run({_quote(workspace / 'run_workflow_case.m')});",
                    nargout=0, stdout=out, stderr=error)
    except Exception as failure:
        detail = (out.getvalue() + '\n' + error.getvalue()).strip()
        primary = RuntimeError(f'WEC-Sim native integration failed: {failure}\n{detail}')
        raise primary from failure
    except BaseException as failure:
        primary = failure
        raise
    finally:
        try:
            _restore_engine(engine, original_directory, original_path)
        except Exception as cleanup_failure:
            if primary is None:
                raise
            # CLI/logging display str(error), not its implicit exception context.
            # Keep the primary cause and buffered MATLAB output visible there,
            # while retaining the cleanup-error type that preserves scratch.
            message = f'{type(primary).__name__}: {primary}'
            transcript = (out.getvalue() + '\n' + error.getvalue()).strip()
            if transcript and transcript not in message:
                message += '\n' + transcript
            combined = EngineCleanupError(
                message + f'\nWEC-Sim shared-engine cleanup also failed: {cleanup_failure}')
            combined.primary_error = primary
            combined.cleanup_error = cleanup_failure
            raise combined from primary
        finally:
            # Never quit an engine shared with another client.
            del engine
    detail = (out.getvalue() + '\n' + error.getvalue()).strip()
    if detail:
        print(detail, file=sys.stderr)


def _restore_engine(engine, original_directory, original_path):
    """Restore settings, destroy the guard, then restore the original path.

    Its class folder must remain reachable through deletion; changing the path
    first invalidates MATLAB's handle-class destructor. Clearing only this class
    also permits the next solve to use a different temporary resource directory.
    """
    failures = []
    guard_exists = None
    try:
        guard_exists = bool(engine.eval("exist('topwaterEngineGuard','var') == 1", nargout=1))
    except Exception as failure:
        failures.append('guard presence: ' + str(failure))
    try:
        if guard_exists is not False:
            setting_failures = []
            settings_restored = False
            for attempt in range(2):
                try:
                    report = json.loads(engine.eval(
                        'jsonencode(topwaterEngineGuard.restore())', nargout=1))
                    if report.get('restored') is True:
                        settings_restored = True
                        break
                    details = report.get('failures', ['unverified non-path settings'])
                    if not isinstance(details, list):
                        details = [str(details)]
                    setting_failures.append(f'settings attempt {attempt+1}: ' + '; '.join(map(str, details)))
                except Exception as failure:
                    setting_failures.append(f'settings attempt {attempt+1}: ' + str(failure))
            if not settings_restored:
                failures.extend(setting_failures)
            # Keep these separate: a throwing destructor must not skip clear.
            for command, label in (
                    ('delete(topwaterEngineGuard)', 'guard deletion'),
                    ('clear topwaterEngineGuard topwaterGuardInput', 'guard release')):
                try:
                    engine.eval(command, nargout=0)
                except Exception as failure:
                    failures.append(label + ': ' + str(failure))
            try:
                guard_exists = bool(engine.eval("exist('topwaterEngineGuard','var') == 1", nargout=1))
                if guard_exists:
                    # An interrupted clear gets an independent release attempt
                    # while the class folder is still on the search path.
                    engine.eval('topwaterEngineGuard=[]; clear topwaterEngineGuard topwaterGuardInput', nargout=0)
                    guard_exists = bool(engine.eval("exist('topwaterEngineGuard','var') == 1", nargout=1))
                    if guard_exists:
                        failures.append('guard release: object remains in the shared workspace')
            except Exception as failure:
                guard_exists = None
                failures.append('guard release verification: ' + str(failure))
        if guard_exists is False:
            try:
                engine.eval('clear WorkflowEngineGuard', nargout=0)
            except Exception as failure:
                failures.append('guard class release: ' + str(failure))
    finally:
        # All cleanup paths attempt these last, including settings/JSON/delete
        # failures. Scratch is retained if any part remains unverifiable.
        try:
            engine.cd(original_directory, nargout=0)
            if engine.pwd(nargout=1) != original_directory:
                failures.append('working directory: original value was not restored')
        except Exception as failure:
            failures.append('working directory: ' + str(failure))
        try:
            _restore_search_path(engine, original_path)
        except Exception as failure:
            failures.append('search path: ' + str(failure))
    if failures:
        raise EngineCleanupError('WEC-Sim: shared-engine restoration failed: ' + '; '.join(failures))


def _restore_search_path(engine, original_path):
    failures = []
    for attempt in range(2):
        try:
            engine.path(original_path, nargout=0)
            if engine.path(nargout=1) == original_path:
                return
            failures.append(f'attempt {attempt+1}: original value was not restored')
        except Exception as failure:
            failures.append(f'attempt {attempt+1}: {failure}')
    raise EngineCleanupError('could not restore the original MATLAB search path after guard deletion: '
                             + '; '.join(failures))


@contextmanager
def _scratch_directory():
    """Retain recovery scratch and never replace a primary stage exception."""
    directory = Path(tempfile.mkdtemp(prefix='topwater-wecsim-'))
    retain = False
    primary = None
    try:
        yield directory
    except BaseException as failure:
        primary = failure
        recovery = getattr(failure, 'recovery_directory', None)
        if isinstance(failure, EngineCleanupError) or recovery is not None:
            retain = True
            reason = 'engine cleanup failed' if isinstance(failure, EngineCleanupError) else 'surface export failed'
            print(f'WEC-Sim: {reason}; temporary recovery directory retained: {directory}', file=sys.stderr)
            if recovery is not None:
                print(f'WEC-Sim: surface replay data: {recovery}', file=sys.stderr)
        raise
    finally:
        if not retain:
            try:
                shutil.rmtree(directory)
            except Exception as failure:
                message = f'WEC-Sim: temporary cleanup failed; recovery directory retained: {directory}: {failure}'
                print(message, file=sys.stderr)
                if primary is None:
                    raise ScratchCleanupError(message) from failure
                primary.add_note(message)


def _check_engine(engine):
    for name in ('Simulink', 'Simscape', 'SimMechanics'):
        if not engine.license('test', name, nargout=1):
            raise BackendUnavailableError(f'WEC-Sim: required MATLAB license {name} is unavailable')
    busy = engine.eval("any(cellfun(@(name) ~strcmp(get_param(name,'SimulationStatus'),'stopped'), "
                       "find_system('type','block_diagram','BlockDiagramType','model')))", nargout=1)
    if busy:
        raise BackendUnavailableError('WEC-Sim: shared MATLAB engine has a model that is not stopped')


def _quote(path):
    return "'" + str(path).replace('\\', '/').replace("'", "''") + "'"


class WarningSummary:
    """Print first events and compact final summaries; retain no public report."""
    def __init__(self, stream=None):
        self.stream = sys.stderr if stream is None else stream
        self.events = {}

    def observe(self, category, stamp, message, *, severity=None):
        stamp = float(stamp)
        if category not in self.events:
            self.events[category] = dict(first=stamp, last=stamp, count=0, peak=None)
            print(f'warning [{category}] at {stamp:g} s: {message}', file=self.stream)
        event = self.events[category]
        event['last'] = stamp
        event['count'] += 1
        if severity is not None:
            value = float(severity)
            if event['peak'] is None or value > event['peak']:
                event['peak'] = value

    def finish(self):
        for category, event in self.events.items():
            peak = '' if event['peak'] is None else f", peak departure {100*event['peak']:.3g}%"
            print(f"warning summary [{category}]: {event['count']} samples, "
                  f"first {event['first']:g} s, last {event['last']:g} s{peak}", file=self.stream)


def warn_restoring_departures(setup, motion, *, stream=None):
    """Compare exact CALM-WATER restoring with the fixed-mean linear wrench.

    This diagnostic is deliberately named for its flat-water reference. It is
    not an actual curved-water wetted-set or nonlinear dynamic accuracy check.
    The 10% threshold is an engineering reporting policy. A 0.1%-weight (or
    weight*length) denominator floor avoids percentages of numerical zero.
    """
    warning = WarningSummary(stream)
    weight = setup.body.mass * GRAVITY
    length = float(np.ptp(setup.body.mesh.vertices, axis=0).max())
    equilibrium = Rotation.from_quat(setup.equilibrium.orientation, scalar_first=True)
    for stamp, q in zip(motion.time_s, motion.q):
        perturbation = Rotation.from_euler('xyz', q[3:])
        state = InitialState(setup.mean_cg_m + q[:3],
                             (perturbation * equilibrium).as_quat(scalar_first=True))
        try:
            h = hydrostatics(setup.body, state, density=WATER_DENSITY, gravity=GRAVITY)
        except (FloodingError, LossOfContactError) as error:
            warning.observe('calm-restoring-unavailable', stamp,
                            f'{error}; mean-sea-reference geometry cannot assess restoring here')
            continue
        exact = np.r_[h.force + [0., 0., -weight], h.moment]
        linear = -setup.hydrostatic.stiffness @ q
        for name, indices, floor in (('force', slice(0, 3), .001*weight),
                                     ('moment', slice(3, 6), .001*weight*length)):
            difference = np.linalg.norm(exact[indices] - linear[indices])
            denominator = max(np.linalg.norm(exact[indices]), np.linalg.norm(linear[indices]), floor)
            fraction = float(difference / denominator)
            if fraction > APPLICABILITY_WARNING_FRACTION:
                warning.observe('calm-restoring-' + name, stamp,
                                f'exact mean-sea-reference calm-water {name} differs from linear restoring by {100*fraction:.3g}% '
                                '(10% reporting threshold; not a dynamic error estimate)', severity=fraction)
    warning.finish()


def _write_stl(path, vertices, faces):
    with Path(path).open('w', encoding='ascii') as stream:
        stream.write('solid topwater\n')
        for face in np.asarray(faces):
            triangle = vertices[face]
            normal = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
            normal /= np.linalg.norm(normal)
            stream.write('facet normal ' + ' '.join(map(str, normal)) + '\nouter loop\n')
            for point in triangle:
                stream.write('vertex ' + ' '.join(map(str, point)) + '\n')
            stream.write('endloop\nendfacet\n')
        stream.write('endsolid topwater\n')
