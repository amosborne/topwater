"""WEC adapter boundaries without a MATLAB or Capytaine simulation."""
from types import SimpleNamespace
import io
import json
from pathlib import Path
import subprocess
import shutil

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.io import savemat
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation

from topwater._wecsim_setup import IncomingTrain, floating_equilibrium, prepare_case
from topwater import wecsim_backend as backend
from topwater.hydrostatics import hydrostatics
from topwater.inputs import InitialState
from topwater.workflow_data import Boat, Sea
from .reference_boat import reference_body


def test_full_equilibrium_balances_force_and_moment_on_physical_boat():
    body = reference_body()
    state, h = floating_equilibrium(body)
    assert_allclose(h.force, [0., 0., body.mass*9.81], rtol=0, atol=1e-6)
    assert_allclose(h.moment, 0., rtol=0, atol=1e-6)
    assert_allclose(state.linear_velocity, 0.)
    assert_allclose(state.angular_velocity, 0.)
    assert h.rim_clearance > 0


def test_reference_boat_preserves_accepted_mean_hull_mass_and_restoring_identity():
    project = Path(__file__).resolve().parents[1]
    physical = json.loads((project/'tests/data/canonical-boat.json').read_text())
    sea = Sea.load(project/'examples/simple-boat/sea.npz')
    # Historical field evidence belongs to the canonical fixture. The public
    # example now uses the standard generator and has its own geometry.
    body = reference_body()
    setup = prepare_case(body, backend._waves(sea), 60.)
    expected = np.asarray(physical['minimum_wetted_triangles_world_m'])
    actual = setup.wet_triangles_m
    assert len(actual) == len(expected) == 17, 'selected-case field evidence requires the same 17 panels'
    distances = np.linalg.norm(actual.mean(axis=1)[:, None]-expected.mean(axis=1)[None], axis=2)
    rows, columns = linear_sum_assignment(distances)
    for i, j in zip(rows, columns):
        vertex_distance = np.linalg.norm(actual[i, :, None]-expected[j, None, :], axis=2)
        a, b = linear_sum_assignment(vertex_distance)
        assert_allclose(actual[i, a], expected[j, b], rtol=0, atol=1e-9)
        normal = np.cross(actual[i, 1]-actual[i, 0], actual[i, 2]-actual[i, 0])
        reference_normal = np.cross(expected[j, 1]-expected[j, 0], expected[j, 2]-expected[j, 0])
        assert normal @ reference_normal > 0, 'wet panel outward winding changed'
    assert_allclose(setup.mean_cg_m, physical['world_cg_m'], rtol=0, atol=1e-10)
    assert_allclose(setup.mass_matrix, physical['dry_mass_matrix_world'], rtol=0, atol=1e-8)
    frozen_state = InitialState(physical['world_cg_m'],
                                physical['minimum_hull']['equilibrium_quaternion_wxyz'])
    expected_hydrostatic = hydrostatics(body, frozen_state, density=1025., gravity=9.81)
    assert_allclose(setup.hydrostatic.stiffness, expected_hydrostatic.stiffness, rtol=0, atol=1e-6)


def test_unqualified_surface_stops_before_any_physical_preparation(monkeypatch, tmp_path):
    monkeypatch.setattr(backend, 'SUPPORTED_FORK_COMMIT', 'qualified-fixture')
    monkeypatch.setattr(backend.metadata, 'version',
                        lambda name: {'capytaine': '3.0.0', 'matlabengine': '25.2.2'}[name])
    monkeypatch.setattr(backend, '_find_fork', lambda: tmp_path)
    monkeypatch.setattr(backend.util, 'find_spec', lambda name: object())

    def unqualified():
        raise ValueError('total surface awaits physical pilot')

    producer = SimpleNamespace(build_hydrodynamics=lambda: None, require_surface_readiness=unqualified)
    monkeypatch.setattr(backend, 'import_module', lambda name: producer)
    monkeypatch.setattr(backend, 'prepare_case', lambda *args: pytest.fail('must not prepare physics'))
    with pytest.raises(ValueError, match='awaits physical pilot'):
        backend.preflight()


def test_held_fork_without_full_directional_bugfixes_fails_before_physics(monkeypatch, tmp_path):
    monkeypatch.setattr(backend, 'SUPPORTED_FORK_COMMIT', '9952634f6f2e2746a1626484e27c3003d1e1d5a9')
    monkeypatch.setattr(backend.metadata, 'version',
                        lambda name: {'capytaine': '3.0.0', 'matlabengine': '25.2.2'}[name])
    monkeypatch.setattr(backend, '_find_fork', lambda: tmp_path)
    monkeypatch.setattr(backend, 'prepare_case', lambda *args: pytest.fail('must not prepare physics'))
    with pytest.raises(backend.BackendUnavailableError, match='lacks corrected FullDir phases'):
        backend.preflight()


@pytest.mark.parametrize('explicit_path', [False, True])
def test_external_fork_discovery_rejects_modified_source(monkeypatch, tmp_path, explicit_path):
    if explicit_path:
        checkout = tmp_path/'configured fork'
        monkeypatch.setenv('TOPWATER_WECSIM_PATH', str(checkout))
    else:
        checkout = tmp_path/'.cache/topwater/wec-sim'
        monkeypatch.delenv('TOPWATER_WECSIM_PATH', raising=False)
        monkeypatch.setattr(Path, 'home', classmethod(lambda cls: tmp_path))
    source=checkout/'source/objects/bodyClass.m'
    source.parent.mkdir(parents=True)
    source.write_text('% frozen producer\n')
    def git(*arguments):
        return subprocess.run([shutil.which('git'),'-C',str(checkout),*arguments],check=True,
                              capture_output=True,text=True).stdout.strip()
    git('init')
    git('add','source')
    git('-c','user.name=Test','-c','user.email=test@example.invalid','commit','-m','fixture')
    monkeypatch.setattr(backend,'SUPPORTED_FORK_COMMIT',git('rev-parse','HEAD'))
    assert backend._find_fork()==checkout.resolve()
    source.write_text('% changed producer\n')
    with pytest.raises(backend.BackendUnavailableError,match='source differs'):
        backend._find_fork()


@pytest.mark.parametrize('water_level', [0., .013, -.017])
def test_native_mean_pose_and_global_embedding_are_applied_once(monkeypatch, tmp_path, water_level):
    equilibrium = Rotation.from_euler('y', -.07)
    embedding = Rotation.from_euler('z', 1.2).as_matrix()
    center = np.array([0., 0., .2])
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    setup = SimpleNamespace(time_s=np.array([0., .02]), mean_hull_vertices_m=vertices,
                            mean_cg_m=center, physical_faces=np.array([[0, 1, 2]]),
                            mass_matrix=np.eye(6), world_from_native=embedding,
                            equilibrium=SimpleNamespace(orientation=equilibrium.as_quat(scalar_first=True)),
                            world_vertical_offset_m=water_level)
    setup.body = SimpleNamespace(mesh=SimpleNamespace(vertices=vertices))
    setup.train = SimpleNamespace(
        directional_elevation_records=lambda times: (np.array([.2, .2+np.pi]),
                                                      np.c_[np.cos(times), .3*np.sin(times)]),
        carriers=SimpleNamespace(angular_frequencies=np.array([1.]), amplitudes=np.array([.2])),
        theta=np.array([.2]), omega=np.array([0., 1.4]))
    setup.train.elevation = lambda xy, times: np.full((len(times), len(xy)), water_level)
    hydro = SimpleNamespace(K=np.eye(6), Ainf=np.eye(6), cb=np.zeros(3), volume=1.,
                            netcdf_path=tmp_path/'not-executed.nc', native_mat_path=tmp_path/'native.mat')
    hydro.excitation_omega = np.array([.1, 3.4])
    hydro.excitation_headings = np.array([0., 2*np.pi])
    hydro.excitation = np.ones((2, 2, 6), dtype=complex)
    q = np.array([[0., 0., 0., 0., 0., 0.], [.3, -.2, .1, .04, -.03, .02]])
    q[:, 2] -= water_level

    def fake_native(workspace, fork):
        savemat(workspace/'motion.mat', dict(time=setup.time_s, q=q,
                velocity=np.zeros((2, 6)), acceleration=np.zeros((2, 6)), excitation=np.ones((2, 6))))

    monkeypatch.setattr(backend, '_run_matlab', fake_native)
    spectral = SimpleNamespace(
        wave=dict(frequencies=np.array([1., 2.]), directions=np.array([-90., 90.]),
                  spectrum=np.zeros(2), spread=np.zeros((2, 2)), phase=np.zeros((2, 2))),
        hydro=dict(simulation_parameters=dict(w=np.array([0., 1., 2.]),
                   direction=np.array([-180., 180.])),
                   hydro_coeffs=dict(excitation=dict(re=np.zeros((6, 2, 3)),
                   im=np.zeros((6, 2, 3))),mean_drift=np.zeros((6, 2, 3)))),
        coefficients=np.zeros((2, 2)), original_omega_rad_s=np.array([1., 2.]),
        native_omega_rad_s=np.array([1., 2.]), omega_roundtrip_error_rad_s=np.zeros(2),
        domega_rad_s=np.ones(2), dtheta_rad=np.ones(2))
    monkeypatch.setattr('topwater._wecsim_spectral.prepare_spectral_inputs', lambda *_: spectral)
    motion = backend.integrate_motion(setup, hydro, tmp_path/'case', tmp_path)
    prepared = backend.loadmat(tmp_path/'case/prepared-case.mat')
    assert 'force' not in prepared
    assert not (tmp_path/'case/compatibility-elevation.txt').exists()
    assert not list((tmp_path/'case').glob('wave-*.txt'))
    assert (tmp_path/'case/spectral-wave.mat').is_file()
    assert 'spectral_hydro' in prepared
    assert 'component_view' not in prepared
    assert_allclose(motion.excitation_native, 1.)
    assert_allclose(motion.position_m, q[:, :3] @ embedding.T + center+[0., 0., water_level])
    assert_allclose(motion.position_m[0], center, atol=1e-15)
    assert_allclose(prepared['n']['initial_q'][0, 0].reshape(-1), [0., 0., -water_level, 0., 0., 0.])
    assert_allclose(prepared['incident_origin_elevation_m'], 0.)
    expected = embedding @ Rotation.from_euler('xyz', q[:, 3:]).as_matrix() @ equilibrium.as_matrix()
    actual = Rotation.from_quat(motion.orientation_wxyz, scalar_first=True).as_matrix()
    assert_allclose(actual, expected, atol=1e-14)


def test_directional_records_reconstruct_origin_field_at_padded_absolute_times():
    train = IncomingTrain.__new__(IncomingTrain)
    train.omega = np.array([0., .7, 1.4])
    train.k = train.omega**2/9.81
    train.theta = np.array([.3, 1.7])
    train.carriers = SimpleNamespace(amplitudes=np.array([.4, .2]))
    train.C = np.array([[[0., .2+.3j, .1-.05j], [0., -.03+.02j, .01+.02j]],
                        [[0., -.1+.2j, .03+.01j], [0., .01-.04j, -.02+.01j]]])
    times = np.array([-20.02, -4., 0., 11., 60., 80.02])
    directions, records = train.directional_elevation_records(times)
    assert_allclose(directions, np.mod([.3, .3+np.pi, 1.7, 1.7+np.pi], 2*np.pi))
    assert_allclose(records.sum(axis=1), train.elevation([[0., 0.]], times)[:, 0],
                    rtol=0, atol=1e-14)
    xy = np.array([1.2, -.4])
    field = np.zeros(len(times))
    for j, theta in enumerate(directions):
        coefficients = train.C.reshape(-1, len(train.omega))[j]
        phase = train.k*(xy @ [np.cos(theta), np.sin(theta)])
        field += (np.exp(-1j*times[:, None]*train.omega+1j*phase) @ coefficients).real
    assert_allclose(field, train.elevation(xy, times)[:, 0], rtol=0, atol=1e-14)


def test_genuinely_calm_sea_has_a_physical_zero_record():
    train = IncomingTrain.__new__(IncomingTrain)
    train.carriers = SimpleNamespace(amplitudes=np.zeros(2))
    directions, records = train.directional_elevation_records(np.array([-1., 0., 1.]))
    assert_allclose(directions, [0.])
    assert_allclose(records, np.zeros((3, 1)))


@pytest.mark.parametrize('packet_max, accepted', [(3.39, True), (3.4, True), (3.41, False)])
def test_native_excitation_checks_whole_packet_not_only_carriers(packet_max, accepted):
    setup = SimpleNamespace(train=SimpleNamespace(omega=np.array([0., packet_max]),
                           carriers=SimpleNamespace(angular_frequencies=np.array([1.]))))
    hydro = SimpleNamespace(excitation_omega=np.array([.1, 3.4]),
                            excitation_headings=np.array([0., 2*np.pi]),
                            excitation=np.ones((2, 2, 6), dtype=complex))
    if accepted:
        backend._validate_excitation_support(setup, hydro)
    else:
        with pytest.raises(ValueError, match='whole packet'):
            backend._validate_excitation_support(setup, hydro)


def test_warning_aggregation_does_not_repeat_every_sample():
    stream = io.StringIO()
    warning = backend.WarningSummary(stream)
    warning.observe('restoring', 1., 'departure', severity=.11)
    warning.observe('restoring', 2., 'departure', severity=.21)
    warning.observe('restoring', 3., 'departure', severity=.12)
    warning.finish()
    lines = stream.getvalue().splitlines()
    assert len(lines) == 2
    assert 'at 1 s' in lines[0]
    assert '3 samples' in lines[1] and 'peak departure 21%' in lines[1]


def test_shared_engine_restoration_runs_when_initialization_throws(monkeypatch, tmp_path):
    class Engine:
        restored = False
        guard_deleted = False

        def pwd(self, **kwargs):
            return str(tmp_path)

        def cd(self, directory, **kwargs):
            assert self.guard_deleted
            assert directory == str(tmp_path)

        def path(self, *args, **kwargs):
            if args:
                assert self.guard_deleted, 'search path must remain intact until the guard is destroyed'
                assert args[0] == 'original path'
            return 'original path'

        def eval(self, command, **kwargs):
            if 'WorkflowEngineGuard(' in command:
                raise RuntimeError('initializeWecSim failed')
            if "exist('topwaterEngineGuard'" in command:
                return not self.guard_deleted
            if 'jsonencode(topwaterEngineGuard.restore())' in command:
                self.restored = True
                return '{"restored":true,"failures":[]}'
            if 'delete(topwaterEngineGuard)' in command:
                self.guard_deleted = True
                return None
            if command in ('clear topwaterEngineGuard topwaterGuardInput', 'clear WorkflowEngineGuard'):
                return None
            pytest.fail('unexpected engine evaluation: ' + command)

    engine = Engine()
    module = SimpleNamespace(find_matlab=lambda: [backend.ENGINE_NAME],
                             connect_matlab=lambda name: engine)
    monkeypatch.setattr(backend, 'import_module', lambda name: module)
    monkeypatch.setattr(backend, '_check_engine', lambda value: None)
    with pytest.raises(RuntimeError, match='initializeWecSim failed'):
        backend._run_matlab(tmp_path, tmp_path)
    assert engine.restored
    assert engine.guard_deleted


@pytest.mark.parametrize('native_fails', [False, True])
@pytest.mark.parametrize('caller_handling_error', [False, True])
def test_matlab_cleanup_failure_keeps_primary_transcript_and_recovery(
        monkeypatch, tmp_path, native_fails, caller_handling_error):
    from topwater.command_log import CommandLog

    native_error = RuntimeError('MATLAB_PRIMARY_MARKER')
    cleanup_error = backend.EngineCleanupError('ENGINE_CLEANUP_MARKER')
    directory = tmp_path/'scratch'
    directory.mkdir()
    real_rmtree = backend.shutil.rmtree
    monkeypatch.setattr(backend.tempfile, 'mkdtemp', lambda **kwargs: str(directory))

    class Engine:
        def pwd(self, **kwargs):
            return str(tmp_path)

        def path(self, **kwargs):
            return 'original MATLAB path'

        def eval(self, command, **kwargs):
            kwargs['stdout'].write('MATLAB_STDOUT_MARKER\n')
            kwargs['stderr'].write('MATLAB_STDERR_MARKER\n')
            if native_fails:
                raise native_error

    engine = Engine()
    module = SimpleNamespace(find_matlab=lambda: [backend.ENGINE_NAME],
                             connect_matlab=lambda name: engine)
    monkeypatch.setattr(backend, 'import_module', lambda name: module)
    monkeypatch.setattr(backend, '_check_engine', lambda value: None)

    def fail_cleanup(*args):
        raise cleanup_error

    monkeypatch.setattr(backend, '_restore_engine', fail_cleanup)
    log_path = tmp_path/'result.log'
    try:
        with pytest.raises(backend.EngineCleanupError) as caught:
            with CommandLog(log_path, 'solve'):
                with backend._scratch_directory():
                    if caller_handling_error:
                        try:
                            raise ValueError('UNRELATED_CALLER_MARKER')
                        except ValueError:
                            backend._run_matlab(directory, tmp_path)
                    else:
                        backend._run_matlab(directory, tmp_path)
        assert directory.is_dir(), 'Unverified engine cleanup must retain recovery scratch'
        text = str(caught.value)
        log = log_path.read_text(encoding='utf-8')
        assert 'UNRELATED_CALLER_MARKER' not in text and 'UNRELATED_CALLER_MARKER' not in log
        assert 'ENGINE_CLEANUP_MARKER' in text and 'ENGINE_CLEANUP_MARKER' in log
        assert 'END status=failed' in log and 'output_published=false' in log
        if native_fails:
            for marker in ('MATLAB_PRIMARY_MARKER', 'MATLAB_STDOUT_MARKER', 'MATLAB_STDERR_MARKER'):
                assert marker in text and marker in log
            assert caught.value.cleanup_error is cleanup_error
            assert caught.value.primary_error.__cause__ is native_error
            assert caught.value.__cause__ is caught.value.primary_error
        else:
            assert caught.value is cleanup_error
    finally:
        real_rmtree(directory)


def test_incomplete_engine_cleanup_retries_and_fails_explicitly():
    class Engine:
        restores = 0
        guard_present = True
        events = []

        def eval(self, command, **kwargs):
            self.events.append(command)
            if command.startswith('exist('):
                return self.guard_present
            if command.startswith('jsonencode('):
                self.restores += 1
                return '{"restored":false,"failures":["file generation folders: unavailable"]}'
            if command.startswith('delete('):
                self.guard_present = False
                return None
            if command.startswith('clear '):
                return None
            pytest.fail('unexpected cleanup evaluation')

        def cd(self, directory, **kwargs):
            self.events.append('cd')

        def pwd(self, **kwargs):
            return 'old directory'

        def path(self, *args, **kwargs):
            self.events.append('path')
            return 'old path'

    engine = Engine()
    with pytest.raises(backend.EngineCleanupError, match='file generation folders'):
        backend._restore_engine(engine, 'old directory', 'old path')
    assert engine.restores == 2
    assert 'delete(topwaterEngineGuard)' in engine.events
    assert 'clear topwaterEngineGuard topwaterGuardInput' in engine.events
    assert engine.events[-1] == 'path'


def test_final_path_failure_is_an_engine_cleanup_error_after_guard_deletion():
    class Engine:
        deleted = False

        def eval(self, command, **kwargs):
            if command.startswith('exist('):
                return not self.deleted
            if command.startswith('jsonencode('):
                return '{"restored":true,"path_deferred":true,"failures":[]}'
            if command.startswith('delete('):
                self.deleted = True
                return None
            if command.startswith('clear '):
                return None
            pytest.fail('unexpected engine evaluation')

        def cd(self, directory, **kwargs):
            pass

        def pwd(self, **kwargs):
            return 'old directory'

        def path(self, *args, **kwargs):
            assert self.deleted
            raise RuntimeError('path restoration unavailable')

    engine = Engine()
    with pytest.raises(backend.EngineCleanupError, match='search path after guard deletion'):
        backend._restore_engine(engine, 'old directory', 'old path')
    assert engine.deleted


@pytest.mark.parametrize('failure_kind', ['json', 'settings-call', 'delete'])
def test_cleanup_exceptions_still_attempt_release_and_final_path(failure_kind):
    class Engine:
        present = True
        events = []

        def eval(self, command, **kwargs):
            self.events.append(command)
            if command.startswith('exist('):
                return self.present
            if command.startswith('jsonencode('):
                if failure_kind == 'json':
                    return 'not JSON'
                if failure_kind == 'settings-call':
                    raise RuntimeError('settings transport failed')
                return '{"restored":true,"failures":[]}'
            if command.startswith('delete('):
                if failure_kind == 'delete':
                    raise RuntimeError('destructor failed')
                self.present = False
                return None
            if command.startswith('clear topwaterEngineGuard'):
                self.present = False
                return None
            if command == 'clear WorkflowEngineGuard':
                return None
            pytest.fail('unexpected cleanup evaluation: ' + command)

        def cd(self, directory, **kwargs):
            self.events.append('cd')

        def pwd(self, **kwargs):
            return 'old directory'

        def path(self, *args, **kwargs):
            self.events.append('path')
            raise RuntimeError('final path failed')

    engine = Engine()
    with pytest.raises(backend.EngineCleanupError) as caught:
        backend._restore_engine(engine, 'old directory', 'old path')
    message = str(caught.value)
    assert ('guard deletion: destructor failed' if failure_kind == 'delete' else 'settings attempt 1:') in message
    assert 'final path failed' in message
    assert 'delete(topwaterEngineGuard)' in engine.events
    assert 'clear topwaterEngineGuard topwaterGuardInput' in engine.events
    assert not engine.present
    assert engine.events[-1] == 'path'


@pytest.mark.parametrize('unreconciled', [False, True])
def test_scratch_is_removed_only_after_verified_engine_cleanup(unreconciled):
    directory = None
    expected = backend.EngineCleanupError if unreconciled else ValueError
    try:
        with pytest.raises(expected):
            with backend._scratch_directory() as directory:
                raise expected('failure after a workflow stage')
        assert directory.exists() == unreconciled
    finally:
        if directory is not None and directory.exists():
            backend.shutil.rmtree(directory)


def test_scratch_cleanup_failure_preserves_primary_exception(monkeypatch, tmp_path, capsys):
    directory = tmp_path/'scratch'
    directory.mkdir()
    real_rmtree = backend.shutil.rmtree
    primary = ValueError('surface geometry failed at frame 217')
    monkeypatch.setattr(backend.tempfile, 'mkdtemp', lambda **kwargs: str(directory))

    def locked(path):
        raise PermissionError('KH.dat is locked')

    monkeypatch.setattr(backend.shutil, 'rmtree', locked)
    try:
        with pytest.raises(ValueError) as caught:
            with backend._scratch_directory():
                raise primary
        assert caught.value is primary
        assert 'KH.dat is locked' in capsys.readouterr().err
        assert any('temporary cleanup failed' in note for note in primary.__notes__)
        assert directory.exists()
    finally:
        real_rmtree(directory)


def test_scratch_cleanup_failure_without_primary_is_fatal(monkeypatch, tmp_path):
    directory = tmp_path/'scratch'
    directory.mkdir()
    real_rmtree = backend.shutil.rmtree
    monkeypatch.setattr(backend.tempfile, 'mkdtemp', lambda **kwargs: str(directory))

    def locked(path):
        raise PermissionError('KH.dat is locked')

    monkeypatch.setattr(backend.shutil, 'rmtree', locked)
    try:
        with pytest.raises(backend.ScratchCleanupError, match='KH.dat is locked'):
            with backend._scratch_directory():
                pass
        assert directory.exists()
    finally:
        real_rmtree(directory)


def test_surface_replay_exception_retains_own_scratch(monkeypatch, tmp_path):
    directory = tmp_path/'scratch'
    directory.mkdir()
    real_rmtree = backend.shutil.rmtree

    class SurfaceExportError(RuntimeError):
        recovery_directory = directory/'surface-recovery'

    primary = SurfaceExportError('invalid waterline')
    monkeypatch.setattr(backend.tempfile, 'mkdtemp', lambda **kwargs: str(directory))
    monkeypatch.setattr(backend.shutil, 'rmtree', lambda path: pytest.fail('replay directory must remain'))
    try:
        with pytest.raises(SurfaceExportError) as caught:
            with backend._scratch_directory():
                raise primary
        assert caught.value is primary
        assert directory.exists()
    finally:
        real_rmtree(directory)
