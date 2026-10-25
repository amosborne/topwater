"""WEC adapter boundaries without a MATLAB or Capytaine simulation."""
from types import SimpleNamespace
import io
import subprocess
import shutil

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.io import savemat
from scipy.spatial.transform import Rotation

from topwater._wecsim_setup import floating_equilibrium
from topwater import wecsim_backend as backend
from .reference_boat import reference_body


def test_full_equilibrium_balances_force_and_moment_on_physical_boat():
    body = reference_body()
    state, h = floating_equilibrium(body)
    assert_allclose(h.force, [0., 0., body.mass*9.81], rtol=0, atol=1e-6)
    assert_allclose(h.moment, 0., rtol=0, atol=1e-6)
    assert_allclose(state.linear_velocity, 0.)
    assert_allclose(state.angular_velocity, 0.)
    assert h.rim_clearance > 0


def test_unqualified_surface_stops_before_any_physical_preparation(monkeypatch, tmp_path):
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


def test_modified_fork_source_is_rejected_even_at_pinned_head(monkeypatch, tmp_path):
    checkout=tmp_path/'fork'
    source=checkout/'source/objects/bodyClass.m'
    source.parent.mkdir(parents=True)
    source.write_text('% frozen producer\n')
    def git(*arguments):
        return subprocess.run([shutil.which('git'),'-C',str(checkout),*arguments],check=True,
                              capture_output=True,text=True).stdout.strip()
    git('init')
    git('add','source')
    git('-c','user.name=Test','-c','user.email=test@example.invalid','commit','-m','fixture')
    monkeypatch.setenv('TOPWATER_WECSIM_PATH',str(checkout))
    monkeypatch.setattr(backend,'SUPPORTED_FORK_COMMIT',git('rev-parse','HEAD'))
    assert backend._find_fork()==checkout.resolve()
    source.write_text('% changed producer\n')
    with pytest.raises(backend.BackendUnavailableError,match='source differs'):
        backend._find_fork()


def test_native_mean_pose_and_global_embedding_are_applied_once(monkeypatch, tmp_path):
    equilibrium = Rotation.from_euler('y', -.07)
    embedding = Rotation.from_euler('z', 1.2).as_matrix()
    center = np.array([0., 0., .2])
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
    setup = SimpleNamespace(time_s=np.array([0., .02]), mean_hull_vertices_m=vertices,
                            mean_cg_m=center, physical_faces=np.array([[0, 1, 2]]),
                            mass_matrix=np.eye(6), world_from_native=embedding,
                            equilibrium=SimpleNamespace(orientation=equilibrium.as_quat(scalar_first=True)))
    hydro = SimpleNamespace(K=np.eye(6), Ainf=np.eye(6), cb=np.zeros(3), volume=1.,
                            netcdf_path=tmp_path/'not-executed.nc')
    q = np.array([[0., 0., 0., 0., 0., 0.], [.3, -.2, .1, .04, -.03, .02]])

    def fake_native(workspace, fork):
        savemat(workspace/'motion.mat', dict(time=setup.time_s, q=q,
                velocity=np.zeros((2, 6)), acceleration=np.zeros((2, 6))))

    monkeypatch.setattr(backend, '_run_matlab', fake_native)
    motion = backend.integrate_motion(setup, hydro, np.zeros((2, 6)), tmp_path/'case', tmp_path)
    assert_allclose(motion.position_m, q[:, :3] @ embedding.T + center)
    expected = embedding @ Rotation.from_euler('xyz', q[:, 3:]).as_matrix() @ equilibrium.as_matrix()
    actual = Rotation.from_quat(motion.orientation_wxyz, scalar_first=True).as_matrix()
    assert_allclose(actual, expected, atol=1e-14)


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
