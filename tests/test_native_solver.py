"""Motion, causal memory, and external feedback contracts without MATLAB."""
from pathlib import Path
from types import SimpleNamespace
import hashlib
import json

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.spatial.transform import Rotation

from topwater.inputs import WaveComponents
from topwater.solve import LoadEvaluation
from topwater.solve.native.dynamics import effective_inertia
from topwater.solve.native.dynamics import integrate_motion
from topwater.solve.native.excitation import Excitation
from topwater.solve.native.radiation import history_force
from topwater.solve.native.radiation import radiation_kernel
from topwater.solve.native.radiation import stage_weights
from topwater.workflow_data import Boat
from topwater.workflow_data import Sea
from topwater.solve.preparation import prepare_case


@pytest.fixture(scope='module')
def reference():
    root = Path(__file__).resolve().parents[1]
    with np.load(root/'tests/data/simple_boat_wecsim.npz', allow_pickle=False) as archive:
        data = {name: archive[name] for name in archive.files}
    provenance = json.loads(str(data['metadata_json']))
    for name, digest in provenance['inputs'].items():
        assert hashlib.sha256((root/'examples/simple-boat-wecsim'/name).read_bytes()).hexdigest() == digest
    boat = Boat.load(root/'examples/simple-boat-wecsim/boat.npz')
    sea = Sea.load(root/'examples/simple-boat-wecsim/sea.npz')
    waves = WaveComponents(sea.amplitudes_m, sea.angular_frequencies_rad_s,
                           sea.directions_rad, sea.phases_rad)
    setup = prepare_case(boat.to_body(), waves, 100., rest_s=10., ramp_s=5.)
    hydro = SimpleNamespace(**{name: data[name] for name in
        ('M', 'K', 'Ainf', 'radiation_omega', 'radiation_damping', 'excitation_omega',
         'excitation_headings', 'excitation', 'density', 'gravity')})
    return setup, hydro, data


def reference_excitation(setup, hydro):
    from scipy.interpolate import RegularGridInterpolator
    from topwater.solve.wecsim.spectral import _coherent_waves
    _, omega, _, directions, _, waves = _coherent_waves(setup.train)
    source = hydro.excitation_omega
    query = np.empty((len(omega), len(directions), 2))
    query[:, :, 0] = np.maximum(omega, source[0])[:, None]
    query[:, :, 1] = directions[None, :]
    transfer = RegularGridInterpolator((source, hydro.excitation_headings), hydro.excitation)(query)
    low = omega < source[0]
    limit = hydro.K[:, 2]
    transfer[low] = limit+(omega[low, None, None]/source[0])**2*(transfer[low]-limit)
    return Excitation(omega, np.sum(transfer*waves[:, :, None], axis=1))


def test_reference_force_operators(reference):
    setup, hydro, data = reference
    excitation = reference_excitation(setup, hydro)
    actual = np.array([excitation.at(t) for t in data['time_s']])
    assert_allclose(actual, data['excitation_native'], rtol=1e-11, atol=1e-8)
    # Frozen WEC-Sim operators intentionally retain their old 2401-node rule.
    # Reproduce it independently; analytical native integration is tested against
    # adaptive quadrature in test_analytical_radiation.py.
    from scipy.integrate import trapezoid
    grid=np.linspace(hydro.radiation_omega[0],hydro.radiation_omega[-1],2401)
    time=np.arange(len(data['radiation_kernel']))*.01
    kernel=np.empty_like(data['radiation_kernel'])
    for i in range(6):
        for j in range(6):
            values=np.interp(grid,hydro.radiation_omega,hydro.radiation_damping[:,i,j]/hydro.radiation_omega)*grid
            kernel[:,i,j]=2/np.pi*trapezoid(np.cos(time[:,None]*grid)*values,grid,axis=1)
    assert_allclose(kernel, data['radiation_kernel'], rtol=1e-11, atol=1e-8)
    solid, residual = effective_inertia(hydro.M, hydro.Ainf)
    assert_allclose(solid, data['solid'], rtol=0, atol=1e-10)
    assert_allclose(residual, data['residual_added_mass'], rtol=0, atol=1e-10)


def test_full_example_converges(reference):
    setup, hydro, data = reference
    excitation = reference_excitation(setup, hydro)
    coarse_kernel = radiation_kernel(hydro.radiation_omega, hydro.radiation_damping, duration_s=10.)[1]
    coarse = integrate_motion(setup, hydro, excitation=excitation, kernel=coarse_kernel)
    fine_kernel = radiation_kernel(hydro.radiation_omega, hydro.radiation_damping, dt_s=.005, duration_s=10.)[1]
    fine = integrate_motion(setup, hydro, dt_s=.005, excitation=excitation, kernel=fine_kernel)
    assert_allclose(coarse.time_s, data['time_s'], rtol=0, atol=1e-9)
    assert np.max(np.linalg.norm(coarse.position_m-fine.position_m[::2], axis=1)) < .001
    error = (Rotation.from_quat(coarse.orientation_wxyz, scalar_first=True)
             *Rotation.from_quat(fine.orientation_wxyz[::2], scalar_first=True).inv()).magnitude()
    assert error.max() < np.deg2rad(.01)
    assert np.max(abs(coarse.force_components['balance_residual'])) < 1e-8
    # Like-for-like regression: this test intentionally retains historical
    # interpolated excitation and 10 s memory. Default native uses a different
    # excitation/history model, so its difference is reported separately.
    assert np.max(np.linalg.norm(coarse.position_m-data['position_m'], axis=1)) < .01
    reference_error = (Rotation.from_quat(coarse.orientation_wxyz, scalar_first=True)
                       *Rotation.from_quat(data['orientation_wxyz'], scalar_first=True).inv()).magnitude()
    assert reference_error.max() < np.deg2rad(.1)


class NoLoad:
    def initial_state(self, boat):
        return np.zeros(0)

    def evaluate(self, boat, state):
        return LoadEvaluation()


class ConstantExcitation:
    def __init__(self, value=None):
        self.value = np.zeros(6) if value is None else np.asarray(value)

    def at(self, stamp):
        return self.value


def plant(end=.2, *, heading=0., equilibrium=None, mass=None, added=None, stiffness=None):
    setup = SimpleNamespace(time_s=np.arange(round(end/.01)+1)*.01,
        world_vertical_offset_m=0., mean_cg_m=np.array([2., -3., 1.]),
        world_from_native=Rotation.from_euler('z', heading).as_matrix(),
        equilibrium=SimpleNamespace(orientation=np.array([1., 0., 0., 0.])
            if equilibrium is None else equilibrium))
    hydro = SimpleNamespace(M=np.eye(6) if mass is None else mass,
        Ainf=np.zeros((6, 6)) if added is None else added,
        K=np.zeros((6, 6)) if stiffness is None else stiffness)
    return setup, hydro


def run_plant(setup, hydro, influence=None, *, excitation=None, dt=.01):
    return integrate_motion(setup, hydro, influence=influence, dt_s=dt,
                            excitation=ConstantExcitation() if excitation is None else excitation,
                            kernel=np.zeros((11, 6, 6)))


def test_calm_equilibrium_and_zero_external_load_are_identical():
    setup, hydro = plant(equilibrium=Rotation.from_euler('xy', [.03, -.02]).as_quat(scalar_first=True))
    empty = run_plant(setup, hydro)
    model = run_plant(setup, hydro, NoLoad())
    assert_allclose(empty.position_m, np.broadcast_to(setup.mean_cg_m, empty.position_m.shape), atol=1e-15)
    np.testing.assert_array_equal(model.position_m, empty.position_m)
    np.testing.assert_array_equal(model.orientation_wxyz, empty.orientation_wxyz)
    assert not empty.velocity_native.any()


def test_torque_is_transformed_and_uses_full_coupled_added_mass():
    added = np.eye(6)*.1
    added[0, 3] = added[3, 0] = .04
    setup, hydro = plant(heading=np.pi/2, added=added)
    class Torque(NoLoad):
        def evaluate(self, boat, state):
            return LoadEvaluation(torque_world_Nm=np.array([0., 2., 0.]))
    motion = run_plant(setup, hydro, Torque())
    expected = np.linalg.solve(hydro.M+added, np.array([0., 0., 0., 2., 0., 0.]))
    assert_allclose(motion.acceleration_native[0], expected, atol=1e-14)
    assert abs(motion.acceleration_native[0, 0]) > 0
    assert_allclose(motion.force_components['additional'][:, 3:],
                    np.broadcast_to([2., 0., 0.], (len(motion.time_s), 3)), atol=1e-14)


def test_auxiliary_state_and_boat_are_integrated_together():
    class Device:
        def initial_state(self, boat):
            return np.array([1., 0.])
        def evaluate(self, boat, state):
            return LoadEvaluation(torque_world_Nm=np.array([state[0], 0., 0.]),
                                  state_derivative=np.array([state[1], -state[0]]))
    setup, hydro = plant(end=1.)
    motion = run_plant(setup, hydro, Device())
    assert_allclose(motion.influence_state[:, 0], np.cos(motion.time_s), atol=1e-10)
    assert_allclose(motion.q[:, 3], 1-np.cos(motion.time_s), atol=1e-10)


def test_sampled_roll_feedback_updates_only_at_accepted_steps():
    class Feedback(NoLoad):
        def __init__(self):
            self.updates, self.evaluations, self.power = [], [], []
            self.command = 0.
        def on_step(self, boat, state):
            self.updates.append(boat.time_s)
            roll = Rotation.from_quat(boat.orientation_wxyz, scalar_first=True).as_euler('xyz')[0]
            self.command = -10*roll
        def evaluate(self, boat, state):
            self.evaluations.append(len(self.updates))
            damping = -2*boat.angular_velocity_rad_s
            self.power.append(damping@boat.angular_velocity_rad_s)
            return LoadEvaluation(torque_world_Nm=damping+np.array([self.command, 0., 0.]))
    setup, hydro = plant(equilibrium=Rotation.from_euler('x', .1).as_quat(scalar_first=True))
    model = Feedback()
    motion = run_plant(setup, hydro, model)
    assert_allclose(model.updates, motion.time_s, atol=0)
    assert len(model.evaluations) == 4*(len(motion.time_s)-1)+1
    assert model.evaluations[:4] == [1, 1, 1, 1]
    assert max(model.power) <= 0 and min(model.power) < 0
    roll = Rotation.from_quat(motion.orientation_wxyz[-1], scalar_first=True).as_euler('xyz')[0]
    assert 0 < roll < .1


@pytest.mark.parametrize('bad', ['nan', 'shape', 'derivative', 'exception'])
def test_external_failures_name_model_and_time(bad):
    class Broken(NoLoad):
        def evaluate(self, boat, state):
            if bad == 'exception':
                raise RuntimeError('controller failure')
            if bad == 'derivative':
                return LoadEvaluation(state_derivative=np.ones(1))
            return LoadEvaluation(force_world_N=np.full(3, np.nan) if bad == 'nan' else np.zeros(2))
    setup, hydro = plant()
    with pytest.raises((ValueError, RuntimeError), match='Broken.*t=0'):
        run_plant(setup, hydro, Broken())


def test_external_snapshots_cannot_change_solver_state():
    class Inspector(NoLoad):
        def evaluate(self, boat, state):
            for array in (boat.position_m, boat.orientation_wxyz, boat.velocity_m_s,
                          boat.angular_velocity_rad_s, state):
                assert not array.flags.writeable
            # Even a deliberate override modifies only the external copy.
            boat.position_m.flags.writeable = True
            boat.position_m[:] = 900
            return LoadEvaluation()
    setup, hydro = plant()
    motion = run_plant(setup, hydro, Inspector())
    assert_allclose(motion.position_m[-1], setup.mean_cg_m)


@pytest.mark.parametrize('initial', [np.array([np.nan]), np.array([1j]), np.ones((1, 1))])
def test_invalid_device_initial_state_has_context(initial):
    class BadInitial(NoLoad):
        def initial_state(self, boat):
            return initial
    setup, hydro = plant()
    with pytest.raises(ValueError, match='BadInitial.*t=0'):
        run_plant(setup, hydro, BadInitial())


@pytest.mark.parametrize('fraction', [0., .5, 1.])
def test_causal_memory_partial_interval_and_cutoff(fraction):
    dt, memory = .1, 1.
    kernel = np.broadcast_to(np.eye(6), (11, 6, 6)).copy()
    weights, current = stage_weights(kernel, dt, fraction)
    # v=t with known history, so the exact integral is memory*(t-memory/2).
    history = np.arange(21)[:, None]*dt*np.ones((1, 6))
    stamp = 2+fraction*dt
    actual = history_force(history, 20, weights)+current@np.full(6, stamp)
    assert_allclose(actual, np.full(6, memory*(stamp-memory/2)), atol=1e-14)
    # At startup the only nonzero interval is [0,fraction*dt].
    actual = history_force(np.zeros((1, 6)), 0, weights)+current@np.full(6, fraction*dt)
    assert_allclose(actual, np.full(6, .5*(fraction*dt)**2), atol=1e-14)


def test_regular_wave_matches_analytic_oscillator():
    K = np.diag([0., 0., 4., 0., 0., 0.])
    setup, hydro = plant(end=2., stiffness=K)
    class Wave:
        def at(self, stamp):
            return np.array([0., 0., np.cos(stamp), 0., 0., 0.])
    motion = run_plant(setup, hydro, excitation=Wave())
    assert_allclose(motion.q[:, 2], (np.cos(motion.time_s)-np.cos(2*motion.time_s))/3, atol=2e-9)



def test_full_history_uses_complete_case_kernel_beyond_ten_seconds():
    setup, hydro = plant(end=12.)
    hydro.radiation_omega = np.array([.1, .5, 1.])
    hydro.radiation_damping = np.broadcast_to(np.eye(6), (3, 6, 6))
    excitation = ConstantExcitation([1., 0., 0., 0., 0., 0.])
    full = integrate_motion(setup, hydro, excitation=excitation, dt_s=.1)
    kernel = radiation_kernel(hydro.radiation_omega, hydro.radiation_damping, duration_s=12., dt_s=.1)[1]
    explicit = integrate_motion(setup, hydro, excitation=excitation, dt_s=.1, kernel=kernel)
    truncated = integrate_motion(setup, hydro, excitation=excitation, dt_s=.1, kernel=kernel[:101])
    assert_allclose(full.position_m, explicit.position_m, rtol=0., atol=0.)
    assert np.max(abs(full.position_m-truncated.position_m)) > 1e-6
