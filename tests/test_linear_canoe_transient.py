"""Independent startup solution and actual-tank frequency checks of the control."""
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .experiments.free_surface.linear_canoe_transient import LinearTank, prepare, integrate
from .experiments.free_surface.linear_canoe_radiation import solve
from .experiments.free_surface.moving_canoe import make_domain


def test_starting_fluid_has_the_exact_forced_undamped_modal_transient():
    # An analytically forced water mode, including its persistent startup mode.
    # A steady harmonic predictor misses the cos(Omega*t) term entirely.
    mass, lam, omega = 2.3, .8, 1.7
    coupling = np.array([[.3, -.5, .8, -.2, .1, .4]])
    added = np.diag([2., 3., 4., 1., 2., 3.])
    restoring = np.diag([.2, .4, 8., 2., 3., .1])
    buoyancy = np.array([0., 0., 100., 0., 0., 0.])
    model = LinearTank(csc_matrix([[mass]]), csc_matrix((1, 1)),
        np.array([[mass*lam]]), coupling, added, restoring, buoyancy)
    direction = np.array([.2, -.1, .4, .05, -.03, .07])

    def motion(time):
        return np.outer([np.cos(omega*time), -omega*np.sin(omega*time),
                         -omega**2*np.cos(omega*time)], direction)

    times = np.linspace(0., 8., 133)
    actual = integrate(model, times, motion, rtol=1e-11, atol=1e-13)
    natural = np.sqrt(9.81*lam)
    coefficient = float((coupling@direction)[0])*omega**2/mass/(natural**2-omega**2)
    eta = coefficient*(np.cos(omega*times)-np.cos(natural*times))
    eta_rate = coefficient*(-omega*np.sin(omega*times)+natural*np.sin(natural*times))
    motion_values = np.array([motion(t) for t in times])
    psi = (eta_rate+motion_values[:, 1]@coupling[0]/mass)/lam
    expected_load = (buoyancy-motion_values[:, 0]@restoring.T-motion_values[:, 2]@added.T
                     +9810.*eta[:, None]*coupling)
    np.testing.assert_allclose(actual['eta'][:, 0], eta, atol=2e-12, rtol=2e-10)
    np.testing.assert_allclose(actual['psi'][:, 0], psi, atol=2e-12, rtol=2e-10)
    np.testing.assert_allclose(actual['hydrodynamic_load'], expected_load, atol=2e-8, rtol=2e-10)


def test_actual_tank_matches_existing_full_frequency_system_with_weak_sponge():
    omega = 4.
    xi = np.array([.4+.2j, -.3+.5j, .7-.4j, .1+.3j, -.2-.1j, .3+.2j])
    velocity, acceleration = -1j*omega*xi, -omega**2*xi
    with threadpool_limits(1):
        domain = make_domain(12, False)
        for strength in (0., 4.):
            model = prepare(domain, sponge_strength=strength)
            report, saved = solve(domain, [omega], sponge_strength=strength)
            psi = saved['potential_per_velocity'][0][domain.free_nodes]@velocity
            eta = saved['height_per_velocity'][0]@velocity
            mass = splu(model.mass.astype(complex).tocsc())
            eta_rate = mass.solve(model.schur@psi-model.coupling@velocity-model.damping@eta)
            psi_rate = -9.81*eta-mass.solve(model.damping@psi)
            np.testing.assert_allclose(eta_rate, -1j*omega*eta, atol=2e-11, rtol=2e-10)
            np.testing.assert_allclose(psi_rate, -1j*omega*psi, atol=2e-11, rtol=2e-10)
            force = -model.restoring@xi-model.added_mass@acceleration-1000.*model.coupling.T@psi_rate
            coefficients = report['cases'][0]
            expected = ((omega**2*np.array(coefficients['added_mass'])
                        +1j*omega*np.array(coefficients['radiation_damping'])-model.restoring)@xi)
            np.testing.assert_allclose(force, expected, atol=3e-8, rtol=3e-10)
            if strength:
                # Consistent projection is not nodal -sigma*psi.
                assert model.damping.nnz > model.mass.shape[0]
