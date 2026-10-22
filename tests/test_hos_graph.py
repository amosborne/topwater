"""Manufactured harmonic continuation, independent of HOS recursion identities."""
import numpy as np
from numpy.testing import assert_allclose

from tests.experiments.free_surface.hos_graph import graph_field


LENGTHS = (96., 80.)


def harmonic_graph(n, steepness=.12):
    x, y = np.meshgrid(np.arange(n)*LENGTHS[0]/n, np.arange(n)*LENGTHS[1]/n)
    a, b = 2*np.pi/LENGTHS[0], 2*np.pi/LENGTHS[1]
    eta = steepness/np.hypot(a, b)*(np.cos(a*x+b*y)+.23*np.sin(2*a*x-b*y))
    ex = steepness/np.hypot(a, b)*(-a*np.sin(a*x+b*y)+.46*a*np.cos(2*a*x-b*y))
    ey = steepness/np.hypot(a, b)*(-b*np.sin(a*x+b*y)-.23*b*np.cos(2*a*x-b*y))
    psi = np.zeros_like(eta)
    velocity = np.zeros((*eta.shape, 3))
    for mx, my, amplitude, phase in [(1, 1, 1.7, .2), (2, -1, .31, -.4)]:
        kx, ky = mx*a, my*b
        k = np.hypot(kx, ky)
        angle = kx*x+ky*y+phase
        decay = amplitude*np.exp(k*eta)
        psi += decay*np.cos(angle)
        velocity[..., 0] -= kx*decay*np.sin(angle)
        velocity[..., 1] -= ky*decay*np.sin(angle)
        velocity[..., 2] += k*decay*np.cos(angle)
    flux = velocity[..., 2]-ex*velocity[..., 0]-ey*velocity[..., 1]
    return x, y, eta, psi, velocity, flux


def continued_trace(field, x, y, eta):
    # Direct harmonic Fourier continuation at sampled graph points: no Taylor
    # recursion or cancellation at z=eta is used in this boundary-residual check.
    plane = sum(field.flat_terms)
    ny, nx = plane.shape
    kx, ky = np.meshgrid(2*np.pi*np.fft.fftfreq(nx, d=LENGTHS[0]/nx),
                         2*np.pi*np.fft.fftfreq(ny, d=LENGTHS[1]/ny))
    spectrum = np.fft.fft2(plane)/(nx*ny)
    return np.array([np.sum(spectrum*np.exp(np.hypot(kx, ky)*z+1j*(kx*a+ky*b))).real
                     for a, b, z in zip(x.ravel()[::31], y.ravel()[::31], eta.ravel()[::31])])


def test_curved_graph_order_gain_against_exact_two_oblique_harmonics():
    x, y, eta, psi, exact, flux = harmonic_graph(32)
    errors, residuals = [], []
    for order in range(1, 6):
        field = graph_field(eta, psi, LENGTHS, order=order, padding=2)
        errors.append(np.linalg.norm(field.velocity-exact)/np.linalg.norm(exact))
        residuals.append(np.linalg.norm(continued_trace(field, x, y, eta)-psi.ravel()[::31])
                         /np.linalg.norm(psi.ravel()[::31]))
        assert np.linalg.norm(field.normal_flux-flux)/np.linalg.norm(flux) < 2*errors[-1]
    print("order 1..5 velocity / independently continued boundary errors:", errors, residuals)
    assert np.all(np.diff(errors) < 0)
    assert np.all(np.diff(residuals) < 0)
    assert errors[-1] < errors[0]/100
    assert residuals[-1] < residuals[0]/100


def test_flat_graph_exact_velocity_and_gauge():
    _, _, eta, psi, velocity, flux = harmonic_graph(24, steepness=0)
    field = graph_field(eta, psi+13.2, LENGTHS, order=5)
    assert_allclose(field.velocity, velocity, atol=2e-14)
    assert_allclose(field.normal_flux, flux, atol=2e-14)


def test_spatial_and_padding_refinement_against_same_physical_graph():
    errors = []
    for n in (8, 16, 32):
        _, _, eta, psi, velocity, _ = harmonic_graph(n)
        errors.append([np.linalg.norm(graph_field(eta, psi, LENGTHS, order=5, padding=p).velocity
                                      -velocity)/np.linalg.norm(velocity) for p in (1, 2, 3)])
    print("grids 8/16/32; padding 1/2/3 velocity errors:", errors)
    assert errors[-1][1] < errors[0][1]/20
    assert abs(errors[-1][2]-errors[-1][1]) < 1e-8
    # An underresolved grid can hide error by alias cancellation: n=8 happens
    # to look better without padding. Padding cannot recover missing input
    # frequencies. At n=16 the resolved product error decreases with padding.
    assert errors[1][1] < errors[1][0]
    assert abs(errors[1][2]-errors[1][1]) < 1e-10
