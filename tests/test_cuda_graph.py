"""Independent SciPy HOS reference and periodic mass/work checks."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')
from tests.experiments.free_surface.cuda_graph import GraphDirichlet, cp, _resize
from tests.experiments.free_surface.hos_graph import graph_field


@pytest.mark.parametrize('order', [2, 4])
def test_resolved_oblique_fields_against_independent_cpu(order):
    n = 33
    x, y = np.meshgrid(2*np.pi*np.arange(n)/n, 2*np.pi*np.arange(n)/n)
    eta = .08*np.cos(x+2*y)+.025*np.sin(2*x-y)
    psi = .13*np.sin(2*x+y)+.04*np.cos(x-2*y)
    graph = GraphDirichlet(eta, (2*np.pi, 2*np.pi), order=order)
    spectra = graph.coefficients(psi)
    velocity, flux = graph.surface(spectra)
    # SciPy ceil(padding*n) must use the SAME odd padded shape, not even 2*n.
    reference = graph_field(eta, psi, (2*np.pi, 2*np.pi), order=order, padding=(2*n+1)/n)
    assert_allclose(cp.asnumpy(velocity), reference.velocity, atol=3e-12, rtol=3e-11)
    assert_allclose(cp.asnumpy(flux), reference.normal_flux, atol=3e-12, rtol=3e-11)
    for spectrum, term in zip(spectra, reference.flat_terms, strict=True):
        assert_allclose(cp.asnumpy(cp.fft.ifft2(spectrum).real), term, atol=3e-13, rtol=3e-11)
    assert abs(float(cp.mean(flux))) < 3e-13
    other = .07*np.cos(x+y)+.03*np.sin(3*x-y)
    other_flux = graph.surface(graph.coefficients(other))[1]
    assert abs(float(cp.mean(cp.asarray(other)*flux-cp.asarray(psi)*other_flux))) < 3e-13


def test_owned_graph_and_odd_fourier_amplitudes():
    n = 17
    x, y = np.meshgrid(2*np.pi*np.arange(n)/n, 2*np.pi*np.arange(n)/n)
    mode = 1+.3*np.cos(3*x+2*y)
    assert_allclose(cp.asnumpy(_resize(_resize(cp.asarray(mode), (35, 35)), (n, n))), mode, atol=2e-15)
    eta = cp.asarray(.04*np.cos(x+y))
    graph = GraphDirichlet(eta, (2*np.pi, 2*np.pi), order=4)
    before = graph.surface(graph.coefficients(mode))[1]
    eta[:] = 10.
    assert_allclose(cp.asnumpy(graph.surface(graph.coefficients(mode))[1]), cp.asnumpy(before), atol=0, rtol=0)
    with pytest.raises(ValueError, match='odd'):
        GraphDirichlet(np.zeros((16, 17)), (1., 1.))
