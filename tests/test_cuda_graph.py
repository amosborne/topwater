"""Independent SciPy HOS reference and periodic mass/work checks."""
import numpy as np
from numpy.testing import assert_allclose
import pytest

pytest.importorskip('cupy')
from tests.experiments.free_surface.cuda_graph import (
    GraphDirichlet, ExteriorTrace, body_maps, harmonic_velocity, cp, _resize,
)
from tests.experiments.free_surface.hos_graph import graph_field


def test_exterior_state_is_preserved_and_hidden_trace_is_discrete_harmonic():
    size = 17
    x, y = np.meshgrid(np.arange(size)-8., np.arange(size)-8.)
    inside = (abs(x) < 3) & (abs(y) < 2)
    fill = ExteriorTrace(inside)
    # Both are exactly discrete-harmonic; the hole is away from periodic seams.
    fields = np.column_stack(((.3+x-.2*y).ravel(), (x*x-y*y+.1*x*y).ravel()))
    exterior = fields[fill.exterior]
    actual = cp.asnumpy(fill.extend(exterior))
    assert_allclose(actual[fill.exterior], exterior, atol=0, rtol=0)
    assert_allclose(actual[fill.hidden], fields[fill.hidden], atol=3e-14, rtol=2e-14)
    assert_allclose(cp.asnumpy(fill.extend(exterior[:, 0])), actual[:, 0], atol=2e-15)
    with pytest.raises(ValueError, match='finite float64'):
        fill.extend(np.zeros(len(fill.exterior), dtype=np.int64))


def test_gpu_body_maps_match_independent_complex_fourier_formula():
    from tests.experiments.free_surface.boat_radiation_time import harmonic_maps
    points = np.array([[-.2, .7, -.03], [1.3, -.5, -.8], [.4, .8, 0.]])
    normals = np.array([[.6, 0., .8], [0., -.8, .6], [1., 0., 0.]])
    for size in [17, 35]:
        actual = body_maps(points, normals, size, .3, [-2., -1.])
        expected = harmonic_maps(points, normals, size, .3, [-2., -1.])
        for gpu, cpu in zip(actual, expected[:2], strict=True):
            assert_allclose(cp.asnumpy(gpu), cpu, atol=3e-16, rtol=2e-13)
        length = size*.3
        x, y = cp.meshgrid(cp.arange(size), cp.arange(size))
        spectrum = cp.fft.fft2(cp.cos(2*np.pi*(2*x+y)/size)).ravel()
        wavevector = 2*np.pi*np.array([2., 1.])/length
        wave = np.linalg.norm(wavevector)
        phase = (points[:, :2]-[-2., -1.])@wavevector
        analytic = np.exp(wave*points[:, 2, None])*np.column_stack((
            -wavevector[0]*np.sin(phase), -wavevector[1]*np.sin(phase), wave*np.cos(phase)))
        assert_allclose(cp.asnumpy(harmonic_velocity(actual[0], spectrum, size, length)),
                        analytic, atol=2e-14, rtol=3e-13)


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


def test_actual_hull_source_traces_against_cpu_reference():
    from topwater.rankine import TrianglePanels
    from tests.reference_boat import simulation_inputs
    from tests.experiments.free_surface.boat_hos_continuation import tilted_geometry
    from tests.experiments.free_surface.boat_image_plane import elevated_trace
    from tests.experiments.free_surface.boat_graph_images import image_height

    inputs = simulation_inputs()
    triangles, _, _ = tilted_geometry(inputs)
    panels = TrianglePanels(triangles, -1)
    n, length = 33, 16.
    origin = inputs.initial_state.position[:2]-length/2
    x, y = np.meshgrid(origin[0]+np.arange(n)*length/n, origin[1]+np.arange(n)*length/n)
    level = triangles[..., 2].max()+.08
    eta = level+.02*np.cos(2*np.pi*x/length)+.02*np.sin(2*np.pi*y/length)
    height, _ = image_height(triangles, eta, [length, length], origin)
    trace, _ = elevated_trace(np.column_stack((x.ravel(), y.ravel(), eta.ravel())), panels, height+.15)
    density = np.column_stack((panels.normals[:, 2],
        np.cross(panels.centroids-inputs.initial_state.position, panels.normals)[:, 0]))
    values = trace@density
    graph = GraphDirichlet(eta-level, (length, length), order=4)
    for col in range(2):
        psi = values[:, col].reshape(n, n)
        spectra = graph.coefficients(psi)
        velocity, flux = graph.surface(spectra)
        reference = graph_field(eta-level, psi, (length, length), order=4, padding=2+.25/n)
        assert_allclose(cp.asnumpy(velocity), reference.velocity, atol=3e-12, rtol=3e-11)
        assert_allclose(cp.asnumpy(flux), reference.normal_flux, atol=3e-12, rtol=3e-11)
        for spectrum, term in zip(spectra, reference.flat_terms, strict=True):
            assert_allclose(cp.asnumpy(cp.fft.ifft2(spectrum).real), term, atol=3e-13, rtol=3e-11)
