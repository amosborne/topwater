"""Independent analytic-field and Duffy-integration checks on real hull data."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial import ConvexHull
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.reference_boat import simulation_inputs
from tests.experiments.free_surface.boat_hos_continuation import tilted_geometry
from tests.experiments.free_surface.body_internal_field import internal_field, single_first_moment


def test_first_moment_and_all_six_loads_against_independent_duffy():
    triangles = np.array([[[0., 0., 0.], [1.2, .1, .2], [.2, .9, -.1]],
                          [[-.2, .3, 1.], [.9, -.1, 1.1], [.5, .8, .8]]])
    panels = TrianglePanels(triangles, -1)
    # Well separated: fixed Gaussian quadrature is appropriate for this control.
    targets = np.array([[.2, -.4, -2.], [-.7, .2, -1.5],
                        [.8, .1, -1.8], [-.3, -.5, -2.2]])
    samples, weights = triangle_samples(triangles, 24)
    kernel = 1/(4*np.pi*np.linalg.norm(samples[None]-targets[:, None, None], axis=-1))
    single = np.einsum('tfq,q,f->tf', kernel, weights, panels.areas)
    first = np.einsum('tfq,fqi,q,f->tfi', kernel, samples, weights, panels.areas)
    actual_s, actual_first = single_first_moment(targets, panels)
    assert_allclose(actual_s, single, atol=2e-15, rtol=3e-13)
    assert_allclose(actual_first, first, atol=2e-15, rtol=3e-13)
    cg = np.array([.3, -.2, .7])
    from tests.experiments.free_surface.body_internal_field import InternalField
    fixture = InternalField(targets, np.array([[1., 0.], [-1., 0.], [0., 1.], [0., -1.]]), targets.mean(axis=0), 1.)
    potential = fixture.field(samples)[0]
    modes = np.concatenate((np.broadcast_to(panels.normals[:, None], samples.shape),
                            np.cross(samples-cg, panels.normals[:, None])), axis=-1)
    expected = -1025*np.einsum('fqc,fqi,q,f->ic', potential, modes, weights, panels.areas)
    assert_allclose(fixture.pressure_load(panels, cg, 1025.), expected, atol=2e-12, rtol=3e-13)


def test_actual_tilted_hull_sources_inside_and_gradient_independent_fd():
    inputs = simulation_inputs()
    triangles, _, _ = tilted_geometry(inputs)
    fixture = internal_field(triangles)
    hull = ConvexHull(np.unique(triangles.reshape(-1, 3), axis=0))
    distances = -(fixture.points@hull.equations[:, :3].T+hull.equations[:, 3])
    assert distances.min() >= .74*fixture.inscribed_radius
    assert_allclose(fixture.strengths.sum(axis=0), 0., atol=0)
    targets = triangles.mean(axis=1)[::9]
    _, gradient = fixture.field(targets)
    for axis in range(3):
        direction = np.eye(3)[axis]
        h = 2e-5
        coarse = (fixture.field(targets+h*direction)[0]-fixture.field(targets-h*direction)[0])/(2*h)
        fine = (fixture.field(targets+h/2*direction)[0]-fixture.field(targets-h/2*direction)[0])/h
        assert np.max(abs(fine-gradient[..., axis])) <= .35*np.max(abs(coarse-gradient[..., axis]))+2e-12
    assert np.isfinite(fixture.pressure_load(TrianglePanels(triangles, -1), inputs.initial_state.position, 1025.)).all()
