"""Independent geometry and manufactured-field derivative controls."""
import numpy as np
from numpy.testing import assert_allclose
from scipy.spatial import ConvexHull

from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.reference_boat import simulation_inputs
from tests.experiments.free_surface.boat_hos_continuation import tilted_geometry, image_field


def test_tilted_cut_remains_on_supplied_hull_below_the_true_graph():
    inputs = simulation_inputs()
    triangles, level, slope = tilted_geometry(inputs)
    points = triangles.reshape(-1, 3)
    height = level+(points[:, :2]-inputs.initial_state.position[:2])@slope
    assert np.max(points[:, 2]-height) < 2e-12
    assert np.any(np.abs(points[:, 2]-height) < 2e-12)
    world = ((inputs.body.mesh.vertices-inputs.body.cg)
             @rotation_matrix(inputs.initial_state.orientation).T+inputs.initial_state.position)
    planes = ConvexHull(world).equations
    assert np.max(points@planes[:, :3].T+planes[:, 3]) < 2e-12
    panels = TrianglePanels(triangles, -1)
    inward = panels.centroids+1e-6*panels.normals
    outward = panels.centroids-1e-6*panels.normals
    assert np.max(inward@planes[:, :3].T+planes[:, 3]) < 2e-12
    assert np.all(np.max(outward@planes[:, :3].T+planes[:, 3], axis=1) > 0)


def test_manufactured_image_field_full_gradient():
    triangles, _, _ = tilted_geometry(simulation_inputs())
    panels = TrianglePanels(triangles, -1)
    points = panels.centroids[:3]
    density = panels.normals[:, [0, 2]]
    height = triangles[..., 2].max()+.2
    _, gradient = image_field(points, panels, density, [height, height+.8])
    for axis in range(3):
        delta = np.eye(3)[axis]*1e-5
        difference = (image_field(points+delta, panels, density, [height, height+.8])[0]
                      -image_field(points-delta, panels, density, [height, height+.8])[0])/2e-5
        assert_allclose(gradient[..., axis], difference, rtol=2e-8, atol=2e-11)
