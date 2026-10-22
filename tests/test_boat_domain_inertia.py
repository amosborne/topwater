"""Independent physical moments for the actual boat's mirrored halfspace."""
import numpy as np
from numpy.testing import assert_allclose

from topwater.rankine import TrianglePanels
from .reference_boat import equilibrium
from .reference_hull import surface_integrals
from .experiments.free_surface.boat_domain_inertia import mirrored_body
from .experiments.free_surface.canoe import body_load


def test_mirrored_boat_preserves_lower_geometry_and_correct_half_volume_work():
    pose, hydro = equilibrium()
    mesh, odd = mirrored_body(hydro.wetted_triangles)
    assert_allclose(mesh.vertices[mesh.faces[:len(hydro.wetted_triangles)]],
                    hydro.wetted_triangles, rtol=0, atol=0)
    volume, center = surface_integrals(mesh.vertices, mesh.faces)
    assert_allclose(volume, 2*hydro.displaced_volume, rtol=2e-14)
    assert_allclose(center[2], 0., atol=2e-16)
    lower = np.flatnonzero(mesh.vertices[:, 2] < 0)
    assert_allclose(odd@mesh.vertices[lower, 2], mesh.vertices[:, 2], rtol=0, atol=0)
    assert np.all(odd[mesh.vertices[:, 2] == 0] == 0)
    load = body_load(mesh, TrianglePanels.from_mesh(mesh, normal_sign=-1),
                     np.arange(len(hydro.wetted_triangles)), pose.position)
    # Divergence theorem on the wetted hull closed by z=0: that cap contributes
    # no work to phi=z. Fluid-outward normals point INTO the body, hence -V.
    offset = hydro.center_of_buoyancy-pose.position
    expected = hydro.displaced_volume*np.array([0., 0., -1., -offset[1], offset[0], 0.])
    assert_allclose(load.T@mesh.vertices[:, 2], expected, atol=3e-15, rtol=2e-14)
