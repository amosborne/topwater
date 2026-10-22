"""Independent similarity identities on the actual canonical hull, GPU optional."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip('cupy')
from .experiments.free_surface import cuda_operator as gpu
from .experiments.free_surface.cuda_directional import energy_direction_action
from .experiments.free_surface.boat_domain_inertia import mirrored_body
from .reference_boat import equilibrium
from topwater.rankine import TrianglePanels


def test_canonical_hull_similarity_derivatives_and_uniform_rule_ownership(monkeypatch):
    try:
        if gpu.cp.cuda.runtime.getDeviceCount() == 0:
            pytest.skip('CUDA device unavailable')
    except gpu.cp.cuda.runtime.CUDARuntimeError:
        pytest.skip('CUDA runtime unavailable')
    mesh, _ = mirrored_body(equilibrium()[1].wetted_triangles)
    points, faces = mesh.vertices, mesh.faces
    panels = TrianglePanels(points[faces], 1)
    phi = points@np.array([.3, -.2, .5])
    assembly = gpu.prepare_energy(faces, panels)
    assert assembly.uniform_quadrature
    # For fixed nodal phi, rigid motion leaves E unchanged. Uniform dilation
    # gives E(s*X)=s*E(X), since gradients scale 1/s and volume scales s^3.
    for label, direction, expected in (
            ('translation', np.broadcast_to([.2, -.3, .4], points.shape), np.zeros(len(points))),
            ('rotation', np.cross([.1, -.2, .3], points), np.zeros(len(points))),
            ('dilation', points, assembly.energy@phi)):
        actual, _ = energy_direction_action(points, faces, direction, phi, assembly=assembly)
        assert_allclose(actual, expected, rtol=2e-10, atol=1e-11, err_msg=label)
    with pytest.raises(ValueError, match='identical current'):
        energy_direction_action(points+1e-6, faces, points, phi, assembly=assembly)
    # An alternate integration kernel cannot borrow the uniform-rule derivative,
    # even if this particular wrapper happens to delegate to the same kernel.
    original = gpu._pairs
    with monkeypatch.context() as context:
        context.setattr(gpu, '_pairs', lambda *args: original(*args))
        alternate = gpu.prepare_energy(faces, panels)
    assert not alternate.uniform_quadrature
    with pytest.raises(ValueError, match='uniform quadrature'):
        energy_direction_action(points, faces, points, phi, assembly=alternate)
    with pytest.raises(ValueError):
        alternate.auxiliary(points, faces, phi)
