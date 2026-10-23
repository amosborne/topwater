"""Independent open-port endpoint work on the actual canonical boat boundary."""
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from topwater import InitialState
from .reference_boat import WATER_DENSITY
from .experiments.free_surface import boat_parent_action as action
from .experiments.free_surface import boat_world_field as world
from .reference_boat import boundary_snapshot
from .experiments.free_surface.boat_geometry_template import evaluate_template


def test_constant_incident_flow_endpoint_work_and_potential_gauge(tmp_path, monkeypatch):
    _, arrays, _ = boundary_snapshot(tmp_path/'port.json', export_template=True)
    velocity = np.array([.2, -.1, 0.])
    def uniform_flow(points, values, time):
        result = np.zeros((len(points), 9))
        result[:, 6:] = velocity
        return result
    monkeypatch.setattr(world, 'fused_cpu', uniform_flow)
    monkeypatch.setattr(action, 'fused_cpu', uniform_flow)
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    eta = arrays['parent_nodes'][:, 2]
    direction = np.eye(len(eta))[0]
    h = 1e-5
    moved = [evaluate_template(arrays, pose=pose, parent_eta=eta+sign*h*direction,
        parent_psi=arrays['parent_psi'])[0] for sign in (1, -1)]
    for current in moved:
        assert_array_equal(current['faces'], arrays['faces'])
        assert_array_equal(current['free_nodes'], arrays['free_nodes'])
    data = [world.boundary_data(current, None, np.zeros(6)) for current in moved]
    db = (data[0]['right']-data[1]['right'])/(2*h)
    moving = (moved[0]['points']-moved[1]['points'])/(2*h)
    mapping = world.boundary_data(arrays, None, np.zeros(6))['parent_map']
    a = np.array([.3, -.2, .4])
    for gauge in (0., 1.3):
        potential = arrays['points']@a+gauge
        field = dict(potential=potential, parent_map=mapping)
        # Parent0 is the SW corner. On its two perimeter edges, R0 goes1->0.
        # For affine psi, int psi*R0 ds = L*(2*psi_corner+psi_other)/6.
        parent_phi = arrays['parent_nodes']@a+gauge
        expected = 0.
        for other, outward in ((1, np.array([-1., 0., 0.])), (3, np.array([0., -1., 0.]))):
            length = np.linalg.norm(arrays['parent_nodes'][other, :2]-arrays['parent_nodes'][0, :2])
            expected -= WATER_DENSITY*(outward@velocity)*length*(2*parent_phi[0]+parent_phi[other])/6
        seam = action.port_seam_load(arrays, field, None, 0.)
        assert_allclose(seam[0], expected, rtol=1e-13, atol=1e-10)
        assert_allclose(seam[[4, 7, 10]], 0., atol=1e-12)
        covector = action.port_transport_covector(arrays, potential, None, 0.)
        transported = -WATER_DENSITY*(potential@db+np.sum(covector*moving))
        assert_allclose(transported, expected, atol=5e-5, rtol=0.)
