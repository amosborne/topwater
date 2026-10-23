"""Canonical chart polynomial derivatives and fixed midpoint prolongation."""
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest
from scipy.spatial.transform import Rotation

from topwater import InitialState
from .reference_boat import GRAVITY, WATER_DENSITY, equilibrium
from .experiments.free_surface.boat_action_tangent import object_derivatives, refine_tangent
from .experiments.free_surface.boat_body_action import fluid_gravity
from .experiments.free_surface.boat_chart_tangent import chart_tangent
from .experiments.free_surface.boat_geometry_template import evaluate_template
from .experiments.free_surface.boat_parent_action import parent_mass
from .reference_boat import boundary_snapshot
from .experiments.free_surface.boat_world_field import boundary_data, refine_boundary


def test_actual_chart_load_mass_gravity_all_directions_and_refinement(tmp_path):
    rest = equilibrium()[0]
    pose = InitialState(np.r_[.123, .087, rest.position[2]], rest.orientation)
    _, reference, _ = boundary_snapshot(tmp_path/'tangent.json', pose=pose, export_template=True)
    coarse, tangent = chart_tangent(reference, pose=pose, parent_eta=reference['parent_nodes'][:, 2],
        parent_psi=reference['parent_psi'])
    # Two levels exercise repeated midpoint lineage, including green free
    # faces and inherited trace weights, on the actual canonical boundary.
    current, dx, ds = refine_tangent(coarse, tangent['point_jacobian'], tangent['parent_map_jacobian'], 2)
    derivative = object_derivatives(current, dx, ds)
    twist = np.array([.03, -.02, .01, .01, -.015, .02])
    for index in range(21):
        direction = np.eye(21)[index]
        moved = []
        h = 2e-6
        for offset in (h, -h):
            q = (Rotation.from_rotvec(offset*direction[3:6])*Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
            value = evaluate_template(reference, pose=InitialState(pose.position+offset*direction[:3], q),
                parent_eta=reference['parent_nodes'][:, 2]+offset*direction[6:], parent_psi=reference['parent_psi'])[0]
            for _ in range(2):
                value = refine_boundary(value, region='body')
            moved.append(value)
        assert_array_equal(current['faces'], moved[0]['faces'])
        loads = [boundary_data(value, None, twist, incident_load=False)['body_load'] for value in moved]
        mass = [parent_mass(value) for value in moved]
        gravity = [fluid_gravity(value) for value in moved]
        assert_allclose(derivative['body_load_jacobian'][:, :, index], (loads[0]-loads[1])/(2*h), atol=1e-8, rtol=2e-6)
        assert_allclose(derivative['parent_mass_jacobian'][index], (mass[0]-mass[1])/(2*h), atol=1e-8, rtol=2e-6)
        # U is a large tank first moment; absolute cancellation floor matters.
        assert_allclose(derivative['gravity_jacobian'][index], (gravity[0]-gravity[1])/(2*h), atol=.002, rtol=2e-6)
    mass = parent_mass(current)
    assert_allclose(derivative['gravity_jacobian'][6:], WATER_DENSITY*GRAVITY*(mass@current['parent_nodes'][:, 2]), atol=1e-8, rtol=2e-12)


def test_rejects_missing_or_nonfinite_directions():
    arrays = dict(points=np.zeros((1, 3)), faces=np.zeros((0, 3), dtype=int),
        parent_nodes=np.zeros((1, 3)), free_nodes=np.array([0]))
    with pytest.raises(ValueError, match='all owned directions'):
        object_derivatives(arrays, np.zeros((1, 3, 6)), np.zeros((1, 1, 6)))
