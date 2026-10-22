"""Independent unreduced-scalar checks of the canonical non-BEM gradient."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits

from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from .experiments.free_surface.audit_canoe_energy import gravity_moments
from .experiments.free_surface.canonical_canoe import body_map, configuration, _surface_maps
from .experiments.free_surface.canonical_geometry import gradient, geometric_jacobians
from .experiments.free_surface.canoe import body_load
from .experiments.free_surface.moving_canoe import make_domain, initial_state, field, surface, geometry
from .experiments.free_surface.pocket import stiffness


@pytest.fixture(scope='module')
def canoe():
    domain = make_domain(12, gpu=False)
    state = initial_state(domain, scale=10.)
    size = len(domain.free_nodes)
    with threadpool_limits(1):
        current = field(domain, state[:3], state[3:7], state[13:13+size],
                        state[13+size:], state[7:13])
        height_rate, _, _ = surface(domain, current, state[13+size:], state[7:13])
        q = configuration(domain, state)
        qdot = np.r_[np.linalg.solve(body_map(q[3:6]), state[7:13]), height_rate]
        actual = gradient(domain, q, qdot, current['potential'])
    return domain, q, qdot, current['potential'], actual


def nonbem_scalar(domain, q, qdot, potential, reference):
    """Direct -L using existing weak B, surface maps and sparse FEM stiffness.

    No BEM assembly or field solve occurs here. Unlike the gradient, gravity
    uses volume moments, body flux uses the weak load, and surface work uses
    the corrected mass/transport matrices.
    """
    orientation = (Rotation.from_rotvec(q[3:6])*Rotation.from_quat(
        domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)
    points, height, body = geometry(domain, q[:3], orientation, q[6:])
    velocity = body_map(q[3:6])@qdot[:6]
    rotation = rotation_matrix(orientation)
    inertia = rotation@domain.body.inertia_cg@rotation.T
    rigid = (.5*domain.body.mass*velocity[:3]@velocity[:3]
        +.5*velocity[3:]@inertia@velocity[3:])
    local_points = points[domain.local_mapping]
    mesh = SimpleNamespace(vertices=local_points, faces=domain.contact.faces)
    panels = TrianglePanels(local_points[domain.contact.faces], 1)
    load = body_load(mesh, panels, domain.contact.body_faces, q[:3])
    body_work = 1000*potential[domain.local_mapping]@load@velocity
    eta_map, body_transport = _surface_maps(domain, points, height, body)
    surface_work = 1000*potential[domain.free_nodes]@(
        eta_map@qdot[6:]+body_transport@velocity)
    fem = 500*potential@stiffness(points, domain.exterior_tets)@potential
    moments = gravity_moments(domain, points)
    gravity = (9810*sum(np.sum(value-rest) for value, rest in zip(moments, reference))
        +domain.body.mass*9.81*(q[2]-domain.rest.position[2]))
    return -rigid-body_work-surface_work+fem+gravity


def test_actual_canoe_gradient_matches_independent_scalar_directions(canoe):
    domain, q, qdot, potential, actual = canoe
    directions = list(np.eye(len(q))[:6])
    xy = domain.reference_points[domain.free_nodes, :2]
    for values in (np.cos(.5*xy[:, 0])*np.cos(xy[:, 1]),
                   np.sin(.7*xy[:, 0]+.4*xy[:, 1])):
        values -= values.mean()
        directions.append(np.r_[np.zeros(6), values/np.max(np.abs(values))])
    reference = gravity_moments(domain, domain.reference_points)
    estimates = []
    with threadpool_limits(1):
        for step in (2e-5, 1e-5):
            estimates.append(np.array([(nonbem_scalar(domain, q+step*v, qdot, potential, reference)
                -nonbem_scalar(domain, q-step*v, qdot, potential, reference))/(2*step)
                for v in directions]))
    # Remove the observed O(h^2) truncation, keeping the tolerance on numerical
    # differentiation rather than allowing a physics-sized gradient discrepancy.
    richardson = (4*estimates[1]-estimates[0])/3
    expected = np.array([actual@direction for direction in directions])
    np.testing.assert_allclose(richardson, expected, rtol=3e-8, atol=5e-8)


def test_exposed_geometry_maps_match_full_configuration_direction(canoe):
    domain, q, qdot, _, _ = canoe
    points, maps = geometric_jacobians(domain, q)
    step = 1e-6
    plus = geometric_jacobians(domain, q+step*qdot)[0]
    minus = geometric_jacobians(domain, q-step*qdot)[0]
    expected = np.column_stack([matrix@qdot for matrix in maps])
    np.testing.assert_allclose((plus-minus)/(2*step), expected, atol=3e-10, rtol=3e-7)
    assert points.shape == domain.reference_points.shape
