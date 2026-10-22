"""Independent uniform-potential-flow checks on the moving canonical canoe."""
import numpy as np
from scipy.spatial.transform import Rotation

from topwater.frames import rotation_matrix
from tests.test_canoe_contact import canoe  # Shared canonical contact fixture.
from tests.experiments.free_surface.ale_surface import rates


def test_flat_constant_potential_is_stationary(canoe):
    _, rest, _, contact = canoe
    free = contact.free_nodes
    geometry = contact.evaluate(rest.position, rest.orientation, np.zeros(len(free)))
    height_rate, potential_rate, velocity = rates(
        geometry.points, contact.faces[contact.free_faces], free,
        np.full(len(free), 2.3), np.zeros(len(free)),
        geometry.eta_jacobian, geometry.body_jacobian, np.zeros(6))
    np.testing.assert_allclose(height_rate, 0., atol=1e-14)
    np.testing.assert_allclose(potential_rate, 0., atol=1e-14)
    np.testing.assert_allclose(velocity, 0., atol=1e-14)


def test_uniform_flow_tilted_water_and_all_six_body_velocities(canoe):
    """A plane and constant 3-D fluid velocity have an exact ALE solution.

    The graph moves vertically at q=W-slope.U. Contact velocity comes directly
    from intersecting moving material edges with that plane. The potential rate
    follows Eulerian Bernoulli plus the full three-dimensional mesh transport;
    neither expected value uses the implementation's corrected mass matrix.
    """
    _, rest, _, contact = canoe
    free = contact.free_nodes
    faces = contact.faces[contact.free_faces]
    position = rest.position+np.array([.003, -.002, .0002])
    orientation = Rotation.from_rotvec([.001, -.0007, .0013]).as_quat(scalar_first=True)
    slope = np.array([.0004, -.0003])
    offset = .0002
    eta = contact.reference_points[free, :2]@slope+offset
    world = (contact.material_vertices-contact.cg)@rotation_matrix(orientation).T+position
    a, b = world[contact.contact_edges[:, 0]], world[contact.contact_edges[:, 1]]
    edge = b-a
    # Exact intersection with z=slope.xy+offset; no iterative geometry solve.
    fraction = (a[:, :2]@slope+offset-a[:, 2])/(edge[:, 2]-edge[:, :2]@slope)
    eta[np.searchsorted(free, contact.contact_nodes)] = (a+fraction[:, None]*edge)[:, 2]
    geometry = contact.evaluate(position, orientation, eta)
    horizontal = np.array([.3, -.2])
    vertical = .07
    q = vertical-slope@horizontal
    psi = geometry.points[free, :2]@horizontal+vertical*eta
    # Constant projected normal flux: integral N_i q dxdy = area*q/3.
    flux = np.zeros(len(free))
    triangles = geometry.points[faces]
    areas = np.cross(triangles[:, 1]-triangles[:, 0],
                     triangles[:, 2]-triangles[:, 0])[:, 2]/2
    for corner in range(3):
        np.add.at(flux, np.searchsorted(free, faces[:, corner]), areas*q/3)
    twist = np.array([.02, -.015, .01, .03, -.02, .015])
    height_rate, potential_rate, mesh_velocity = rates(
        geometry.points, faces, free, psi, flux,
        geometry.eta_jacobian, geometry.body_jacobian, twist)

    exact_velocity = np.zeros((len(free), 3))
    exact_velocity[:, 2] = q
    p = geometry.points[contact.contact_nodes]
    rigid_velocity = twist[:3]+np.cross(twist[3:], p-position)
    slide_rate = (q-rigid_velocity[:, 2]+rigid_velocity[:, :2]@slope)/(
        edge[:, 2]-edge[:, :2]@slope)
    exact_velocity[np.searchsorted(free, contact.contact_nodes)] = (
        rigid_velocity+slide_rate[:, None]*edge)
    fluid_speed_squared = horizontal@horizontal+vertical**2
    expected_potential = (-9.81*eta-.5*fluid_speed_squared
        +exact_velocity[:, :2]@horizontal+vertical*exact_velocity[:, 2])
    np.testing.assert_allclose(height_rate, exact_velocity[:, 2], atol=3e-13, rtol=3e-12)
    np.testing.assert_allclose(mesh_velocity[free], exact_velocity, atol=3e-13, rtol=3e-12)
    np.testing.assert_allclose(potential_rate, expected_potential, atol=3e-13, rtol=3e-12)
