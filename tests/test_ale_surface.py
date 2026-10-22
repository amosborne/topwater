"""Independent uniform-potential-flow checks on the moving canonical canoe."""
import numpy as np
import pytest
from scipy.spatial.transform import Rotation
from scipy.optimize import root

from tests.laplace_reference import triangle_samples
from tests.test_canoe_contact import canoe  # Shared canonical contact fixture.
from tests.experiments.free_surface.ale_surface import rates


def test_flat_constant_potential_is_stationary(canoe):
    _, rest, _, contact = canoe
    free = contact.free_nodes
    geometry = contact.evaluate(rest.position, rest.orientation, np.zeros(len(free)))
    height_rate, potential_rate, velocity = rates(
        geometry.points, contact.faces[contact.free_faces], free,
        np.full(len(free), 2.3), np.zeros(len(free)),
        geometry.height_jacobians, geometry.body_jacobian, np.zeros(6))
    np.testing.assert_allclose(height_rate, 0., atol=1e-14)
    np.testing.assert_allclose(potential_rate, 0., atol=1e-14)
    np.testing.assert_allclose(velocity, 0., atol=1e-14)


@pytest.mark.parametrize('damping',[0.,.4])
def test_uniform_flow_tilted_water_and_all_six_body_velocities(canoe,damping):
    """A plane and constant 3-D fluid velocity have an exact ALE solution.

    The graph has exact Eulerian rate q=W-slope.U. An independent difference of
    actual mesh coordinates supplies its ALE transport, including all nonlocal
    contact and interior-hull motion. No production transport matrix constructs
    the expected kinematic constraint or Bernoulli rate.
    """
    _, rest, _, contact = canoe
    free = contact.free_nodes
    faces = contact.faces[contact.free_faces]
    position = rest.position+np.array([.003, -.002, .0002])
    orientation = Rotation.from_rotvec([.001, -.0007, .0013]).as_quat(scalar_first=True)
    slope = np.array([.0004, -.0003])
    offset = .0002
    eta = contact.reference_points[free, :2]@slope+offset
    controls=np.searchsorted(free,contact.contact_nodes)
    def plane_error(values):
        heights=eta.copy();heights[controls]=values
        points=contact.evaluate(position,orientation,heights).points[contact.contact_nodes]
        return points[:,2]-points[:,:2]@slope-offset
    solved=root(plane_error,eta[controls])
    assert solved.success
    eta[controls]=solved.x
    geometry = contact.evaluate(position, orientation, eta)
    np.testing.assert_allclose(geometry.points[free,2],
        geometry.points[free,:2]@slope+offset,atol=2e-13,rtol=0)
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
        geometry.height_jacobians, geometry.body_jacobian, twist,
        damping=np.full(len(free),damping))
    exact_velocity=_geometry_difference(contact,position,orientation,eta,height_rate,twist)
    np.testing.assert_allclose(height_rate-exact_velocity[free,:2]@slope,
        q-damping*eta,atol=3e-9,rtol=0)
    fluid_speed_squared = horizontal@horizontal+vertical**2
    expected_potential = (-9.81*eta-.5*fluid_speed_squared
        +exact_velocity[free, :2]@horizontal+vertical*exact_velocity[free, 2])
    # Damping both canonical variables introduces artificial pressure in the
    # absorption zone. Its Eulerian correction is sigma*(eta*u_z-psi); the
    # physical uniform fluid velocity itself is unchanged at this instant.
    expected_potential+=damping*(eta*vertical-psi)
    np.testing.assert_allclose(height_rate, exact_velocity[free, 2], atol=3e-9, rtol=0)
    np.testing.assert_allclose(mesh_velocity, exact_velocity, atol=3e-8, rtol=2e-7)
    np.testing.assert_allclose(potential_rate, expected_potential, atol=3e-9, rtol=0)


def _geometry_difference(contact,position,orientation,eta,height_rate,twist):
    step=1e-6
    rotation=Rotation.from_quat(orientation,scalar_first=True)
    moved=[]
    for dt in (step,-step):
        q=(Rotation.from_rotvec(dt*twist[3:])*rotation).as_quat(scalar_first=True)
        moved.append(contact.evaluate(position+dt*twist[:3],q,eta+dt*height_rate).points)
    return (moved[0]-moved[1])/(2*step)


def test_coupled_height_transport_recovers_independent_geometric_flux(canoe):
    """Nonflat graph, independent height rates and six-axis rigid motion."""
    _,rest,_,contact=canoe
    free=contact.free_nodes;faces=contact.faces[contact.free_faces]
    position=rest.position+np.array([.002,-.001,.001])
    orientation=Rotation.from_rotvec([.003,-.002,.001]).as_quat(scalar_first=True)
    x,y=contact.reference_points[free,:2].T
    eta=.001*np.cos(.6*x+.4*y)
    prescribed=.002*np.cos(.5*x)+.001*np.sin(.3*y)
    twist=np.array([.01,-.015,.02,.03,-.02,.015])
    geometry=contact.evaluate(position,orientation,eta)
    moving=_geometry_difference(contact,position,orientation,eta,prescribed,twist)
    points=geometry.points[faces]
    normal=np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0])
    slope=-normal[:,:2]/normal[:,2,None]
    bary,weights=triangle_samples(np.eye(3)[None],5);bary=bary[0]
    indices=np.searchsorted(free,faces)
    graph_rate=np.einsum('qc,fc->fq',bary,prescribed[indices])
    mesh_motion=np.einsum('qc,fci->fqi',bary,moving[faces,:2])
    eulerian_rate=graph_rate-np.einsum('fi,fqi->fq',slope,mesh_motion)
    flux=np.zeros(len(free))
    for corner in range(3):
        np.add.at(flux,indices[:,corner],np.einsum('f,q,q,fq->f',
            normal[:,2]/2,weights,bary[:,corner],eulerian_rate))
    height_rate,_,velocity=rates(geometry.points,faces,free,np.zeros(len(free)),flux,
        geometry.height_jacobians,geometry.body_jacobian,twist)
    np.testing.assert_allclose(height_rate,prescribed,atol=3e-9,rtol=0)
    np.testing.assert_allclose(velocity,moving,atol=3e-8,rtol=2e-7)
