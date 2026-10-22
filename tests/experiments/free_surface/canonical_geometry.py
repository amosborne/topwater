"""Non-BEM configuration gradient of fixed-field, on-shell reduced Luke -L.

Nodal potential and canonical configuration velocity stay fixed. Surface mesh
velocity depends on the configuration through its geometry Jacobian. Schwarz
symmetry contracts that Jacobian's derivative with qdot using two geometry
evaluations, rather than a full mesh Hessian or one field solve per coordinate.
"""
import numpy as np
from scipy.sparse import csr_matrix, hstack
from scipy.spatial.transform import Rotation

from topwater.frames import rotation_matrix
from .canonical_canoe import body_map
from .moving_canoe import geometry

RHO = 1000.
GRAVITY = 9.81


def _orientation(domain, rotvec):
    return (Rotation.from_rotvec(rotvec)*Rotation.from_quat(
        domain.rest.orientation, scalar_first=True)).as_quat(scalar_first=True)


def geometric_jacobians(domain, q):
    """Return global points and three CSR dX_axis/dq maps; q=(position,r,eta)."""
    orientation = _orientation(domain, q[3:6])
    points, height, body = geometry(domain, q[:3], orientation, q[6:])
    mapping = body_map(q[3:6])
    jacobians = tuple(hstack((csr_matrix(body[:, axis]@mapping), height[axis]),
        format='csr') for axis in range(3))
    return points, jacobians


def _rigid_kinetic(domain, q, qdot):
    velocity = body_map(q[3:6])@qdot[:6]
    rotation = rotation_matrix(_orientation(domain, q[3:6]))
    inertia = rotation@domain.body.inertia_cg@rotation.T
    kinetic = (.5*domain.body.mass*velocity[:3]@velocity[:3]
        +.5*velocity[3:]@inertia@velocity[3:])
    return velocity, kinetic


def _area_gradient(triangles, vector):
    """Vertex gradient of .5*cross(x1-x0,x2-x0).vector at fixed vector."""
    first = .5*np.cross(triangles[:, 2]-triangles[:, 0], vector)
    second = .5*np.cross(vector, triangles[:, 1]-triangles[:, 0])
    return np.stack((-first-second, first, second), axis=1)


def gradient(domain, q, qdot, potential, step=1e-5):
    """Return non-BEM d(-L)/dq, with full nodal potential and qdot fixed.

    Includes rigid kinetic, body-flux and free-surface one-form terms, exterior
    FEM kinetic energy, and total gravity. The separate BEM fixed-trace kinetic
    shape gradient must be added. ``step`` controls only two geometry-Jacobian
    evaluations along qdot and three cheap explicit rotation differences.
    """
    if not np.isfinite(step) or step <= 0:
        raise ValueError('geometry difference step must be finite and positive')
    points, jacobians = geometric_jacobians(domain, q)
    velocity, _ = _rigid_kinetic(domain, q, qdot)
    moving = np.column_stack([matrix@qdot for matrix in jacobians])

    def scatter(faces, values):
        nodal = np.zeros_like(points)
        for corner in range(3):
            np.add.at(nodal, faces[:, corner], values[:, corner])
        return nodal

    def pullback(nodal):
        return sum(matrix.T@nodal[:, axis] for axis, matrix in enumerate(jacobians))

    # Explicit rotation changes both world velocity K(r)*ydot and inertia.
    result = np.zeros(len(q))
    velocity_derivative = np.empty((6, 3))
    for axis in range(3):
        offset = np.zeros(len(q))
        offset[3+axis] = step
        plus = _rigid_kinetic(domain, q+offset, qdot)
        minus = _rigid_kinetic(domain, q-offset, qdot)
        velocity_derivative[:, axis] = (plus[0]-minus[0])/(2*step)
        result[3+axis] = -(plus[1]-minus[1])/(2*step)

    # -rho*phi.T*B*V: differentiate area, weighted position, center and V.
    body_faces = domain.local_mapping[domain.contact.faces[domain.contact.body_faces]]
    triangles = points[body_faces]
    area = .5*np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    trace = potential[body_faces]
    mean = trace.mean(axis=1)
    coefficient = (trace.sum(axis=1)[:, None]+trace)/12
    weighted = np.einsum('fc,fci->fi', coefficient, triangles)-q[:3]*mean[:, None]
    vector = velocity[:3]*mean[:, None]+np.cross(velocity[3:], weighted)
    vertex = (_area_gradient(triangles, vector)
        +coefficient[:, :, None]*np.cross(area, velocity[3:])[:, None])
    result += pullback(scatter(body_faces, -RHO*vertex))
    result[:3] += RHO*np.sum(mean[:, None]*np.cross(area, velocity[3:]), axis=0)
    velocity_gradient = -RHO*np.r_[np.sum(area*mean[:, None], axis=0),
        np.sum(np.cross(weighted, area), axis=0)]
    result[3:6] += velocity_derivative.T@velocity_gradient

    # -rho*integral_FS psi*w.n. dS/dq = J.T*gX + dot(J).T*gW.
    _, plus = geometric_jacobians(domain, q+step*qdot)
    _, minus = geometric_jacobians(domain, q-step*qdot)
    dot_jacobians = tuple((p-m)/(2*step) for p, m in zip(plus, minus))
    faces = domain.free_faces
    triangles = points[faces]
    area = .5*np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    trace = potential[faces]
    coefficient = (trace.sum(axis=1)[:, None]+trace)/12
    vector = np.einsum('fc,fci->fi', coefficient, moving[faces])
    position_gradient = scatter(faces, _area_gradient(triangles, vector))
    moving_gradient = scatter(faces, coefficient[:, :, None]*area[:, None])
    result -= RHO*(pullback(position_gradient)+sum(matrix.T@moving_gradient[:, axis]
        for axis, matrix in enumerate(dot_jacobians)))

    # Fixed nodal potential FEM energy: volume*(.5*u^2*I-u outer u)*grad(N).
    tets = domain.exterior_tets
    coordinates = points[tets]
    edges = coordinates[:, 1:]-coordinates[:, :1]
    shape_gradient = np.empty((len(tets), 4, 3))
    shape_gradient[:, 1:] = np.linalg.inv(edges).transpose(0, 2, 1)
    shape_gradient[:, 0] = -shape_gradient[:, 1:].sum(axis=1)
    volume = np.abs(np.linalg.det(edges))/6
    fluid = np.einsum('tci,tc->ti', shape_gradient, potential[tets])
    stress = (.5*np.sum(fluid*fluid, axis=1)[:, None, None]*np.eye(3)
        -fluid[:, :, None]*fluid[:, None, :])
    element = RHO*volume[:, None, None]*np.einsum('tij,tcj->tci', stress, shape_gradient)
    nodal = np.zeros_like(points)
    for corner in range(4):
        np.add.at(nodal, tets[:, corner], element[:, corner])
    result += pullback(nodal)

    # Reynolds gravity work on hull and all FS faces. Fixed tank walls have
    # zero normal mesh motion; the internal BEM/FEM interface cancels exactly.
    faces = np.vstack((body_faces, domain.free_faces))
    triangles = points[faces]
    area = .5*np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    z = triangles[:, :, 2]
    vertex = RHO*GRAVITY*area[:, None]*(z.sum(axis=1)[:, None, None]+z[:, :, None])/12
    result += pullback(scatter(faces, vertex))
    result[2] += domain.body.mass*GRAVITY
    return result
