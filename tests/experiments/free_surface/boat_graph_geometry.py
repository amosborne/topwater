"""Shared moving-plane geometry, rigid kinematics and Fourier dry certificates."""
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
from topwater import InitialState
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics, FloodingError
from topwater.rankine import TrianglePanels
from .boat_spatial_work import subdivide
from .graph_geometry_reference import fourier_points, tensor_fourier


def constant_twist_pose(initial, time):
    """Zero world linear/angular acceleration for the coupled known RHS."""
    initial_rotation = Rotation.from_quat(initial.orientation, scalar_first=True)
    omega = initial_rotation.apply(initial.angular_velocity)
    orientation = Rotation.from_rotvec(omega*time)*initial_rotation
    return InitialState(initial.position+initial.linear_velocity*time,
        orientation.as_quat(scalar_first=True), initial.linear_velocity,
        orientation.inv().apply(omega))


def wet_snapshot(body, pose, plane_origin, frame, grid_xy, length, refinement):
    orientation = Rotation.from_matrix(frame)*Rotation.from_quat(pose.orientation, scalar_first=True)
    state = InitialState(frame@(pose.position-plane_origin), orientation.as_quat(scalar_first=True))
    wet = hydrostatics(body, state, density=WATER_DENSITY, gravity=GRAVITY)
    triangles = wet.wetted_triangles@frame+plane_origin
    for _ in range(refinement):
        triangles = subdivide(triangles)
    footprint = ConvexHull((wet.waterline_segments@frame+plane_origin)[..., :2].reshape(-1, 2))
    distance = grid_xy@footprint.equations[:, :2].T+footprint.equations[:, 2]
    inside = np.all(distance <= 64*np.finfo(float).eps*length, axis=1)
    relative = (body.mesh.vertices-body.cg)@rotation_matrix(pose.orientation).T
    signs = np.sign(((relative+pose.position-plane_origin)@frame.T)[:, 2])
    if np.any(signs == 0):
        raise ValueError('stencil needs persistent clipping recipes, no on-plane input vertex')
    return TrianglePanels(triangles, -1), inside, wet.source_faces, signs


def fourier_dry_certificate(body, world, eta, length, origin):
    """Sufficient Lipschitz dry-rim/cap bounds, with explicit roundoff allowance."""
    size = eta.shape[0]
    k = 2*np.pi*np.fft.fftfreq(size, d=length/size)
    spectrum = np.fft.fft2(eta)/eta.size
    derivative = np.array([np.sum(abs(spectrum*k[None, :])),
                           np.sum(abs(spectrum*k[:, None]))])
    allowance = 64*np.finfo(float).eps*eta.size*max(1., float(np.sum(abs(spectrum))), float(abs(world).max()))
    rim = body.mesh.opening_rim
    lower, sampled = [], []
    for i, j in zip(rim, np.roll(rim, -1)):
        a, b = world[i], world[j]
        points = a+((np.arange(64)+.5)/64)[:, None]*(b-a)
        height = fourier_points(eta, points[:, :2], length, origin)[0]
        gaps = points[:, 2]-height
        bound = .5/64*(abs(b[2]-a[2])+derivative@abs(b[:2]-a[:2]))+allowance
        sampled.append(float(gaps.min()))
        lower.append(float(gaps.min()-bound))
    if sampled and min(sampled) <= 0:
        raise FloodingError('actual Fourier graph reaches a sampled physical rim point')
    if lower and min(lower) <= 0:
        raise ValueError('cannot certify dry rim against Fourier graph; not proven flooding')
    cap_bound = None
    if len(body.mesh.artificial_faces):
        cap = world[np.unique(body.mesh.faces[body.mesh.artificial_faces])]
        x = np.linspace(cap[:, 0].min(), cap[:, 0].max(), 129)
        y = np.linspace(cap[:, 1].min(), cap[:, 1].max(), 65)
        height = tensor_fourier(eta, x, y, length, origin)[0]
        error = .5*derivative@np.array([np.diff(x).max(), np.diff(y).max()])+allowance
        cap_bound = float(cap[:, 2].min()-height.max()-error)
        if cap_bound <= 0:
            raise ValueError('cannot certify artificial closure dry against Fourier graph; not proven flooding')
    return dict(rim_sample_minimum_m=min(sampled) if sampled else None,
        rim_clearance_lower_bound_m=min(lower) if lower else None,
        artificial_closure_clearance_lower_bound_m=cap_bound,
        derivative_absolute_sum_bounds=derivative.tolist(), rounding_allowance_m=allowance,
        scope='Sufficient analytic Lipschitz bounds with engineering floating allowance; not interval arithmetic')


