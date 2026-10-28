"""Convex standard V-hull and a distributed structure plus two payload masses.

Shape controls come from marine-gyro's offshore-warped-v/2, revision
651d77d169844295b7c9abea1f8475ed4d064019 (src/+hull/standardShape.m).
The convex envelope, rather than interpolated deadrise, defines this surface.
The mass law and effective payload locations are Topwater assumptions, not
manufacturer measurements. All integration uses this same convex geometry.
"""

import warnings

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from ._validation import array
from .workflow_data import Boat


_XI = np.array([0., .20, .40, .60, .75, .85, .95, 1.])
_SHEER = np.array([.91, .99, 1., .92, .70, .47, .18, 0.])
_CHINE = np.array([.83, .86, .88, .80, .61, .40, .142, 0.])
_KEEL = np.array([.025, .010, 0., .008, .040, .120, .358, .580])
_TOP = np.array([.47, .49, .52, .55, .57, .58, .58, .58])
_WARP_DEG = np.array([0., .5, 2., 5., 10., 17., 27., 29.])
_PAYLOAD_POSITIONS = np.array([[.05, 0., .03], [.50, 0., .20]])


def _geometry(beam_over_length, deadrise_deg):
    """Return the closed envelope and its nonphysical upper facets, scaled by L."""
    beam = beam_over_length
    ys, yc = beam*_SHEER/2, beam*_CHINE/2
    zk, zs = beam*_KEEL, beam*_TOP
    beta = deadrise_deg+_WARP_DEG
    zc = zk+yc*np.tan(np.deg2rad(beta))
    if np.any(beta >= 90) or np.any(zc > zs+1e-12*beam):
        raise ValueError('deadrise is incompatible with the standard hull shape')
    ring = np.stack((np.stack((np.zeros(8), zk), axis=1),
                     np.stack((yc, zc), axis=1), np.stack((ys, zs), axis=1),
                     np.stack((-ys, zs), axis=1), np.stack((-yc, zc), axis=1)), axis=1)
    cloud = np.concatenate((np.broadcast_to(_XI[:, None, None], (8, 5, 1)), ring), axis=2)
    # Every profile converges at the bow. Give Qhull one unique bow vertex.
    points = np.vstack((cloud[:-1].reshape(-1, 3), cloud[-1, 0]))
    sheer = np.r_[np.tile([False, False, True, True, False], 7), True]
    try:
        envelope = ConvexHull(points)
    except QhullError as error:
        raise ValueError('specifications do not define a three-dimensional standard hull') from error
    faces = envelope.simplices.copy()
    triangles = points[faces]
    normal = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    flip = np.einsum('ij,ij->i', normal, envelope.equations[:, :3]) < 0
    faces[flip] = faces[flip][:, [0, 2, 1]]
    caps = np.all(sheer[faces], axis=1)
    if not np.any(caps) or np.any(envelope.equations[caps, 2] <= 0):
        raise ValueError('standard hull must have an upward sheer closure')
    return points, faces, envelope.equations[caps], caps


def _section(points, edges, xi):
    """Intersect the convex surface with x=xi, returning its ordered yz polygon."""
    a, b = points[edges[:, 0]], points[edges[:, 1]]
    crossing = (a[:, 0] < xi) & (b[:, 0] > xi) | (b[:, 0] < xi) & (a[:, 0] > xi)
    on = points[points[:, 0] == xi, 1:]
    a, b = a[crossing], b[crossing]
    cuts = a[:, 1:]+((xi-a[:, 0])/(b[:, 0]-a[:, 0]))[:, None]*(b[:, 1:]-a[:, 1:])
    yz = np.unique(np.vstack((on, cuts)), axis=0)
    if len(yz) < 3:
        return yz  # The bow converges to one point.
    return yz[ConvexHull(yz).vertices]


def _structure_moments(points, faces, cap_planes, *, stations=41):
    """Unit-mass first/raw-second moments, using swept section-line quadrature.

    As in marine-gyro, 70% of structural mass follows bottom/sides and 30%
    follows a deck at 60% of local section depth. This is not a laminate model.
    Coordinates and moments are normalized by L and L squared, respectively.
    """
    edges = np.unique(np.sort(np.concatenate((faces[:, [0, 1]], faces[:, [1, 2]],
                                             faces[:, [2, 0]])), axis=1), axis=0)
    regular = np.linspace(0., 1., stations)
    # Retain exact controls; np.unique alone can retain near-identical stations.
    regular = regular[np.min(abs(regular[:, None]-_XI), axis=1) > 64*np.finfo(float).eps]
    xi = np.sort(np.r_[_XI, regular])
    dx = np.diff(xi)
    weights = np.r_[dx[0], dx[:-1]+dx[1:], dx[-1]]/2
    groups = [[], []]
    for x, weight in zip(xi, weights):
        polygon = _section(points, edges, x)
        if len(polygon) < 3:
            continue
        following = np.roll(polygon, -1, axis=0)
        for a, b in zip(polygon, following):
            endpoints = np.column_stack((np.full(2, x), [a, b]))
            # The artificial upper closure carries no structural mass.
            on_cap = abs(endpoints @ cap_planes[:, :3].T+cap_planes[:, 3]) <= 1e-10
            if not np.any(np.all(on_cap, axis=0)):
                groups[0].append((endpoints[0], endpoints[1], weight))
        z = polygon[:, 1].min()+.60*np.ptp(polygon[:, 1])
        crossing = ((polygon[:, 1] <= z) & (following[:, 1] > z)
                    | (following[:, 1] <= z) & (polygon[:, 1] > z))
        a, b = polygon[crossing], following[crossing]
        y = a[:, 0]+(z-a[:, 1])/(b[:, 1]-a[:, 1])*(b[:, 0]-a[:, 0])
        if len(y) != 2:
            raise ValueError('standard hull section must contain one internal deck segment')
        groups[1].append((np.array([x, y.min(), z]), np.array([x, y.max(), z]), weight))
    first, second = np.zeros(3), np.zeros((3, 3))
    for fraction, group in zip((.70, .30), groups):
        a, b, w = zip(*group)
        a, b = np.asarray(a), np.asarray(b)
        measure = np.linalg.norm(b-a, axis=1)*w
        dm = fraction*measure/measure.sum()
        first += np.sum(dm[:, None]*(a+b)/2, axis=0)
        raw = ((np.einsum('ni,nj->nij', a, a)+np.einsum('ni,nj->nij', b, b))/3
               +(np.einsum('ni,nj->nij', a, b)+np.einsum('ni,nj->nij', b, a))/6)
        second += np.einsum('n,nij->ij', dm, raw)
    # The authored model is exactly port/starboard balanced.
    first[1] = 0.
    second[0, 1] = second[1, 0] = second[1, 2] = second[2, 1] = 0.
    return first, second


def generate_boat(length_m, beam_m, deadrise_deg, mass_kg):
    """Build a standard boat from four finite SI specifications (angle in degrees).

    x is forward from transom, y port and z up from the profile's lowest keel
    datum. Structural mass is 0.1*L_ft**3+70 kg. The remaining operating mass is
    split equally between (0.05L,0,0.03L) and (0.50L,0,0.20L), with no intrinsic
    lump inertia. Rejects a nonpositive payload budget or unstable/flooded
    calm-water equilibrium. No MATLAB, catalog, or external checkout is needed.
    """
    length, beam, deadrise, mass = [array(value, (), name).item() for value, name in
        ((length_m, 'length_m'), (beam_m, 'beam_m'),
         (deadrise_deg, 'deadrise_deg'), (mass_kg, 'mass_kg'))]
    if length <= 0 or beam <= 0 or beam >= length:
        raise ValueError('length and beam must be positive, with beam smaller than length')
    if not 0 <= deadrise < 60:
        raise ValueError('deadrise must be in [0, 60) degrees')
    if mass <= 0:
        raise ValueError('mass must be positive')
    length_ft = length/.3048
    with np.errstate(over='raise', invalid='raise'):
        structure_mass = .1*np.float64(length_ft)**3+70.
    payload_mass = mass-structure_mass
    if payload_mass <= 0:
        raise ValueError(f'total mass must exceed hull/deck mass ({structure_mass:g} kg)')
    points, faces, cap_planes, caps = _geometry(beam/length, deadrise)
    first, second = _structure_moments(points, faces, cap_planes)
    positions = _PAYLOAD_POSITIONS
    cg = (structure_mass*first+payload_mass/2*positions.sum(axis=0))/mass
    raw = structure_mass*second+payload_mass/2*(positions.T @ positions)
    central = raw-mass*np.outer(cg, cg)
    inertia = length**2*(np.trace(central)*np.eye(3)-central)
    physical = faces[~caps]
    used, inverse = np.unique(physical, return_inverse=True)
    boat = Boat(points[used]*length, inverse.reshape(-1, 3), mass, cg*length, inertia)
    # Import only the shared, Python-only static check, never a solver backend.
    from ._wecsim_setup import floating_equilibrium
    floating_equilibrium(boat.to_body())
    if not 20 <= length_ft <= 26 or not 18 <= deadrise <= 25:
        warnings.warn('standard boat is outside the nominal 20-26 ft, 18-25 degree family; '
                      'geometry and mass assumptions are extrapolated', UserWarning, stacklevel=2)
    return boat
