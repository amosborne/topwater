"""Contained finite-depth P4 strip measurement of the existing surface sponge.

No canoe, nonlinear model, backend API or outgoing-boundary qualification.
Exp(-i*omega*t) and exp(i*ky*y) leave the FULL elliptic x-z problem. Both
surface equations remain in the augmented system. Modal ports retain the
propagating and a reported/converged number of evanescent vertical modes.
Port integration accuracy and truncation are separate controls.
"""
import argparse
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .p4_tetra import basis, indices, key

GRAVITY = 9.81
DEPTH = 4.


@dataclass(frozen=True)
class Strip:
    points: np.ndarray
    triangles: np.ndarray
    elements: np.ndarray
    top: np.ndarray
    left: np.ndarray
    right: np.ndarray
    xlevels: np.ndarray
    zlevels: np.ndarray
    omega: float
    angle: float
    k: float
    ky: float
    kappa: float


def _gauss(order):
    x, w = leggauss(order)
    return (x+1)/2, w/2


def _line(order):
    t, w = _gauss(order)
    value, _ = basis(np.column_stack((1-t, t)))
    return t, w, value


def make_strip(omega=2., angle=0., width=1.5943951023931955, *,
               spacing=np.pi/6, layers=4, quiet_wavelengths=1.):
    """Fixed affine P4 triangles; theta is measured from the x-normal port."""
    if (not all(np.isfinite(v) for v in (omega, angle, width, spacing, quiet_wavelengths))
            or omega <= 0 or width <= 0 or spacing <= 0 or quiet_wavelengths <= 0
            or not 0 <= angle <= np.pi/4 or layers not in (4, 8)):
        raise ValueError('positive strip parameters, angle0..pi/4 and4/8 graded layers required')
    k = brentq(lambda v: GRAVITY*v*np.tanh(DEPTH*v)-omega**2, 1e-9,
               max(10., 2*omega**2/GRAVITY))
    ky, kappa = k*np.sin(angle), k*np.cos(angle)
    quiet = quiet_wavelengths*2*np.pi/kappa
    x = np.r_[np.linspace(-quiet, 0., int(np.ceil(quiet/spacing))+1),
              np.linspace(0., width, int(np.ceil(width/spacing))+1)[1:]]
    z = -DEPTH*(1-np.linspace(0, 1, layers+1))**2
    points = np.array([(a, b) for a in x for b in z])
    ids = np.arange(len(points)).reshape(len(x), len(z))
    triangles = []
    for i in range(len(x)-1):
        for j in range(layers):
            a, b, c, d = ids[i, j], ids[i+1, j], ids[i+1, j+1], ids[i, j+1]
            triangles.extend(((a, b, c), (a, c, d)))
    triangles = np.asarray(triangles)
    lookup = {((i, 4),): i for i in range(len(points))}
    coordinates = list(points)
    elements = np.empty((len(triangles), 15), dtype=int)
    for i, corners in enumerate(triangles):
        for j, alpha in enumerate(indices(3)):
            identity = key(corners, alpha)
            if identity not in lookup:
                lookup[identity] = len(coordinates)
                coordinates.append(alpha@points[corners]/4)
            elements[i, j] = lookup[identity]
    def edge_nodes(edges):
        return np.array([[lookup[key(edge, alpha)] for alpha in indices(2)] for edge in edges])
    top = edge_nodes(np.column_stack((ids[:-1, -1], ids[1:, -1])))
    left = edge_nodes(np.column_stack((ids[0, :-1], ids[0, 1:])))
    right = edge_nodes(np.column_stack((ids[-1, :-1], ids[-1, 1:])))
    return Strip(np.asarray(coordinates), triangles, elements, top, left, right,
                 x, z, omega, angle, k, ky, kappa)


def volume_operator(strip):
    """Exact polynomial integral grad(v).grad(u)+ky^2*v*u on affine triangles."""
    u, w = _gauss(5)
    a, b = np.meshgrid(u, u, indexing='ij');wa, wb = np.meshgrid(w, w, indexing='ij')
    bary = np.column_stack(((1-a-(1-a)*b).ravel(), a.ravel(), ((1-a)*b).ravel()))
    value, derivative = basis(bary)
    derivative = derivative[:, :, 1:]-derivative[:, :, :1]
    weight = (wa*wb*(1-a)).ravel()
    corner = strip.points[strip.triangles]
    jacobian = corner[:, 1:]-corner[:, :1]
    determinant = np.linalg.det(jacobian)
    if np.any(determinant <= 0):
        raise ValueError('strip triangles require positive orientation')
    gradient = np.einsum('qaj,tij->tqai', derivative, np.linalg.inv(jacobian))
    local = np.einsum('tqai,tqbi,q,t->tab', gradient, gradient, weight, determinant)
    local += strip.ky**2*np.einsum('qa,qb,q,t->tab', value, value, weight, determinant)
    row = np.repeat(strip.elements, 15, axis=1).ravel()
    col = np.tile(strip.elements, (1, 15)).ravel()
    return coo_matrix((local.ravel(), (row, col)), shape=(len(strip.points),)*2).tocsc()


def surface_matrices(strip, strength, *, order=12):
    """M exact degree8; Q integrates physical nonnegative C2 sigma, not nodal sigma."""
    if not np.isfinite(strength) or strength < 0 or order < 7:
        raise ValueError('nonnegative finite sponge and quadrature>=7 required')
    free = np.unique(strip.top)
    rowmap = np.full(len(strip.points), -1);rowmap[free] = np.arange(len(free))
    local_ids = rowmap[strip.top]
    length = np.diff(strip.xlevels)
    _, w5, n5 = _line(5)
    mass_local = length[:, None, None]*np.einsum('q,qi,qj->ij', w5, n5, n5)
    t, w, shape = _line(order)
    sample_x = strip.xlevels[:-1, None]+length[:, None]*t
    s = np.clip(sample_x/strip.xlevels[-1], 0., 1.)
    sigma = strength*s**3*(10-15*s+6*s*s)
    weighted_local = length[:, None, None]*np.einsum('q,qi,qj,fq->fij', w, shape, shape, sigma)
    row = np.repeat(local_ids, 5, axis=1).ravel();col = np.tile(local_ids, (1, 5)).ravel()
    mass = coo_matrix((mass_local.ravel(), (row, col)), shape=(len(free),)*2).tocsc()
    weighted = coo_matrix((weighted_local.ravel(), (row, col)), shape=mass.shape).tocsc()
    return free, mass, weighted


def vertical_modes(omega, count):
    """Finite-depth Robin eigenmodes; count is the EVANESCENT count."""
    if count < 1 or count > 512:
        raise ValueError('use1..512 evanescent port modes')
    beta = omega**2/GRAVITY
    k = brentq(lambda v: v*np.tanh(DEPTH*v)-beta, 1e-9, max(10., 2*beta))
    roots = np.array([brentq(lambda v: -v*np.tan(DEPTH*v)-beta,
        ((m+.5)*np.pi+1e-9)/DEPTH, ((m+1)*np.pi-1e-9)/DEPTH)
        for m in range(count)])
    norm = np.r_[DEPTH/(2*np.cosh(k*DEPTH)**2)+np.tanh(k*DEPTH)/(2*k),
                 DEPTH/2+np.sin(2*roots*DEPTH)/(4*roots)]
    return k, roots, norm


def port(strip, side, count=64, *, quadrature_multiplier=1):
    """Continuum modal DtN projected onto the FE trace, with resolved integration.

    The formula alone is not an exact finite-dimensional outgoing port: modal
    truncation and projection must be independently converged. Quadrature order
    grows with the HIGHEST mode's oscillations on EACH vertical segment.
    """
    if side not in ('left', 'right') or quadrature_multiplier not in (1, 2):
        raise ValueError('left/right port and quadrature multiplier1/2 required')
    edges = getattr(strip, side);nodes = np.unique(edges)
    lookup = np.full(len(strip.points), -1);lookup[nodes] = np.arange(len(nodes))
    k, mu, norm = vertical_modes(strip.omega, count)
    # Propagating mode is column0; evanescent root mu[0] is column1.
    projection = np.zeros((len(nodes), count+1))
    orders = []
    for edge in edges:
        # edges' P4 index ordering is alpha=(0,4)..(4,0), not endpoint-first.
        z = strip.points[edge, 1];low, high = z.min(), z.max();length = high-low
        order = int(quadrature_multiplier*(np.ceil(mu[-1]*length/2)+12))
        t, w, shape = _line(order);sample = low+length*t
        vertical = np.column_stack((np.cosh(k*(sample+DEPTH))/np.cosh(k*DEPTH),
                                     np.cos((sample[:, None]+DEPTH)*mu)))
        local = length*np.einsum('q,qi,qm->im', w, shape, vertical)
        # Original edges are oriented low->high; _line's alpha ordering agrees.
        np.add.at(projection, lookup[edge], local)
        orders.append(order)
    eigenvalue = np.r_[1j*strip.kappa, -np.sqrt(mu*mu+strip.ky**2)]
    matrix = (projection*(eigenvalue/norm))@projection.T
    return dict(nodes=nodes, projection=projection, matrix=matrix, norm=norm,
                k=k, mu=mu, quadrature_orders=orders, evanescent_modes=count)


def _scatter_port(value, size):
    ids = value['nodes'];n = len(ids)
    return coo_matrix((value['matrix'].ravel(),
        (np.repeat(ids, n), np.tile(ids, n))), shape=(size, size)).tocsc()


def depth_mode_amplitude(strip, potential, stations, *, order=12):
    """Independent sampling/projection onto the physical propagating depth mode."""
    t, w = _gauss(order)
    norm = DEPTH/(2*np.cosh(strip.k*DEPTH)**2)+np.tanh(strip.k*DEPTH)/(2*strip.k)
    answer = []
    for x in stations:
        i = np.searchsorted(strip.xlevels, x, side='right')-1
        if not 0 <= i < len(strip.xlevels)-1:
            raise ValueError('modal extraction stations must be strictly inside the strip')
        u = (x-strip.xlevels[i])/(strip.xlevels[i+1]-strip.xlevels[i])
        z, weight, layer = [], [], []
        # Every arbitrary vertical cut crosses the diagonal inside each brick.
        # Split there: an unsplit quadrature would integrate across a kink in
        # the FE derivative/trace polynomial and contaminate reflected modes.
        for j, (low, high) in enumerate(zip(strip.zlevels[:-1], strip.zlevels[1:])):
            crossing = low+u*(high-low)
            for a, b in ((low, crossing), (crossing, high)):
                if b > a:
                    z.extend(a+(b-a)*t);weight.extend((b-a)*w);layer.extend([j]*len(t))
        z, weight, j = np.asarray(z), np.asarray(weight), np.asarray(layer)
        v = (z-strip.zlevels[j])/np.diff(strip.zlevels)[j]
        lower = v <= u
        bary = np.where(lower[:, None], np.column_stack((np.full(len(z), 1-u), u-v, v)),
                        np.column_stack((1-v, np.full(len(z), u), v-u)))
        element = 2*(i*(len(strip.zlevels)-1)+j)+(~lower)
        shape, _ = basis(bary)
        sampled = np.einsum('qi,qi->q', shape, potential[strip.elements[element]])
        zmode = np.cosh(strip.k*(z+DEPTH))/np.cosh(strip.k*DEPTH)
        answer.append(np.sum(weight*zmode*sampled)/norm)
    return np.asarray(answer)



def solve_case(strip, strength=4., *, right='wall', modes=64,
               port_quadrature_multiplier=1, surface_order=12):
    """One small strip solve; source review precedes any execution by design."""
    if right not in ('wall', 'outgoing'):
        raise ValueError('right end must be wall or outgoing')
    energy = volume_operator(strip)
    free, mass, weighted = surface_matrices(strip, strength, order=surface_order)
    left = port(strip, 'left', modes, quadrature_multiplier=port_quadrature_multiplier)
    energy = energy.astype(complex)-_scatter_port(left, len(strip.points))
    ports = {'left': left}
    if right == 'outgoing':
        ports['right'] = port(strip, 'right', modes, quadrature_multiplier=port_quadrature_multiplier)
        energy -= _scatter_port(ports['right'], len(strip.points))
    derivative = (weighted-1j*strip.omega*mass).tocoo()
    coupling = coo_matrix((-derivative.data, (free[derivative.row], derivative.col)),
                         shape=(len(strip.points), len(free))).tocsc()
    system = bmat([[energy, coupling], [coupling.T, -GRAVITY*mass]], format='csc')
    input_bytes = system.data.nbytes+system.indices.nbytes+system.indptr.nbytes
    if system.shape[0] > 30000 or system.nnz > 2000000 or input_bytes > 60000000:
        raise RuntimeError('strip input preflight exceeded; these caps do not bound LU fill')
    source = np.zeros(system.shape[0], complex)
    source[left['nodes']] = -2j*strip.kappa*left['projection'][:, 0]
    factor = splu(system);solution = factor.solve(source)
    potential, eta = solution[:len(strip.points)], solution[len(strip.points):]
    residual = float(np.linalg.norm(system@solution-source)/np.linalg.norm(source))
    left_trace = left['projection'][:, 0]@potential[left['nodes']]/left['norm'][0]
    port_reflection = left_trace-1
    transmission = 0j
    if right == 'outgoing':
        rp = ports['right'];transmission = rp['projection'][:, 0]@potential[rp['nodes']]/rp['norm'][0]
    dissipation = float((2*strip.omega/GRAVITY*np.vdot(potential[free], weighted@potential[free])).real)
    incoming_flux = strip.kappa*left['norm'][0]
    balance = 1-abs(port_reflection)**2-abs(transmission)**2-dissipation/incoming_flux
    # Fixed stations and a narrower alternate window expose contamination and
    # phase/mesh errors. No correction or subtraction of control reflection.
    fits = []
    extraction_errors = []
    quiet = -strip.xlevels[0]
    for low, high in ((-.85, -.15), (-.7, -.3)):
        station = quiet*np.linspace(low, high, 17)
        measured = depth_mode_amplitude(strip, potential, station)
        checked = depth_mode_amplitude(strip, potential, station, order=16)
        extraction_errors.append(float(np.linalg.norm(measured-checked)/np.linalg.norm(checked)))
        design = np.column_stack((np.exp(1j*strip.kappa*station), np.exp(-1j*strip.kappa*station)))
        coefficient = np.linalg.lstsq(design, measured, rcond=None)[0]
        scale = np.linalg.norm(coefficient)
        if abs(coefficient[0]) < 1e-8*scale:
            raise FloatingPointError('reflection fit has unresolved incident amplitude')
        reflection = coefficient[1]/coefficient[0]
        fits.append(dict(reflection=[float(reflection.real), float(reflection.imag)],
            reflection_amplitude=float(abs(reflection)), reflected_power=float(abs(reflection)**2),
            incident_amplitude=float(abs(coefficient[0])), design_condition=float(np.linalg.cond(design)),
            relative_fit_residual=float(np.linalg.norm(design@coefficient-measured)/np.linalg.norm(measured))))
    reference_reflection = port_reflection*np.exp(2j*strip.kappa*strip.xlevels[0])
    expected_transmission = np.exp(1j*strip.kappa*(strip.xlevels[-1]-strip.xlevels[0]))
    return dict(qualified=False, omega=strip.omega, angle_degrees=float(np.rad2deg(strip.angle)),
        width=float(strip.xlevels[-1]), quiet_length=quiet, strength=strength, right_boundary=right,
        zlevels=strip.zlevels.tolist(), maximum_horizontal_spacing=float(np.diff(strip.xlevels).max()),
        potential_nodes=len(strip.points), surface_nodes=len(free), evanescent_port_modes=modes,
        port_quadrature_multiplier=port_quadrature_multiplier,
        port_quadrature_order_range=[min(left['quadrature_orders']), max(left['quadrature_orders'])],
        surface_quadrature_order=surface_order, complex_system_rows=system.shape[0],
        complex_input_nonzeros=system.nnz, complex_input_array_bytes=int(input_bytes),
        factor_memory_bound=False, relative_system_residual=residual,
        port_reflection_at_sponge_start=[float(reference_reflection.real), float(reference_reflection.imag)],
        port_reflection_amplitude=float(abs(port_reflection)),
        port_transmission=[float(transmission.real), float(transmission.imag)],
        transmission_amplitude=float(abs(transmission)),
        transmission_phase_error_radians=float(np.angle(transmission/expected_transmission)) if right=='outgoing' else None,
        relative_flux_balance_error=float(abs(balance)), dissipated_incoming_power_fraction=dissipation/incoming_flux,
        quiet_region_fits=fits,
        depth_mode_q12_q16_maximum_relative_difference=max(extraction_errors),
        scope='Body-free finite-depth strip only. Modal truncation, port integration and spatial convergence '
              'must be checked separately. No reflected-wave or canoe qualification from an algebraic power balance.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--omega', type=float, default=2.)
    parser.add_argument('--angle', type=float, default=0., help='degrees0..45')
    parser.add_argument('--width', type=float, default=1.5943951023931955)
    parser.add_argument('--strength', type=float, default=0.)
    parser.add_argument('--right', choices=['wall', 'outgoing'], default='outgoing')
    parser.add_argument('--spacing', type=float, default=np.pi/6)
    parser.add_argument('--layers', type=int, choices=[4, 8], default=4)
    parser.add_argument('--modes', type=int, default=64)
    parser.add_argument('--port-quadrature-multiplier', type=int, choices=[1, 2], default=1)
    parser.add_argument('--output', type=Path, default=Path('artifacts/sponge-strip-control.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        strip = make_strip(args.omega, np.deg2rad(args.angle), args.width,
                           spacing=args.spacing, layers=args.layers)
        report = solve_case(strip, args.strength, right=args.right, modes=args.modes,
                            port_quadrature_multiplier=args.port_quadrature_multiplier)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report), flush=True)
