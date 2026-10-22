"""Finite-depth matched-coordinate PML strip falsifier, not solver adoption.

The frequency-domain complex coordinate changes the FULL elliptic x-z problem,
its ky^2 mass and its surface measure. This is not the physical Q sponge.
Reuses only the independently checked strip geometry/port/extraction primitives.
Coordinate viewpoint: https://arxiv.org/html/2310.04127v1 (finite-depth waves),
https://arxiv.org/html/2108.05348v1 (PML discretization/angle limitations).
No time-domain nonlinear boundary or canoe accuracy follows from this test.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .p4_tetra import basis
from .sponge_strip import (GRAVITY, DEPTH, make_strip, port, _scatter_port,
                           _gauss, _line, depth_mode_amplitude)


def stretching(x, width, strength, omega, *, real_terminal=1.):
    """Return f, f'=a+i*sigma/omega; the default retains the old coordinate."""
    if (not all(np.isfinite(v) for v in (width, strength, omega, real_terminal))
            or width <= 0 or strength < 0 or omega <= 0 or real_terminal < 1):
        raise ValueError('positive width/frequency, nonnegative sigma, real_terminal>=1 required')
    x = np.asarray(x)
    fraction = np.clip(x/width, 0., 1.)
    sigma = strength*fraction**3*(10-15*fraction+6*fraction**2)
    integral = strength*width*(2.5*fraction**4-3*fraction**5+fraction**6)
    # This diagnostic's mesh stops at width. The continuation past it has the
    # terminal constant stretch, so f'=s also holds outside the mesh.
    integral += strength*np.maximum(x-width, 0.)
    if real_terminal == 1.:
        # Preserve the previously checked default arithmetic exactly.
        return x+1j*integral/omega, 1+1j*sigma/omega, sigma
    real_integral = width*(2.5*fraction**4-3*fraction**5+fraction**6)+np.maximum(x-width, 0.)
    real_coordinate = x+(real_terminal-1)*real_integral
    real_derivative = 1+(real_terminal-1)*fraction**3*(10-15*fraction+6*fraction**2)
    return real_coordinate+1j*integral/omega, real_derivative+1j*sigma/omega, sigma


def exact_field(strip, points, strength, *, right='wall', real_terminal=1.):
    """Analytic full field/derivatives with UNIT incoming mode at the left port."""
    if right not in ('wall', 'outgoing'):
        raise ValueError('right end must be wall or outgoing')
    points = np.asarray(points)
    x, z = points[..., 0], points[..., 1]
    width = strip.xlevels[-1];left = strip.xlevels[0]
    coordinate, stretch, _ = stretching(x, width, strength, strip.omega, real_terminal=real_terminal)
    end, _, _ = stretching(width, width, strength, strip.omega, real_terminal=real_terminal)
    forward = np.exp(1j*strip.kappa*(coordinate-left))
    backward = (np.exp(1j*strip.kappa*(2*end-coordinate-left))
                if right == 'wall' else np.zeros_like(forward))
    vertical = np.cosh(strip.k*(z+DEPTH))/np.cosh(strip.k*DEPTH)
    vertical_derivative = strip.k*np.sinh(strip.k*(z+DEPTH))/np.cosh(strip.k*DEPTH)
    potential = vertical*(forward+backward)
    derivative = np.stack((1j*strip.kappa*stretch*vertical*(forward-backward),
                           vertical_derivative*(forward+backward)), axis=-1)
    reflection = np.exp(2j*strip.kappa*end) if right == 'wall' else 0j
    return potential, derivative, complex(reflection)


def _triangle_rule(order):
    u, w = _gauss(order)
    a, b = np.meshgrid(u, u, indexing='ij');wa, wb = np.meshgrid(w, w, indexing='ij')
    bary = np.column_stack(((1-a-(1-a)*b).ravel(), a.ravel(), ((1-a)*b).ravel()))
    value, derivative = basis(bary)
    return value, derivative[:, :, 1:]-derivative[:, :, :1], (wa*wb*(1-a)).ravel()


def pml_matrices(strip, strength, *, volume_order=12, surface_order=12, real_terminal=1.):
    """Current physical mesh, complex coordinate coefficients; no physical Q.

    1/s is rational, so polynomial order5 exactness does NOT certify this
    volume integral. Report/control12 versus18 independently of mesh changes.
    """
    value, derivative, weight = _triangle_rule(volume_order)
    corner = strip.points[strip.triangles]
    jacobian = corner[:, 1:]-corner[:, :1]
    determinant = np.linalg.det(jacobian)
    if np.any(determinant <= 0):
        raise ValueError('positive affine strip triangles required')
    gradient = np.einsum('qaj,tij->tqai', derivative, np.linalg.inv(jacobian))
    sample = np.einsum('qa,tai->tqi', value, strip.points[strip.elements])
    _, stretch, _ = stretching(sample[..., 0], strip.xlevels[-1], strength, strip.omega, real_terminal=real_terminal)
    local = np.einsum('tqa,tqb,tq,q,t->tab', gradient[..., 0], gradient[..., 0],
                      1/stretch, weight, determinant)
    local += np.einsum('tqa,tqb,tq,q,t->tab', gradient[..., 1], gradient[..., 1],
                       stretch, weight, determinant)
    local += strip.ky**2*np.einsum('qa,qb,tq,q,t->tab', value, value, stretch, weight, determinant)
    row = np.repeat(strip.elements, 15, axis=1).ravel();col = np.tile(strip.elements, (1, 15)).ravel()
    energy = coo_matrix((local.ravel(), (row, col)), shape=(len(strip.points),)*2).tocsc()
    free = np.unique(strip.top)
    lookup = np.full(len(strip.points), -1);lookup[free] = np.arange(len(free))
    ids = lookup[strip.top];length = np.diff(strip.xlevels)
    t, weight, value = _line(surface_order)
    sample_x = strip.xlevels[:-1, None]+length[:, None]*t
    _, surface_stretch, _ = stretching(sample_x, strip.xlevels[-1], strength, strip.omega, real_terminal=real_terminal)
    local = length[:, None, None]*np.einsum('q,qi,qj,tq->tij', weight, value, value, surface_stretch)
    row = np.repeat(ids, 5, axis=1).ravel();col = np.tile(ids, (1, 5)).ravel()
    surface_mass = coo_matrix((local.ravel(), (row, col)), shape=(len(free),)*2).tocsc()
    return energy, free, surface_mass


def field_check(strip, potential, strength, *, right='wall', order=18, real_terminal=1.):
    """Independent full analytic field error and transformed-form power.

    The gradient norm inside PML is a COMPUTATIONAL derivative error, not a
    physical water velocity/energy. Quiet-region gradients ARE physical, with
    ky*u included. No positive-Q physical-sponge identity is used.
    """
    value, derivative, weight = _triangle_rule(order)
    error = np.zeros(4);reference = np.zeros(4);imaginary_form = 0.
    for start in range(0, len(strip.elements), 16):
        nodes = strip.elements[start:start+16]
        corner = strip.points[strip.triangles[start:start+16]]
        jacobian = corner[:, 1:]-corner[:, :1]
        determinant = np.linalg.det(jacobian)
        gradient = np.einsum('qaj,tij->tqai', derivative, np.linalg.inv(jacobian))
        points = np.einsum('qa,tai->tqi', value, strip.points[nodes])
        actual = np.einsum('qa,ta->tq', value, potential[nodes])
        actual_derivative = np.einsum('tqai,ta->tqi', gradient, potential[nodes])
        expected, expected_derivative, _ = exact_field(strip, points, strength, right=right, real_terminal=real_terminal)
        measure = determinant[:, None]*weight
        e0 = abs(actual-expected)**2;r0 = abs(expected)**2
        e1 = np.sum(abs(actual_derivative-expected_derivative)**2, axis=-1)+strip.ky**2*e0
        r1 = np.sum(abs(expected_derivative)**2, axis=-1)+strip.ky**2*r0
        quiet = points[..., 0] < 0
        for i, (e, r) in enumerate(((e0, r0), (e1, r1), (e0*quiet, r0*quiet), (e1*quiet, r1*quiet))):
            error[i] += np.sum(measure*e);reference[i] += np.sum(measure*r)
        _, stretch, _ = stretching(points[..., 0], strip.xlevels[-1], strength, strip.omega, real_terminal=real_terminal)
        density = (1/stretch).imag*abs(actual_derivative[..., 0])**2
        density += stretch.imag*(abs(actual_derivative[..., 1])**2+strip.ky**2*abs(actual)**2)
        imaginary_form += np.sum(measure*density)
    # Independent physical-coordinate FS integral of the transformed Robin term.
    t, w, shape = _line(20)
    length = np.diff(strip.xlevels)
    x = strip.xlevels[:-1, None]+length[:, None]*t
    trace = np.einsum('qa,ta->tq', shape, potential[strip.top])
    _, stretch, _ = stretching(x, strip.xlevels[-1], strength, strip.omega, real_terminal=real_terminal)
    imaginary_form -= strip.omega**2/GRAVITY*np.einsum('t,q,tq->', length, w, stretch.imag*abs(trace)**2)
    names = ('full_strip_potential_relative_l2', 'full_strip_computational_gradient_relative_l2',
             'quiet_potential_relative_l2', 'quiet_physical_velocity_relative_l2')
    return dict(zip(names, np.sqrt(error/reference).tolist())), float(imaginary_form)


def solve_pml_case(strip, strength=4., *, right='wall', modes=128,
                   volume_order=12, surface_order=12, formulation='augmented',
                   port_quadrature_multiplier=1, real_terminal=1.):
    """One reviewed small frequency-domain calculation; no family auto-run."""
    if formulation not in ('augmented', 'scalar') or right not in ('wall', 'outgoing'):
        raise ValueError('augmented/scalar and wall/outgoing required')
    energy, free, mass = pml_matrices(strip, strength, volume_order=volume_order,
                                     surface_order=surface_order, real_terminal=real_terminal)
    # Scalar equivalent before ports, used only for the transformed-form balance.
    mass_entries = mass.tocoo()
    robin = coo_matrix((mass_entries.data, (free[mass_entries.row], free[mass_entries.col])),
                       shape=energy.shape).tocsc()
    scalar = energy-strip.omega**2/GRAVITY*robin
    left = port(strip, 'left', modes, quadrature_multiplier=port_quadrature_multiplier)
    port_operator = _scatter_port(left, len(strip.points))
    rp = None
    if right == 'outgoing':
        rp = port(strip, 'right', modes, quadrature_multiplier=port_quadrature_multiplier)
        port_operator += _scatter_port(rp, len(strip.points))
    if formulation == 'augmented':
        derivative = (1j*strip.omega*mass).tocoo()
        coupling = coo_matrix((derivative.data, (free[derivative.row], derivative.col)),
                             shape=(len(strip.points), len(free))).tocsc()
        system = bmat([[energy-port_operator, coupling], [coupling.T, -GRAVITY*mass]], format='csc')
    else:
        system = scalar-port_operator
    input_bytes = int(system.data.nbytes+system.indices.nbytes+system.indptr.nbytes)
    if system.shape[0] > 30000 or system.nnz > 2000000 or input_bytes > 60000000:
        raise RuntimeError('strip input preflight exceeded; not a factor-memory bound')
    source = np.zeros(system.shape[0], complex)
    source[left['nodes']] = -2j*strip.kappa*left['projection'][:, 0]
    factor = splu(system);solution = factor.solve(source)
    if not np.all(np.isfinite(solution)):
        raise FloatingPointError('nonfinite matched-coordinate strip solution')
    potential = solution[:len(strip.points)]
    eta = solution[len(strip.points):] if formulation=='augmented' else 1j*strip.omega/GRAVITY*potential[free]
    eta_error = float(np.linalg.norm(eta-1j*strip.omega/GRAVITY*potential[free])/np.linalg.norm(eta))
    residual = float(np.linalg.norm(system@solution-source)/np.linalg.norm(source))
    left_trace = left['projection'][:, 0]@potential[left['nodes']]/left['norm'][0]
    reflection = (left_trace-1)*np.exp(2j*strip.kappa*strip.xlevels[0])
    transmission = rp['projection'][:, 0]@potential[rp['nodes']]/rp['norm'][0] if rp else 0j
    _, _, exact_reflection = exact_field(strip, strip.points[:1], strength, right=right, real_terminal=real_terminal)
    incoming_flux = strip.kappa*left['norm'][0]
    discrete_absorption = -float(np.vdot(potential, scalar@potential).imag)/incoming_flux
    physical_errors, independently_integrated_imaginary_form = field_check(strip, potential, strength, right=right, real_terminal=real_terminal)
    independent_absorption = -independently_integrated_imaginary_form/incoming_flux
    fits = []
    for lo, hi in ((-.85, -.15), (-.7, -.3)):
        station = -strip.xlevels[0]*np.linspace(lo, hi, 17)
        measured = depth_mode_amplitude(strip, potential, station)
        checked = depth_mode_amplitude(strip, potential, station, order=16)
        design = np.column_stack((np.exp(1j*strip.kappa*station), np.exp(-1j*strip.kappa*station)))
        coefficient = np.linalg.lstsq(design, measured, rcond=None)[0]
        if (not np.all(np.isfinite(coefficient))
                or abs(coefficient[0]) < 1e-8*max(1., np.linalg.norm(coefficient))):
            raise FloatingPointError('unresolved unit-port incident amplitude')
        fitted = coefficient[1]/coefficient[0]
        fits.append(dict(reflection=[float(fitted.real), float(fitted.imag)],
            incident_coefficient=[float(coefficient[0].real), float(coefficient[0].imag)],
            reflected_coefficient=[float(coefficient[1].real), float(coefficient[1].imag)],
            absolute_complex_reflection_error=float(abs(fitted-exact_reflection)),
            relative_fit_residual=float(np.linalg.norm(design@coefficient-measured)/np.linalg.norm(measured)),
            depth_mode_q12_q16_relative_difference=float(np.linalg.norm(measured-checked)/np.linalg.norm(checked))))
    end, _, _ = stretching(strip.xlevels[-1], strip.xlevels[-1], strength, strip.omega, real_terminal=real_terminal)
    expected_transmission = np.exp(1j*strip.kappa*(end-strip.xlevels[0])) if right=='outgoing' else 0j
    return dict(qualified=False, case='finite-depth complex-coordinate PML strip only',
        omega=strip.omega, angle_degrees=float(np.rad2deg(strip.angle)), width=float(strip.xlevels[-1]),
        strength=strength, real_terminal=real_terminal, right_boundary=right, formulation=formulation,
        layer_coordinate_endpoint=[float(end.real),float(end.imag)],
        represented_real_layer_length=float(end.real),
        zlevels=strip.zlevels.tolist(), maximum_horizontal_spacing=float(np.diff(strip.xlevels).max()),
        potential_nodes=len(strip.points), surface_nodes=len(free), evanescent_modes=modes,
        layer_horizontal_cell_widths=np.diff(strip.xlevels[strip.xlevels>=0]).tolist(),
        port_quadrature_multiplier=port_quadrature_multiplier, volume_quadrature_order=volume_order,
        surface_quadrature_order=surface_order, field_check_volume_order=18, field_check_surface_order=20,
        complex_system_rows=system.shape[0], complex_input_nonzeros=system.nnz,
        complex_input_array_bytes=input_bytes, factor_memory_bound=False,
        relative_system_residual=residual, dynamic_eta_identity_relative_error=eta_error,
        reflection_at_layer_start=[float(reflection.real), float(reflection.imag)],
        exact_reflection_at_layer_start=[float(exact_reflection.real), float(exact_reflection.imag)],
        reflection_amplitude=float(abs(reflection)), exact_reflection_amplitude=float(abs(exact_reflection)),
        absolute_complex_reflection_error=float(abs(reflection-exact_reflection)),
        transmission=[float(transmission.real), float(transmission.imag)],
        absolute_complex_transmission_error=float(abs(transmission-expected_transmission)),
        transformed_form_absorbed_incoming_fraction=discrete_absorption,
        independent_transformed_form_absorbed_incoming_fraction=independent_absorption,
        discrete_modal_flux_balance_error=float(abs(1-abs(reflection)**2-abs(transmission)**2-discrete_absorption)),
        independent_modal_flux_balance_error=float(abs(1-abs(reflection)**2-abs(transmission)**2-independent_absorption)),
        fullfield_errors=physical_errors, quiet_fits=fits,
        scope='Analytic full-field/reflection and transformed-form accuracy falsifier only. '
              'This is NOT physicalQ damping and is NOT a nonlinear time-domain boundary. '
              'Coefficient quadrature, spatial, port and angle convergence remain explicit. '
              'No canoe/Capytaine/outgoing architecture or performance qualification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--omega', type=float, default=2.)
    parser.add_argument('--angle', type=float, default=0., help='degrees0..45')
    parser.add_argument('--width', type=float, default=7.877580409572782)
    parser.add_argument('--strength', type=float, default=4.)
    parser.add_argument('--real-terminal', type=float, default=1.)
    parser.add_argument('--spacing', type=float, default=np.pi/6)
    parser.add_argument('--layers', type=int, choices=[4, 8], default=4)
    parser.add_argument('--volume-order', type=int, default=12)
    parser.add_argument('--modes', type=int, default=128)
    parser.add_argument('--right', choices=['wall', 'outgoing'], default='wall')
    parser.add_argument('--formulation', choices=['augmented', 'scalar'], default='augmented')
    parser.add_argument('--output', type=Path, default=Path('artifacts/pml-strip-control.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        strip = make_strip(args.omega, np.deg2rad(args.angle), args.width,
                           spacing=args.spacing, layers=args.layers)
        result = solve_pml_case(strip, args.strength, right=args.right, modes=args.modes,
                                volume_order=args.volume_order, formulation=args.formulation, real_terminal=args.real_terminal)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result), flush=True)
