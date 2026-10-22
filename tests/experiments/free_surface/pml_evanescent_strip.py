"""Lowest-evanescent-mode PML strip accuracy falsifier; no physical flux ratio.

Imaginary coordinate stretching changes the phase of evanescent modes, not
exp(-alpha*x) spatial decay. A finite wall returns exp(-2*alpha*L), regardless
of PML strength. This is a necessary nearfield control, not canoe adoption.
The scalar operator is the exact elimination of the reviewed Ms augmented
system. No physical Q-sponge power identity applies to this control.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .sponge_strip import (DEPTH, GRAVITY, make_strip, vertical_modes, port, _scatter_port)
from .pml_strip import stretching, pml_matrices, _triangle_rule


def make_evanescent_strip(omega=2., angle=0., width=1.5943951023931955, *,
                           spacing=np.pi/6, layers=4):
    """Use half an EVANESCENT decay length as quiet interval, avoiding bad scaling.

    This changes only the independent diagnostic's launch-port location. It
    does not shorten the PML or alter any canoe domain or contact geometry.
    """
    k, roots, _ = vertical_modes(omega, 1)
    ky, kappa = k*np.sin(angle), k*np.cos(angle)
    alpha = np.sqrt(roots[0]**2+ky**2)
    quiet_length = .5/alpha
    return make_strip(omega, angle, width, spacing=spacing, layers=layers,
                      quiet_wavelengths=quiet_length*kappa/(2*np.pi))


def exact_evanescent_field(strip, points, strength, *, right='wall'):
    """Unit incoming UNNORMALIZED Z1=cos(mu0(z+h)) at the left reference plane."""
    if right not in ('wall', 'outgoing'):
        raise ValueError('wall/outgoing right boundary required')
    _, roots, _ = vertical_modes(strip.omega, 1)
    mu = roots[0];alpha = np.sqrt(mu*mu+strip.ky**2)
    points = np.asarray(points);x, z = points[..., 0], points[..., 1]
    width, left = strip.xlevels[-1], strip.xlevels[0]
    f, s, _ = stretching(x, width, strength, strip.omega)
    end, _, _ = stretching(width, width, strength, strip.omega)
    incoming = np.exp(-alpha*(f-left))
    returning = np.exp(-alpha*(2*end-f-left)) if right=='wall' else np.zeros_like(incoming)
    vertical = np.cos(mu*(z+DEPTH))
    value = vertical*(incoming+returning)
    derivative = np.stack((-alpha*s*vertical*(incoming-returning),
        -mu*np.sin(mu*(z+DEPTH))*(incoming+returning)), axis=-1)
    reflection = np.exp(-2*alpha*end) if right=='wall' else 0j
    return value, derivative, complex(reflection)


def _field_errors(strip, potential, strength, right, *, order=18):
    """Independent real(dx dz) L2 field/error; PML derivatives are computational."""
    value, derivative, weights = _triangle_rule(order)
    error = np.zeros(4);reference = np.zeros(4)
    for start in range(0, len(strip.elements), 16):
        ids = strip.elements[start:start+16]
        corners = strip.points[strip.triangles[start:start+16]]
        jacobian = corners[:, 1:]-corners[:, :1]
        measure = np.linalg.det(jacobian)[:, None]*weights
        gradient = np.einsum('qaj,tij->tqai', derivative, np.linalg.inv(jacobian))
        sample = np.einsum('qa,tai->tqi', value, strip.points[ids])
        actual = np.einsum('qa,ta->tq', value, potential[ids])
        actual_gradient = np.einsum('tqai,ta->tqi', gradient, potential[ids])
        expected, expected_gradient, _ = exact_evanescent_field(strip, sample, strength, right=right)
        e0 = abs(actual-expected)**2;r0 = abs(expected)**2
        e1 = np.sum(abs(actual_gradient-expected_gradient)**2, axis=-1)+strip.ky**2*e0
        r1 = np.sum(abs(expected_gradient)**2, axis=-1)+strip.ky**2*r0
        quiet = sample[..., 0] < 0
        for i, (e, r) in enumerate(((e0, r0), (e1, r1), (e0*quiet, r0*quiet), (e1*quiet, r1*quiet))):
            error[i] += np.sum(measure*e);reference[i] += np.sum(measure*r)
    names = ('full_strip_potential_relative_l2', 'full_strip_computational_gradient_relative_l2',
             'quiet_potential_relative_l2', 'quiet_physical_velocity_relative_l2')
    return dict(zip(names, np.sqrt(error/reference).tolist()))


def solve_evanescent_case(strip, strength=4., *, right='wall', modes=128,
                           volume_order=12, surface_order=12, port_quadrature_multiplier=1):
    """One small reviewed solve. A single evanescent mode has ZERO incident flux."""
    if right not in ('wall', 'outgoing'):
        raise ValueError('wall/outgoing right boundary required')
    energy, free, mass = pml_matrices(strip, strength, volume_order=volume_order,
                                     surface_order=surface_order)
    entries = mass.tocoo()
    robin = coo_matrix((entries.data, (free[entries.row], free[entries.col])), shape=energy.shape).tocsc()
    system = energy-strip.omega**2/GRAVITY*robin
    left = port(strip, 'left', modes, quadrature_multiplier=port_quadrature_multiplier)
    system -= _scatter_port(left, len(strip.points))
    if right=='outgoing':
        system -= _scatter_port(port(strip, 'right', modes,
            quadrature_multiplier=port_quadrature_multiplier), len(strip.points))
    size = int(system.data.nbytes+system.indices.nbytes+system.indptr.nbytes)
    if system.shape[0]>30000 or system.nnz>2000000 or size>60000000:
        raise RuntimeError('strip input caps exceeded; not a factor-memory bound')
    mu = left['mu'][0];alpha = np.sqrt(mu*mu+strip.ky**2)
    source = np.zeros(len(strip.points), complex)
    # Propagating mode is column0. LOWEST EVANESCENT ROOT mu[0] is column1.
    # Incoming e^-alpha(x-xleft) gives (partial_n-Lambda)u=+2alpha*Z1.
    source[left['nodes']] = 2*alpha*left['projection'][:, 1]
    factor = splu(system);potential = factor.solve(source)
    if not np.all(np.isfinite(potential)):
        raise FloatingPointError('nonfinite evanescent-strip field')
    residual = float(np.linalg.norm(system@potential-source)/np.linalg.norm(source))
    trace = left['projection'].T@potential[left['nodes']]/left['norm']
    amplification = np.exp(-2*alpha*strip.xlevels[0])
    reflected_left = trace[1]-1
    reflection = reflected_left*amplification
    _, _, exact = exact_evanescent_field(strip, strip.points[:1], strength, right=right)
    norm = left['norm'][1]
    return dict(qualified=False, case='lowest evanescent finite-depth PML strip only',
        omega=strip.omega, angle_degrees=float(np.rad2deg(strip.angle)), width=float(strip.xlevels[-1]),
        strength=strength, right_boundary=right, mu0=float(mu), alpha=float(alpha),
        selected_vertical_mode='port column1=cos(mu[0]*(z+h)); column0 remains propagating',
        quiet_length=float(-strip.xlevels[0]), left_to_layer_reflection_amplification=float(amplification),
        exact_physical_spatial_roundtrip_decay=float(np.exp(-2*alpha*strip.xlevels[-1])),
        reflection_at_layer_start=[float(reflection.real), float(reflection.imag)],
        exact_reflection_at_layer_start=[float(exact.real), float(exact.imag)],
        reflection_amplitude=float(abs(reflection)), exact_reflection_amplitude=float(abs(exact)),
        absolute_complex_reflection_error=float(abs(reflection-exact)),
        selected_left_modal_trace=[float(trace[1].real), float(trace[1].imag)],
        reflected_coefficient_at_left=[float(reflected_left.real), float(reflected_left.imag)],
        exact_selected_left_modal_trace=[float((1+exact/amplification).real),float((1+exact/amplification).imag)],
        spurious_propagating_mode_amplitude_at_left=float(abs(trace[0])),
        selected_vertical_mode_l2_norm_squared=float(norm),
        fullfield_errors=_field_errors(strip, potential, strength, right),
        relative_system_residual=residual, potential_nodes=len(strip.points),
        complex_input_nonzeros=system.nnz, complex_input_array_bytes=size, factor_memory_bound=False,
        evanescent_port_modes=modes, volume_quadrature_order=volume_order,
        surface_quadrature_order=surface_order, field_check_quadrature_order=18,
        port_quadrature_multiplier=port_quadrature_multiplier,
        maximum_horizontal_spacing=float(np.diff(strip.xlevels).max()), zlevels=strip.zlevels.tolist(),
        incident_physical_power_flux=0., incident_power_normalization_used=False,
        scope='Unit-amplitude evanescent field/reflection only; no normalization by zero physical incident '
              'power. Imaginary stretching does not improve exp(-2alpha L) amplitude decay. '
              'No canoe, nonlinear boundary, nearfield or outgoing qualification from this necessary test.')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--omega',type=float,default=2.)
    parser.add_argument('--angle',type=float,default=0.)
    parser.add_argument('--width',type=float,default=1.5943951023931955)
    parser.add_argument('--strength',type=float,choices=[0.,4.,8.],default=4.)
    parser.add_argument('--spacing',type=float,default=np.pi/6)
    parser.add_argument('--layers',type=int,choices=[4,8],default=4)
    parser.add_argument('--volume-order',type=int,default=12)
    parser.add_argument('--modes',type=int,default=128)
    parser.add_argument('--right',choices=['wall','outgoing'],default='wall')
    parser.add_argument('--output',type=Path,default=Path('artifacts/pml-evanescent-control.json'))
    args=parser.parse_args()
    with threadpool_limits(1):
        strip=make_evanescent_strip(args.omega,np.deg2rad(args.angle),args.width,
            spacing=args.spacing,layers=args.layers)
        result=solve_evanescent_case(strip,args.strength,right=args.right,modes=args.modes,
            volume_order=args.volume_order)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
