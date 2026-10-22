"""Contained finite-frequency six-mode control of the whole-box P4 candidate.

The linear P4 candidate has no P1 interface restriction. It is not the current
P1 nonlinear solver's exact linearization, nor a qualified outgoing-wave model.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse import bmat, coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from tests.benchmark_rankine import peak_rss_bytes
from tests.laplace_reference import triangle_samples
from .moving_canoe import make_domain, geometry
from .p4_whole_canoe import make_volume
from .p4_canoe import harmonic_quartic, loads
from . import p4_tetra
from .validate_radiation import coefficient_gates, MODES


def physical_sponge(points, outer, strength=4., buffer=.5):
    """The existing C2 physical ramp, using fixed outer vertex coordinates."""
    inner = np.array([2*np.pi/3, np.pi/3])+buffer
    widths = np.asarray(outer)-inner
    if strength < 0 or buffer < 0 or np.any(widths <= 0):
        raise ValueError('invalid physical sponge settings')
    s = np.clip((np.abs(points[..., :2])-inner)/widths, 0., 1.)
    ramp = s**3*(10+s*(-15+6*s))
    return strength*(ramp[..., 0]+ramp[..., 1]-ramp[..., 0]*ramp[..., 1]), widths


def surface_matrices(mesh, faces, order, *, strength=4., buffer=.5):
    """Exact P4 mass; positive-quadrature physical sigma-weighted mass.

    Sigma is evaluated physically, never interpolated in the sign-changing
    P4 nodal basis. The same full P4 space is used for phi and eta.
    """
    triangles = mesh.points[faces]
    normal = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    area = normal[:, 2]/2
    if np.any(area <= 0) or np.any(triangles[:, :, 2] != 0):
        raise ValueError('the rest free surface must be flat with positive projected triangles')
    nodes = p4_tetra.face_nodes(mesh, faces)
    free = np.unique(nodes)
    lookup = np.full(len(mesh.points), -1, dtype=int)
    lookup[free] = np.arange(len(free))
    local = lookup[nodes]
    rows = np.repeat(local, 15, axis=1).ravel()
    cols = np.tile(local, (1, 15)).ravel()
    barycentric, weights = triangle_samples(np.eye(3)[None], 5)
    basis, _ = p4_tetra.basis(barycentric[0])
    elements = area[:, None, None]*np.einsum('q,qi,qj->ij', weights, basis, basis)
    mass = coo_matrix((elements.ravel(), (rows, cols)), shape=(len(free), len(free))).tocsc()
    barycentric, weights = triangle_samples(np.eye(3)[None], order)
    barycentric = barycentric[0]
    basis, _ = p4_tetra.basis(barycentric)
    samples = barycentric@triangles
    outer = np.max(np.abs(mesh.points[free, :2]), axis=0)
    sigma, widths = physical_sponge(samples, outer, strength, buffer)
    if np.min(sigma) < 0 or np.max(sigma) > strength*(1+1e-14):
        raise RuntimeError('physical sponge lost nonnegative boundedness')
    elements = area[:, None, None]*np.einsum('q,qi,qj,fq->fij', weights, basis, basis, sigma)
    weighted = coo_matrix((elements.ravel(), (rows, cols)), shape=mass.shape).tocsc()
    return free, mass, weighted, widths


def read_reference(path, wet_hash, center, frequencies):
    data = json.loads(path.read_text())
    if (data['wetted_triangles_sha256'] != wet_hash or data['wetted_panels'] != 704
            or data['water_depth'] != 4. or not np.array_equal(data['rotation_center'], center)
            or data['mode_order'] != list(MODES)):
        raise ValueError('Capytaine report does not match current wet geometry/CG/depth/modes')
    if [case['omega'] for case in data['cases']] != list(frequencies):
        raise ValueError('reference frequency family differs')
    matrices = [case['capytaine'] for case in data['cases']]
    for matrix in matrices:
        for name in ('added_mass', 'radiation_damping'):
            value = np.asarray(matrix[name])
            if value.shape != (6, 6) or not np.all(np.isfinite(value)):
                raise ValueError('invalid complete reference matrix')
    return matrices, data['capytaine_green_function_settings']


def solve_frequency(energy, load, active, free, mass, weighted, omega):
    free_active = np.searchsorted(active, free)
    if np.any(free_active >= len(active)) or not np.array_equal(active[free_active], free):
        raise ValueError('inactive free-surface potential node')
    derivative = (weighted-1j*omega*mass).tocoo()
    coupling = coo_matrix((-derivative.data,
                          (free_active[derivative.row], derivative.col)),
                         shape=(len(active), len(free))).tocsc()
    system = bmat([[energy[active][:, active].astype(complex), coupling],
                   [coupling.T, -9.81*mass]], format='csc')
    # Input storage is a measured count, NOT a bound on pivoted LU fill.
    input_bytes = int(system.data.nbytes+system.indices.nbytes+system.indptr.nbytes)
    if len(active)+len(free) > 100000 or system.nnz > 12000000 or input_bytes > 300000000:
        raise RuntimeError('complex system exceeded the declared input preflight caps')
    print(json.dumps(dict(stage='complex input preflight', omega=omega,
                          rows=system.shape[0], nonzeros=system.nnz,
                          input_array_bytes=input_bytes,
                          factor_memory_bound=False)), flush=True)
    right = np.vstack((load[active], np.zeros((len(free), 6)))).astype(complex)
    factor = splu(system)
    solution = factor.solve(right)
    # Capture actual factor-array counts after the authorized factorization.
    lower, upper = factor.L, factor.U
    factor_storage = dict(lower_nonzeros=lower.nnz, upper_nonzeros=upper.nnz,
        exposed_factor_array_bytes=int(sum(value.data.nbytes+value.indices.nbytes+value.indptr.nbytes
                                           for value in (lower, upper))))
    response = load[active].T@solution[:len(active)]
    trace = solution[free_active]
    dissipation = (2*omega/9.81)*(trace.conj().T@weighted@trace)
    residual = float(np.linalg.norm(system@solution-right)/np.linalg.norm(right))
    reciprocity = float(np.linalg.norm(response-response.T)/np.linalg.norm(response))
    passivity = float(np.linalg.norm(response.imag-dissipation)/np.linalg.norm(response.imag))
    if max(residual, reciprocity, passivity) > 1e-8:
        raise RuntimeError(f'finite-frequency algebra failed: {residual}, {reciprocity}, {passivity}')
    report = dict(omega=omega, added_mass=(1000*response.real).tolist(),
        radiation_damping=(1000*omega*response.imag).tolist(),
        relative_system_residual=residual, reciprocity_relative=reciprocity,
        passivity_identity_relative=passivity,
        minimum_damping_eigenvalue=float(np.linalg.eigvalsh(500*omega*(response.imag+response.imag.T)).min()),
        input_array_bytes=input_bytes, factor_storage=factor_storage,
        process_peak_working_set_bytes=int(peak_rss_bytes()))
    del factor, lower, upper, system, solution
    return report


def run(reference_path):
    frequencies = [2., 4., 6.]
    domain = make_domain(12, False, 2)
    mesh, faces, free_faces, body_faces, preflight = make_volume(domain)
    points, _, _ = geometry(domain, domain.rest.position, domain.rest.orientation,
                             np.zeros(len(domain.free_nodes)))
    wet = points[domain.local_mapping][domain.contact.faces[domain.contact.body_faces]][:, ::-1]
    wet_hash = sha256(np.asarray(wet, dtype='<f8').tobytes()).hexdigest()
    reference, settings = read_reference(reference_path, wet_hash, domain.rest.position, frequencies)
    energy = p4_tetra.stiffness(mesh)
    _, patch = loads(mesh, faces, domain.rest.position)
    polynomial, _ = harmonic_quartic(mesh.points, domain.rest.position)
    patch_error = float(np.linalg.norm(energy@polynomial-patch)/np.linalg.norm(patch))
    if patch_error > 1e-9:
        raise RuntimeError('whole-domain harmonic-quartic patch failed')
    load, _ = loads(mesh, body_faces, domain.rest.position)
    active = np.unique(mesh.elements)
    free, mass, q8, widths = surface_matrices(mesh, free_faces, 8)
    free12, mass12, q12, widths12 = surface_matrices(mesh, free_faces, 12)
    if not np.array_equal(free, free12) or np.linalg.norm((mass-mass12).data) != 0 or not np.array_equal(widths, widths12):
        raise RuntimeError('quadrature control changed the physical trace space or mass')
    q_relative = float(np.linalg.norm((q12-q8).data)/np.linalg.norm(q12.data))
    families = {}
    for order, weighted in ((8, q8), (12, q12)):
        cases = []
        for i, omega in enumerate(frequencies):
            case = solve_frequency(energy, load, active, free, mass, weighted, omega)
            case['capytaine'] = reference[i]
            d = np.array([1., 1., 1., 3.5, 3.5, 3.5])
            scale = 1000*3.5**3*d[:, None]*d[None, :]
            case.update(added_mass_nondimensional=(np.asarray(case['added_mass'])/scale).tolist(),
                reference_added_mass_nondimensional=(np.asarray(reference[i]['added_mass'])/scale).tolist(),
                damping_nondimensional=(np.asarray(case['radiation_damping'])/(omega*scale)).tolist(),
                reference_damping_nondimensional=(np.asarray(reference[i]['radiation_damping'])/(omega*scale)).tolist())
            cases.append(case)
            print(json.dumps(dict(stage='frequency completed', quadrature_order=order,
                omega=omega, relative_system_residual=case['relative_system_residual'],
                reciprocity_relative=case['reciprocity_relative'],
                passivity_identity_relative=case['passivity_identity_relative'],
                factor_storage=case['factor_storage'],
                process_peak_working_set_bytes=case['process_peak_working_set_bytes'])), flush=True)
        families[order] = cases
    response_changes = {}
    for key in ('added_mass_nondimensional', 'damping_nondimensional'):
        coarse = np.array([case[key] for case in families[8]])
        fine = np.array([case[key] for case in families[12]])
        response_changes[key] = dict(relative_frobenius=float(np.linalg.norm(fine-coarse)/np.linalg.norm(fine)),
                                      maximum_absolute=float(np.max(abs(fine-coarse))))
    return dict(qualified=False, case='finite-frequency whole-tank tetra-P4 candidate',
        n=12, wet_subdivisions=2, exterior_scale=1., water_depth=4., density=1000., gravity=9.81,
        mode_order=list(MODES), sponge_strength=4., sponge_buffer=.5, sponge_widths=widths.tolist(),
        current_wet_triangles_sha256=wet_hash, rotation_center=domain.rest.position.tolist(),
        active_p4_nodes=len(active), free_surface_p4_nodes=len(free), preflight=preflight,
        harmonic_quartic_patch_relative_error=patch_error,
        weighted_mass_q8_q12_relative_difference=q_relative,
        quadrature_response_changes=response_changes, cases=families[12], q8_cases=families[8],
        coefficient_gates=coefficient_gates(families[12]), reference_report=str(reference_path),
        capytaine_green_function_settings=settings,
        scope='Single finite-box/rest P4 candidate with unconstrained P4 interface. '
              'Not the current P1 nonlinear implementation. Spatial propagation, '
              'waterline fields, reference, finite-domain/absorber convergence and '
              'nonlinear moving geometry remain unqualified. No tuning or performance qualification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=Path('artifacts/p4-current-radiation.json'))
    args = parser.parse_args()
    with threadpool_limits(1):
        report = run(args.reference)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    gates = {name: {key: gate[key] for key in (
        'passed', 'maximum_tolerance_ratio', 'worst_omega',
        'worst_force_mode', 'worst_radiating_mode')}
        for name, gate in report['coefficient_gates'].items()}
    print(json.dumps(dict(qualified=False, output=str(args.output), coefficient_gates=gates,
        weighted_mass_q8_q12_relative_difference=report['weighted_mass_q8_q12_relative_difference'],
        quadrature_response_changes=report['quadrature_response_changes'])), flush=True)
