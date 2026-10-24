"""Bounded actual fixed-hull Green-energy Ritz verification, not a runtime."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import eigh
from topwater.constants import GRAVITY
from topwater.rankine import TrianglePanels, single_layer_trace
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import RegularizedSnapshot, owned
from .physical_weak_flux import PhysicalWeakFlux
from .cut_hat_geometry import cut_hat
from .compact_p1_maps import body_maps
from .owned_odd_body_trace import owned_odd_maps
from .odd_body_reference import source_flux
from .odd_variational_energy import OddRitzEnergy
from .panel_gradient_reference import checked_influence


def record(value):
    return dict(real=float(np.real(value)), imag=float(np.imag(value)))


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    old, gamma = owned(args.body16), owned(args.gamma)
    for name in ('water', 'phase', 'mode', 'body16'):
        digest = sha256(getattr(args, name).read_bytes()).hexdigest()
        if gamma['parents'][name] != digest or (name != 'body16' and old['parents'][name] != digest):
            raise ValueError('Ritz/raw Gram parent identity differs')
    for name in ('plane_p1_field.py', 'plane_single_layer.py', 'cut_hat_geometry.py'):
        relative = Path('tests/experiments/free_surface')/name
        if sha256(relative.read_bytes()).hexdigest() != snapshot.water['producer_sha256'].get(str(relative)):
            raise ValueError('physical contact source owner changed: '+name)
    patch = snapshot.patch
    with np.load(args.phase.with_suffix('.npz')) as z:
        footprint = z['polygon'].copy()
    with np.load(args.body16.with_suffix('.npz')) as z:
        psi_old, sigma_old = z['psi'].copy(), z['sigma'].copy()
        gram_old = z['q16_gram'].copy()
        cached_h, cached_n = z['q16_H'].copy(), z['q16_NH'].copy()
        cached_r, cached_j = z['q16_R'].copy(), z['q16_JR'].copy()
    with np.load(args.gamma.with_suffix('.npz')) as z:
        gram_reference = z['updated_gram'].copy()
        contact_reference = z['contact_q8_action'].copy()
    if sha256((patch.mass@psi_old).tobytes()).hexdigest() != old['physical_psi_moment_sha256']:
        raise ValueError('old physical state moments differ')
    panels = TrianglePanels(snapshot.triangles, -1)
    n, nb, ns = len(psi_old), len(panels), len(patch.panels)
    affected = snapshot.flux.rows
    touched = np.full(len(patch.nodes), False)
    touched[patch.active_nodes[affected]] = True
    area_count = int(np.any(touched[patch.faces[patch.parent_id]], axis=1).sum())
    hx, hy = snapshot.flux.spacing
    line_count = sum(len(cut_hat(patch.nodes[node], hx, hy, footprint)['waterline']['segments'])
                     for node in patch.active_nodes[affected])
    nq = nb*16**2
    pair_bound = dict(Kcontact8=(area_count*64+line_count*8)*ns,
        F16=len(patch.cut_triangles)*256*nb, body_PN=nq*ns, body_RJ=nq*2*nb)
    pair_bound['small_controls'] = 2*ns+8*nb
    if sum(pair_bound.values()) > args.maximum_pairs:
        raise RuntimeError('full Ritz source-pair preflight exceeds budget')
    # Covers raw/reduced dense blocks, Cholesky/eigh copies/work and streamed
    # pair/maps. This array estimate excludes Python/JIT/library baseline.
    memory = 64*n*n+32*n*nb+32*nq*nb+16*nq+32*64*ns+32*1024**2
    if memory > args.maximum_bytes:
        raise MemoryError('Ritz dense verification array/workspace preflight exceeds budget')
    report.update(complete=False, qualified=False, active_nodes=n, body_panels=nb,
        source_pair_preflight=pair_bound, source_pairs_so_far=0,
        array_workspace_bound_bytes=memory, quadrature=dict(Kcontact=8, F=16, body=16), stages=[])
    arrays = dict(old_psi=psi_old, old_sigma=sigma_old)
    checkpoint(args, report, arrays)

    # The refactored matrix owner must retain the independently public owned
    # centroid trace, including its source-side jump and reflected chain.
    ids = np.array([0, nb//2])
    cp, cn, _ = owned_odd_maps(panels, panels.centroids[ids], ids, np.full((2, 3), 1/3))
    public = single_layer_trace(panels, ids)
    image = checked_influence(panels.centroids[ids]*[1., 1., -1.], panels)
    expected_p = public.single_layer-image.single_layer
    expected_n = np.einsum('mni,mi->mn', public.gradient_single_layer
        -image.gradient_single_layer*[1., 1., -1.], panels.normals[ids])
    np.testing.assert_allclose(cp, expected_p, rtol=1e-10, atol=1e-13)
    np.testing.assert_allclose(cn, expected_n, rtol=1e-10, atol=1e-13)
    report['owned_matrix_centroid_parity_max'] = float(max(abs(cp-expected_p).max(), abs(cn-expected_n).max()))
    report['prior_intrinsic_self_scalar_controls'] = old['intrinsic_self_scalar_controls']

    started = perf_counter()
    flux = PhysicalWeakFlux(patch, footprint, snapshot.flux, order=8,
        maximum_pairs=pair_bound['Kcontact8'], maximum_bytes=args.maximum_bytes)
    np.testing.assert_array_equal(flux.affected_rows, affected)
    np.testing.assert_allclose(flux.contact_rows@psi_old, contact_reference, rtol=1e-10, atol=1e-12)
    K = np.empty((n, n))
    for first in range(0, n, 32):
        last = min(first+32, n)
        trial = np.zeros((n, last-first))
        trial[np.arange(first, last), np.arange(last-first)] = 1.
        K[:, first:last] = flux.apply(trial)
    np.testing.assert_allclose(K@psi_old, flux.apply(psi_old.real)+1j*flux.apply(psi_old.imag),
                               rtol=1e-10, atol=1e-12)
    pairs = flux.diagnostics['source_pairs']
    report['stages'].append(dict(stage='K8', seconds=perf_counter()-started,
        diagnostics=flux.diagnostics))
    del flux
    F, diagnostics = source_flux(patch, panels, 16, maximum_pairs=pair_bound['F16'],
                                 maximum_bytes=args.maximum_bytes)
    pairs += diagnostics['source_pairs']
    report['stages'].append(dict(stage='F16', **diagnostics))
    arrays.update(K=K, F=F)
    report['source_pairs_so_far'] = pairs
    checkpoint(args, report, arrays)

    points, qw = triangle_samples(panels.triangles, 16)
    points = points.reshape(-1, 3)
    owners = np.repeat(np.arange(nb), len(qw))
    bary, _ = triangle_samples(np.eye(3)[None], 16)
    bary = np.tile(bary[0], (nb, 1))
    weights = (panels.areas[:, None]*qw).ravel()
    normals = panels.normals[owners]
    T = np.column_stack((normals, np.cross(points-snapshot.inputs.initial_state.position, normals)))
    started = perf_counter()
    R, J, diagnostics = owned_odd_maps(panels, points, owners, bary)
    pairs += diagnostics['source_pairs']
    report['stages'].append(dict(stage='R_J16', seconds=perf_counter()-started, diagnostics=diagnostics))
    HH, HR = K.copy(), F.copy()
    RH, RR = np.zeros((nb, n)), R.T@(weights[:, None]*J)
    BH, BR = np.zeros((6, n)), T.T@(weights[:, None]*R)
    hp, hn = np.empty(nq, complex), np.empty(nq, complex)
    started = perf_counter()
    for first in range(0, nq, 64):
        last = min(first+64, nq)
        p, normal = body_maps(patch, points[first:last], normals[first:last])
        pairs += (last-first)*ns
        if pairs > args.maximum_pairs:
            raise RuntimeError('Ritz source-pair budget exhausted')
        p, normal = patch.condense_maps(p), patch.condense_maps(normal)
        w = weights[first:last, None]
        HH += p.T@(w*normal)
        HR += p.T@(w*J[first:last])
        RH += R[first:last].T@(w*normal)
        BH += T[first:last].T@(w*p)
        hp[first:last], hn[first:last] = p@psi_old, normal@psi_old
        if last % 512 == 0 or last == nq:
            print(json.dumps(dict(stage='streamed_H16', completed_points=last,
                total_points=nq, seconds=perf_counter()-started, source_pairs=pairs)), flush=True)
    report['stages'].append(dict(stage='streamed_H16', seconds=perf_counter()-started))
    rp, rn = R@sigma_old, J@sigma_old
    for actual, cached in ((hp, cached_h), (hn, cached_n), (rp, cached_r), (rn, cached_j)):
        np.testing.assert_allclose(actual, cached, rtol=1e-9, atol=1e-11)
    gram = np.array([[np.vdot(psi_old, HH@psi_old), np.vdot(psi_old, HR@sigma_old)],
                     [np.vdot(sigma_old, RH@psi_old), np.vdot(sigma_old, RR@sigma_old)]])
    np.testing.assert_allclose(gram, gram_reference, rtol=1e-8, atol=1e-11)
    arrays.update(HH=HH, HR=HR, RH=RH, RR=RR, BH=BH, BR=BR,
        body_R=R, body_J=J, body_H_old=hp, body_NH_old=hn,
        body_weights=weights, body_T=T, body_points=points, source_areas=panels.areas,
        old_gram=gram, prior_body16_gram=gram_old)
    report.update(source_pairs_so_far=pairs, raw_gram_reference_max=float(abs(gram-gram_reference).max()))
    checkpoint(args, report, arrays)

    energy = OddRitzEnergy(HH, HR, RH, RR, BH, BR, panels.areas)
    report['raw_diagnostics'] = energy.raw_diagnostics
    report['reduction_symmetry_max'] = energy.reduction_symmetry_max
    arrays['source_energy_eigenvalues'] = energy.source_eigenvalues
    # Independent evaluation of the raw physical Green integrand, not the
    # already-contracted symmetric block quadratic. Only the fixed old-mode
    # real/imaginary span is used; all body fields are streamed/cached above.
    def field_functional(a, sigma, u, *, with_scale=False):
        psi = a[0]*psi_old.real+a[1]*psi_old.imag
        H = a[0]*hp.real+a[1]*hp.imag
        N = a[0]*hn.real+a[1]*hn.imag
        source, derivative = R@sigma, J@sigma
        k, f = K@psi, F@sigma
        value = .5*(psi@k+psi@f+weights@((H+source)*(N+derivative)))-u@(T.T@(weights*(H+source)))
        scale = .5*(abs(psi*k).sum()+abs(psi*f).sum()
            +weights@(abs(H+source)*abs(N+derivative)))+abs(u)@(abs(T).T@(weights*abs(H+source)))
        return (value, float(scale)) if with_scale else value

    a, sigma, u = np.array([1., 0.]), sigma_old.real, np.zeros(6)
    psi = psi_old.real
    gp, gs, gu = energy.gradients(psi, sigma, u)
    direction = (np.array([0., 1.]), sigma_old.imag, np.array([.1, -.2, .3, -.04, .05, -.06]))
    expected = gp@psi_old.imag+gs@direction[1]+gu@direction[2]
    directional_scale = abs(gp*psi_old.imag).sum()+abs(gs*direction[1]).sum()+abs(gu*direction[2]).sum()
    # Dimension-dependent binary64 arithmetic allowance; this gates ownership
    # of the quadratic functional, not continuum or physical force accuracy.
    arithmetic_factor = 64*np.finfo(float).eps*(n+nq+nb)
    fd = []
    report['raw_field_FD'] = fd
    for h in (1e-4, 5e-5):
        plus, plus_scale = field_functional(a+h*direction[0], sigma+h*direction[1], u+h*direction[2], with_scale=True)
        minus, minus_scale = field_functional(a-h*direction[0], sigma-h*direction[1], u-h*direction[2], with_scale=True)
        error = float(abs((plus-minus)/(2*h)-expected))
        allowance = arithmetic_factor*(directional_scale+(plus_scale+minus_scale)/(2*h))
        fd.append(dict(h=h, derivative=float((plus-minus)/(2*h)), absolute_error=error,
            absolute_directional_gradient_scale=float(directional_scale), arithmetic_allowance=float(allowance)))
        if error > allowance:
            raise ValueError('raw physical-field functional derivative ownership gate failed')
    # Old physical psi with zero velocity and every unit body velocity: these
    # are algebraic reduction/derivative controls, not radiation validation.
    case_psi = np.repeat(psi_old[:, None], 7, axis=1)
    case_u = np.column_stack((np.zeros(6), np.eye(6)))
    case_sigma = energy.solve_density(case_psi, case_u)
    gp_cases, gs_cases, gu_cases = energy.gradients(case_psi, case_sigma, case_u)
    source_scale = (abs(energy.cross_energy).T@abs(case_psi)
        +abs(energy.source_energy)@abs(case_sigma)+abs(energy.BR).T@abs(case_u)).max()
    wave_scale = (abs(energy.wave_energy)@abs(case_psi)
        +abs(energy.cross_energy)@abs(case_sigma)+abs(energy.BH).T@abs(case_u)).max()
    port_scale = (abs(energy.BH)@abs(case_psi)+abs(energy.BR)@abs(case_sigma)).max()
    algebraic_controls = []
    report['algebraic_reduction_controls'] = algebraic_controls
    for name, error, scale in (
            ('source_stationarity', abs(gs_cases).max(), source_scale),
            ('wave_flux_unreduced_gradient', abs(energy.wave_flux(case_psi, case_u)-gp_cases).max(), wave_scale),
            ('body_port_unreduced_gradient', abs(energy.body_port(case_psi, case_u)+gu_cases).max(), port_scale)):
        allowance = arithmetic_factor*scale
        algebraic_controls.append(dict(name=name, maximum_absolute=float(error), arithmetic_allowance=float(allowance)))
        if error > allowance:
            raise ValueError('algebraic Ritz reduction ownership gate failed: '+name)
    sigma_ritz = energy.solve_density(psi_old, np.zeros(6))
    physical_flux = K@psi_old+F@sigma_ritz
    variational_flux = energy.wave_flux(psi_old, np.zeros(6))
    potential_port = -T.T@(weights*(hp+R@sigma_ritz))
    body_port = energy.body_port(psi_old, np.zeros(6))
    port_error = float(abs(potential_port+body_port).max())
    port_scale = float((abs(T).T@(weights*(abs(hp)+abs(R)@abs(sigma_ritz)))).max())
    port_allowance = arithmetic_factor*port_scale
    report.update(body_potential_port_quadrature_parity_max=port_error,
                  body_potential_port_arithmetic_allowance=float(port_allowance))
    if port_error > port_allowance:
        raise ValueError('quadrature body potential-port ownership gate failed')
    qbody = hn+J@sigma_ritz
    flux_difference = physical_flux-variational_flux
    flux_mass_norm = np.sqrt(np.real(np.vdot(flux_difference, snapshot.mass_solve(flux_difference))))
    arrays.update(ritz_sigma=sigma_ritz, physical_flux=physical_flux,
        variational_flux=variational_flux, raw_body_neumann=qbody)
    report.update(raw_field_FD=fd, raw_field_derivative_expected=float(expected),
        raw_functional_parity_absolute=float(abs(field_functional(a, sigma, u)-energy.mixed_functional(psi, sigma, u))),
        body_potential_port_quadrature_parity_max=port_error,
        body_potential_port_arithmetic_allowance=float(port_allowance),
        algebraic_reduction_controls=algebraic_controls,
        algebraic_control_scope='Old physical psi with zero and all six unit body velocities; not physics validation',
        raw_body_neumann_work=record(np.vdot(hp+R@sigma_ritz, weights*qbody)),
        old_nystrom_body_neumann_work=record(np.vdot(hp+rp, weights*(hn+rn))),
        old_physical_flux_mass_norm=float(np.sqrt(np.real(np.vdot(K@psi_old+F@sigma_old,
            snapshot.mass_solve(K@psi_old+F@sigma_old))))),
        physical_variational_flux_difference_mass_norm=float(flux_mass_norm),
        old_to_ritz_sigma_relative=float(np.linalg.norm(sigma_ritz-sigma_old)/np.linalg.norm(sigma_old)))
    checkpoint(args, report, arrays)

    started = perf_counter()
    g_error = float(abs(energy.G-energy.G.T).max())
    g_scale = float(2*abs(energy.wave_energy).max()+abs(energy.G).max())
    g_allowance = 64*np.finfo(float).eps*(nb+1)*g_scale
    report['G_transpose_gate'] = dict(maximum_absolute=g_error, arithmetic_allowance=g_allowance,
                                     scope='Floating arithmetic of the same symmetric energy Gram; no repair')
    if g_error > g_allowance:
        raise ValueError('reduced Green energy exceeds transpose roundoff gate')
    mu, vectors = eigh(energy.G, patch.mass.toarray(), check_finite=True)
    growth = np.sqrt((-GRAVITY*mu).astype(complex))
    if not np.isfinite(mu).all() or not np.isfinite(growth).all():
        raise FloatingPointError('nonfinite unmodified generalized Ritz spectrum')
    arrays.update(mu=mu, lambda_plus=growth)
    checkpoint(args, report, arrays)
    selected = np.argsort(mu)[:8]
    arrays['selected_psi'] = vectors[:, selected]
    residuals = []
    for index in selected:
        v = vectors[:, index]
        left, right = energy.G@v, mu[index]*(patch.mass@v)
        numerator, denominator = np.linalg.norm(left-right), np.linalg.norm(left)+np.linalg.norm(right)
        residuals.append(float(numerator/denominator) if denominator else (0. if numerator == 0 else None))
    report.update(complete=True, spectrum_seconds=perf_counter()-started,
        minimum_mu=float(mu.min()), maximum_mu=float(mu.max()), negative_mu_count=int((mu < 0).sum()),
        positive_reduced_energy=bool(mu.min() > 0),
        energy_gate_status=('positive quadrature-defined reduced energy' if mu.min() > 0
                            else 'rejected: reduced energy is not strictly positive; no repair'),
        maximum_growth_per_s=float(growth.real.max()), selected_relative_residuals=residuals,
        process_peak_rss=peak_rss_bytes(), limits=[
            'Symmetric Green-integrand energy defines a new Ritz discretization; no posthoc operator symmetrization',
            'Physical water mass only; covered coefficients remain algebraic extension',
            'Raw body leakage/work and physical-versus-variational flux differences are retained',
            'Zero exterior trace outside compact patch is not open-water propagation',
            'No nonlinear, moving-domain, transient, 60s motion or runtime qualification'])
    checkpoint(args, report, arrays)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'body16', 'gamma', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--maximum-pairs', type=int, default=500_000_000)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.maximum_pairs, args.maximum_bytes) <= 0:
        parser.error('fresh archive and positive work/storage bounds required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True)
    hashes = {}
    for source in [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'), *Path(__file__).parent.glob('*.py')]:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    report = dict(complete=False, producer_sha256=hashes, parents={name: sha256(getattr(args, name).read_bytes()).hexdigest()
        for name in ('water', 'phase', 'mode', 'body16', 'gamma')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
