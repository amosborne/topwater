"""Same-domain body-density refinement, frozen physical psi and weak flux.

No Capytaine, transient integration, filter, eigenvalue repair or changed hull.
Off-centroid work remains a sampled one-sided diagnostic, not a boundary norm.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import lu_solve
from scipy.sparse.linalg import LinearOperator, eigs, ArpackNoConvergence
from scipy.spatial import ConvexHull
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from .boat_spatial_work import subdivide
from .regularized_snapshot import RegularizedSnapshot


def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    if args.order != snapshot.phase['orders'][-1]['order']:
        raise ValueError('body refinement must retain the archived physical Cw quadrature order')
    mass = snapshot.patch.mass
    area = float(mass.sum())
    scale = args.eta_rms/np.sqrt(np.vdot(snapshot.eta, mass@snapshot.eta).real/area)
    psi = snapshot.psi*scale
    report.update(cpu_saved_operator_audit_error=snapshot.flux.audit_error,
        physical_psi_moment_sha256=sha256((mass@psi).tobytes()).hexdigest(),
        physical_inputs_sha256=snapshot.mode['input_sha256'], old_panel_sha256=snapshot.mode['panel_sha256'],
        quadrature_order=args.order, fixed_far_port='zero homogeneous; no prescribed far forcing',
        body_levels=[], complete=False, qualified=False)
    # These are OLD physical points/weights and OLD physical h/h2 offsets for
    # BOTH body-density spaces. Refinement does not relocate the samples.
    old = TrianglePanels(snapshot.triangles, -1)
    points, weights = triangle_samples(snapshot.triangles, 4)
    parents = np.repeat(np.arange(len(old)), len(weights))
    base = points.reshape(-1, 3)
    normals = old.normals[parents]
    w = (old.areas[:, None]*weights).ravel()
    h = args.offset_fraction*old.diameters[parents]
    rise = normals[:, 2] < 0
    h[rise] = np.minimum(h[rise], .25*abs(base[rise, 2])/(-normals[rise, 2]))
    solid = ConvexHull(snapshot.triangles.reshape(-1, 3))
    probe_pairs = 2*sum(len(base)*(len(snapshot.patch.panels)+2*nb) for nb in (len(old), 4*len(old)))
    probe_bytes = 8*(2*len(base)*len(snapshot.patch.nodes)+1024*len(snapshot.patch.panels)+4096*4*len(old))
    if probe_pairs > args.maximum_pairs or probe_bytes > args.maximum_bytes:
        raise RuntimeError('same old-point body probe exceeds bounded work/storage preflight')
    report['probe_pair_work_upper_bound'] = probe_pairs
    report['probe_array_bytes_upper_bound'] = probe_bytes
    report['probe_sha256'] = {name: sha256(a.tobytes()).hexdigest() for name, a in
        [('old_gauss4_points', base), ('normals', normals), ('weights', w), ('physical_offsets', h)]}
    saved = dict(physical_psi=psi, old_points=base, old_normals=normals, old_weights=w, old_offsets=h)
    previous_rate = None
    refined_system = None
    for level, triangles in enumerate((snapshot.triangles, subdivide(snapshot.triangles))):
        started = perf_counter()
        system = snapshot.coarse if level == 0 else snapshot.body_system(triangles, args.order,
            maximum_pairs=args.maximum_pairs, maximum_bytes=args.maximum_bytes)
        chi, sigma, yc, factor = snapshot.solve(system, psi)
        load = snapshot.flux(chi)
        rate = snapshot.mass_solve(load)
        body_residual = system['n']@chi+system['j']@sigma
        trace_residual = mass@chi+system['cross']@sigma-mass@psi
        row = dict(level=level, body_panels=len(triangles),
            physical_cross_quadrature_order=args.order,
            triangles_sha256=sha256(triangles.tobytes()).hexdigest(),
            body_normal_collocation_residual_max=float(abs(body_residual).max()),
            physical_trace_moment_residual_max=float(abs(trace_residual).max()),
            physical_trace_moment_residual_relative=float(np.linalg.norm(trace_residual)/np.linalg.norm(mass@psi)),
            frozen_state_rate_mass_norm=float(np.sqrt(np.vdot(rate, mass@rate).real/area)),
            effective_rayleigh_real=float((np.vdot(psi, load)/np.vdot(psi, mass@psi)).real),
            effective_rayleigh_imag=float((np.vdot(psi, load)/np.vdot(psi, mass@psi)).imag),
            probes=[])
        if previous_rate is not None:
            difference = rate-previous_rate
            row['same_physical_state_rate_difference_mass_norm'] = float(np.sqrt(np.vdot(difference, mass@difference).real/area))
        previous_rate = rate
        report['body_levels'].append(row)
        saved[f'level{level}_chi'] = chi
        saved[f'level{level}_sigma'] = sigma
        saved[f'level{level}_eta_rate'] = rate
        print(json.dumps(dict(level=level, panels=len(triangles), frozen_state_solved=True)), flush=True)
        for fraction in (1., .5):
            target = base-fraction*h[:, None]*normals
            if (np.any(target[:, 2] >= 0) or np.any(np.all(target == base, axis=1))
                    or np.any(np.max(target@solid.equations[:, :3].T+solid.equations[:, 3], axis=1) <= 0)):
                raise ValueError('fixed old normal probes fail physical exterior certificate')
            p, q = snapshot.body_field(triangles, chi, sigma, target, normals)
            record = dict(offset_factor=fraction,
                potential_rms=float(np.sqrt(w@abs(p)**2/w.sum())),
                normal_rms=float(np.sqrt(w@abs(q)**2/w.sum())), normal_max=float(abs(q).max()),
                green_body_work_J=float(WATER_DENSITY*np.real(w@(np.conj(p)*q))))
            row['probes'].append(record)
            saved[f'level{level}_h{fraction}_potential'] = p
            saved[f'level{level}_h{fraction}_normal'] = q
            print(json.dumps(dict(level=level, **record)), flush=True)
        row['seconds'] = perf_counter()-started
        if level == 1:
            refined_system = (system, yc, factor)
        checkpoint(args, report, saved)
    # Optional bounded refined spectrum AFTER same-state leakage/action data.
    if args.maximum_matvecs:
        system, yc, factor = refined_system
        count = len(psi)
        scaling = np.sqrt(mass.diagonal())
        velocity = np.sqrt(GRAVITY*np.ptp(snapshot.inputs.body.mesh.vertices[:, 1]))
        evaluations = 0
        def action(v):
            nonlocal evaluations
            evaluations += 1
            if evaluations > args.maximum_matvecs:
                raise RuntimeError('refined eigen matvec budget exhausted')
            eta, physical = v[:count]/scaling, velocity*v[count:]/scaling
            chi = physical+yc@lu_solve(factor, system['n']@physical)
            return np.r_[scaling*snapshot.mass_solve(snapshot.flux(chi)), -GRAVITY/velocity*v[:count]]
        operator = LinearOperator((2*count, 2*count), matvec=action, dtype=float)
        try:
            values, vectors = eigs(operator, k=4, which='LR', ncv=40, tol=1e-8, maxiter=250,
                v0=np.sin(np.arange(2*count)*.371)+np.cos(np.arange(2*count)*.613))
            status = 'converged'
        except ArpackNoConvergence as error:
            values, vectors, status = error.eigenvalues, error.eigenvectors, 'partial convergence'
        except RuntimeError as error:
            values, vectors, status = np.empty(0, complex), np.empty((2*count, 0), complex), str(error)
        modes = []
        for value, vector in zip(values, vectors.T, strict=True):
            eta, physical = vector[:count]/scaling, velocity*vector[count:]/scaling
            sigma = -lu_solve(factor, system['n']@physical)
            chi = physical-yc@sigma
            load = snapshot.flux(chi)
            weak_kinematic = load-value*(mass@eta)
            weak_dynamic = -GRAVITY*(mass@eta)-value*(mass@physical)
            modes.append(dict(real=float(value.real), imag=float(value.imag),
                weak_kinematic_relative_residual=float(np.linalg.norm(weak_kinematic)/
                    (np.linalg.norm(load)+abs(value)*np.linalg.norm(mass@eta))),
                weak_dynamic_relative_residual=float(np.linalg.norm(weak_dynamic)/
                    (GRAVITY*np.linalg.norm(mass@eta)+abs(value)*np.linalg.norm(mass@physical))),
                body_normal_residual_max=float(abs(system['n']@chi+system['j']@sigma).max()),
                physical_trace_moment_residual_max=float(abs(mass@chi+system['cross']@sigma-mass@physical).max())))
        report['refined_spectrum'] = dict(status=status, matvecs=evaluations, values=modes)
        saved['refined_eigenvalues'] = values
        saved['refined_eigenvectors'] = vectors
    report.update(complete=True, limits=['Only body-density/collocation partition changes; old physical psi moments and flux operator fixed',
        'Gauss4 and h/h2 body work are sampled diagnostics; quadrature convergence is not established',
        'No transient, nonlinear motion, passivity repair or physical accuracy qualification'])
    checkpoint(args, report, saved)


def checkpoint(args, report, saved):
    np.savez(args.output.with_suffix('.npz'), **saved)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--order', type=int, choices=[4, 8, 16], default=8)
    parser.add_argument('--eta-rms', type=float, default=.001)
    parser.add_argument('--offset-fraction', type=float, default=1e-4)
    parser.add_argument('--maximum-pairs', type=int, default=500_000_000)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    parser.add_argument('--maximum-matvecs', type=int, default=800)
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or not 0 < args.eta_rms <= .01 or not 0 < args.offset_fraction < .01):
        parser.error('fresh output and positive bounded state/offset scales required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True)
    hashes = {}
    for source in [*Path('src/topwater').glob('*.py'), *Path('tests').glob('*.py'), *Path(__file__).parent.glob('*.py')]:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    report = dict(complete=False, producer_sha256=hashes,
        parents={name: sha256(getattr(args, name).read_bytes()).hexdigest() for name in ('water', 'phase', 'mode')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
