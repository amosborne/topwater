"""Odd-image physical-P1 trace: frozen action/leak and bounded spectrum gate.

Zero trace outside the patch is homogeneous truncated Dirichlet data, NOT an
open-ocean propagation model. No transient, filter, symmetrization or6DOF claim.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.sparse.linalg import LinearOperator, eigs, ArpackNoConvergence
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.rankine import TrianglePanels
from .boat_spatial_work import source_matrix
from .regularized_snapshot import RegularizedSnapshot, owned
from .odd_body_reference import source_flux, harmonic_field, source_field


def checkpoint(args, report, saved):
    np.savez(args.output.with_suffix('.npz'), **saved)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    reference = owned(args.body_probe)
    expected_parents = {name: sha256(getattr(args, name).read_bytes()).hexdigest() for name in ('water', 'phase', 'mode')}
    if reference['parents'] != expected_parents or reference['physical_inputs_sha256'] != snapshot.mode['input_sha256']:
        raise ValueError('old body-probe/input ownership differs')
    with np.load(args.body_probe.with_suffix('.npz')) as z:
        psi = z['physical_psi'].copy()
        base, normals, weights, offsets = (z[key].copy() for key in
            ('old_points', 'old_normals', 'old_weights', 'old_offsets'))
    for name, value in [('old_gauss4_points', base), ('normals', normals), ('weights', weights), ('physical_offsets', offsets)]:
        if sha256(value.tobytes()).hexdigest() != reference['probe_sha256'][name]:
            raise ValueError('old physical probe ownership differs: '+name)
    patch, mass = snapshot.patch, snapshot.patch.mass
    if sha256((mass@psi).tobytes()).hexdigest() != reference['physical_psi_moment_sha256']:
        raise ValueError('frozen physical psi moments changed')
    panels = TrianglePanels(snapshot.triangles, -1)
    count, nb = len(psi), len(panels)
    jd = source_matrix(panels)
    factor = lu_factor(jd)
    nmap = snapshot.coarse['n']
    sigma = -lu_solve(factor, nmap@psi)
    # Same old psi coefficient field is now the ACTUAL plane trace; unlike
    # the even-image model there is no projected Cw source trace or Schur.
    report.update(physical_psi_moment_sha256=reference['physical_psi_moment_sha256'],
        physical_inputs_sha256=snapshot.mode['input_sha256'], panel_sha256=snapshot.mode['panel_sha256'],
        active_nodes=count, body_panels=nb, body_collocation_residual_max=float(abs(nmap@psi+jd@sigma).max()),
        cpu_saved_operator_audit_error=snapshot.flux.audit_error, probes=[], orders=[])
    saved = dict(physical_psi=psi, sigma=sigma, jd=jd)
    preflight = [len(patch.cut_triangles)*q*q*nb for q in args.orders]
    if any(work > args.maximum_pairs for work in preflight) or 24*len(patch.nodes)*nb > args.maximum_bytes:
        raise RuntimeError('odd source-flux assembly preflight exceeds work/storage budget')
    probe_pairs = 2*len(base)*(len(patch.panels)+3*nb)
    probe_bytes = 8*(2*len(base)*len(patch.nodes)+1024*len(patch.panels)+4096*nb)
    if probe_pairs > args.maximum_pairs or probe_bytes > args.maximum_bytes:
        raise RuntimeError('odd body probe exceeds work/storage preflight')
    report.update(F_source_pair_preflight=preflight, probe_source_pair_bound=probe_pairs,
        probe_array_byte_bound=probe_bytes,
        eigen_coordinate_scaling='sqrt physical mass; psi also divided sqrt(g*beam); similarity only')
    # SAME physical old Gauss4 points and h/h2 offsets as the even-source test.
    for fraction in (1., .5):
        target = base-fraction*offsets[:, None]*normals
        hp, hn = harmonic_field(patch, psi, target, normals)
        sp, sn = source_field(panels, sigma, target, normals)
        value, normal = hp+sp, hn+sn
        record = dict(offset_factor=fraction, normal_rms=float(np.sqrt(weights@abs(normal)**2/weights.sum())),
            normal_max=float(abs(normal).max()), potential_rms=float(np.sqrt(weights@abs(value)**2/weights.sum())),
            green_body_work_J=float(WATER_DENSITY*np.real(weights@(np.conj(value)*normal))))
        report['probes'].append(record)
        saved[f'h{fraction}_potential'] = value
        saved[f'h{fraction}_normal'] = normal
        print(json.dumps(record), flush=True)
    previous_load = previous_F = None
    for order in args.orders:
        F, flux_diagnostics = source_flux(patch, panels, order,
            maximum_pairs=args.maximum_pairs, maximum_bytes=args.maximum_bytes)
        load = snapshot.flux(psi)+F@sigma
        rate = snapshot.mass_solve(load)
        rayleigh = np.vdot(psi, load)/np.vdot(psi, mass@psi)
        row = dict(order=order, **flux_diagnostics,
            effective_rayleigh_real=float(rayleigh.real), effective_rayleigh_imag=float(rayleigh.imag),
            frozen_state_rate_mass_norm=float(np.sqrt(np.vdot(rate, mass@rate).real/float(mass.sum()))),
            spectrum=dict(status='not started'))
        if previous_load is not None:
            difference = snapshot.mass_solve(load-previous_load)
            row['source_flux_qchange_rate_mass_norm'] = float(np.sqrt(np.vdot(difference, mass@difference).real/float(mass.sum())))
            row['F_qchange_relative_frobenius'] = float(np.linalg.norm(F-previous_F)/np.linalg.norm(F))
        previous_load, previous_F = load, F
        report['orders'].append(row)
        saved[f'q{order}_F'] = F
        saved[f'q{order}_eta_rate'] = rate
        checkpoint(args, report, saved)
        print(json.dumps(dict(order=order, frozen_action=row)), flush=True)
        scaling = np.sqrt(mass.diagonal())
        velocity = np.sqrt(GRAVITY*np.ptp(snapshot.inputs.body.mesh.vertices[:, 1]))
        evaluations = 0
        def action(v):
            nonlocal evaluations
            evaluations += 1
            if evaluations > args.maximum_matvecs:
                raise RuntimeError('odd eigen matvec budget exhausted')
            physical = velocity*v[count:]/scaling
            source = -lu_solve(factor, nmap@physical)
            load = snapshot.flux(physical)+F@source
            return np.r_[scaling*snapshot.mass_solve(load), -GRAVITY/velocity*v[:count]]
        op = LinearOperator((2*count, 2*count), matvec=action, dtype=float)
        try:
            values, vectors = eigs(op, k=4, which='LR', ncv=40, tol=1e-8, maxiter=250,
                v0=np.sin(np.arange(2*count)*.371)+np.cos(np.arange(2*count)*.613))
            status = 'converged'
        except ArpackNoConvergence as error:
            values, vectors, status = error.eigenvalues, error.eigenvectors, 'partial convergence'
        except RuntimeError as error:
            values, vectors, status = np.empty(0, complex), np.empty((2*count, 0), complex), str(error)
        modes = []
        for value, vector in zip(values, vectors.T, strict=True):
            eta, physical = vector[:count]/scaling, velocity*vector[count:]/scaling
            source = -lu_solve(factor, nmap@physical)
            weak = snapshot.flux(physical)+F@source
            modes.append(dict(real=float(value.real), imag=float(value.imag),
                weak_kinematic_relative_residual=float(np.linalg.norm(weak-value*(mass@eta))/
                    (np.linalg.norm(weak)+abs(value)*np.linalg.norm(mass@eta))),
                weak_dynamic_relative_residual=float(np.linalg.norm(-GRAVITY*(mass@eta)-value*(mass@physical))/
                    (GRAVITY*np.linalg.norm(mass@eta)+abs(value)*np.linalg.norm(mass@physical))),
                body_collocation_residual_max=float(abs(nmap@physical+jd@source).max())))
        row['spectrum'] = dict(status=status, matvecs=evaluations, modes=modes)
        saved[f'q{order}_eigenvalues'] = values
        saved[f'q{order}_eigenvectors'] = vectors
        checkpoint(args, report, saved)
    report.update(complete=True, qualified=False, limits=[
        'Changes both source representation and physical trace/far closure; not isolated causal attribution',
        'Physical trace is P1 with zero exterior-patch Dirichlet data, not open-ocean propagation',
        'Gauss4 body work and source-flux q4/q8 changes are sampled sensitivity, not error certificates',
        'Bounded spectral search can be incomplete; no claim of global passivity or nonlinear6DOF qualification'])
    checkpoint(args, report, saved)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'body-probe', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--orders', type=int, nargs='+', default=[4, 8])
    parser.add_argument('--maximum-pairs', type=int, default=500_000_000)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    parser.add_argument('--maximum-matvecs', type=int, default=800)
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists() or args.orders != [4, 8]
            or min(args.maximum_pairs, args.maximum_bytes, args.maximum_matvecs) <= 0):
        parser.error('fresh output/archive and declared q4/q8 orders required')
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
        for name in ('water', 'phase', 'mode', 'body_probe')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
