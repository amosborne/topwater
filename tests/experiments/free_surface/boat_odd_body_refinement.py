"""One same-domain odd-source68-to272 refinement at the worst owned odd mode.

Physical psi, water M/K/E, old body Gauss4 points and physical offsets remain
fixed. Only the body density/collocation partition changes. No transient.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.spatial import ConvexHull
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.rankine import TrianglePanels
from .regularized_snapshot import RegularizedSnapshot, owned
from .boat_spatial_work import source_matrix, subdivide
from .compact_p1_maps import body_maps
from .odd_body_reference import source_flux, harmonic_field, source_field
from .odd_spectrum_engine import checkpoint, run_spectrum


def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    spectrum, odd, probe = owned(args.spectrum), owned(args.odd), owned(args.body_probe)
    parents = {name: sha256(getattr(args, name).read_bytes()).hexdigest()
               for name in ('water', 'phase', 'mode', 'odd', 'body_probe')}
    if spectrum['parents'] != parents or odd['parents'] != {k: parents[k] for k in ('water', 'phase', 'mode', 'body_probe')}:
        raise ValueError('spectrum/odd physical input ownership differs')
    if probe['parents'] != {k: parents[k] for k in ('water', 'phase', 'mode')}:
        raise ValueError('old physical probe ownership differs')
    patch, mass = snapshot.patch, snapshot.patch.mass
    area = float(mass.sum())
    with np.load(args.spectrum.with_suffix('.npz')) as z:
        index = int(z['selected_indices'][0])
        psi = z['selected_physical_psi'][:, 0].copy()
        lam = complex(z['lambda_plus'][index])
    eta = -lam*psi/GRAVITY
    rms = np.sqrt(np.vdot(eta, mass@eta).real/area)
    if not np.isfinite(rms) or rms <= 0:
        raise ValueError('selected mode cannot be normalized to physical eta RMS')
    psi *= args.eta_rms/rms
    with np.load(args.odd.with_suffix('.npz')) as z:
        coarse = dict(n=snapshot.coarse['n'], j=z['jd'].copy(), F4=z['q4_F'].copy(), F8=z['q8_F'].copy())
    with np.load(args.body_probe.with_suffix('.npz')) as z:
        base, normals, weights, offsets = (z[key].copy() for key in ('old_points', 'old_normals', 'old_weights', 'old_offsets'))
    for name, value in [('old_gauss4_points', base), ('normals', normals), ('weights', weights), ('physical_offsets', offsets)]:
        if sha256(value.tobytes()).hexdigest() != probe['probe_sha256'][name]:
            raise ValueError('old physical probe changed: '+name)
    old = TrianglePanels(snapshot.triangles, -1)
    fine_triangles = subdivide(snapshot.triangles)
    fine_panels = TrianglePanels(fine_triangles, -1)
    np.testing.assert_allclose(fine_panels.areas, np.tile(old.areas/4, 4), rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(fine_panels.normals, np.tile(old.normals, (4, 1)), rtol=1e-12, atol=1e-13)
    total_queries = 2*len(base)*len(patch.panels)
    total_queries += 2*len(base)*2*(len(old)+len(fine_panels))
    probe_bytes = 8*(2*len(base)*len(patch.nodes)+1024*len(patch.panels)+4096*len(fine_panels))
    if total_queries > args.maximum_pairs or probe_bytes > args.maximum_bytes:
        raise RuntimeError('fixed-probe work/storage bound exceeds budget')
    report.update(physical_psi_moment_sha256=sha256((mass@psi).tobytes()).hexdigest(),
        physical_eta_rms_m=args.eta_rms, frozen_mode_index=index,
        frozen_mode_lambda_real=lam.real, frozen_mode_lambda_imag=lam.imag,
        physical_inputs_sha256=snapshot.mode['input_sha256'], panel_sha256=snapshot.mode['panel_sha256'],
        body_levels=[], harmonic_probe_evaluations=2,
        probe_sha256=probe['probe_sha256'],
        probe_source_pair_bound=total_queries,
        probe_array_bytes_bound=probe_bytes,
        refinement_area_difference=float(fine_panels.areas.sum()-old.areas.sum()),
        refinement_first_moment_difference=((fine_panels.areas[:, None]*fine_panels.centroids).sum(axis=0)
            -(old.areas[:, None]*old.centroids).sum(axis=0)).tolist(),
        candidate_limits=['Only one body-density refinement; failure stops blind body refinement',
            'Zero physical potential outside the patch is truncated Dirichlet data, not open-ocean propagation',
            'Body Gauss4/h-halving and Fq4/q8 are sampled sensitivity, not certified continuum convergence'])
    saved = dict(frozen_physical_psi=psi, old_points=base, old_normals=normals, old_weights=weights, old_offsets=offsets)
    # Hpsi is unchanged across both body spaces. Compute each physical offset
    # exactly ONCE; only odd body-source contributions are subsequently rebuilt.
    solid = ConvexHull(snapshot.triangles.reshape(-1, 3))
    harmonic = []
    targets = []
    for fraction in (1., .5):
        target = base-fraction*offsets[:, None]*normals
        if (np.any(target[:, 2] >= 0) or np.any(np.all(target == base, axis=1))
                or np.any(np.max(target@solid.equations[:, :3].T+solid.equations[:, 3], axis=1) <= 0)):
            raise ValueError('fixed old normal probes fail exterior certificate')
        targets.append(target)
        harmonic.append(harmonic_field(patch, psi, target, normals))
    fine = None
    for level, panels in enumerate((old, fine_panels)):
        started = perf_counter()
        if level == 0:
            system = coarse
            construction = dict(scope='SHA-owned archived68 N/JD/F maps')
        else:
            map_pairs = len(panels)*len(patch.panels)+2*len(panels)**2
            map_bytes = 8*(4*len(panels)*len(patch.nodes)+1024*len(patch.panels))
            if map_pairs > args.maximum_pairs or map_bytes > args.maximum_bytes:
                raise RuntimeError('refined N/JD construction exceeds bounded work/storage')
            full_p, full_n = body_maps(patch, panels.centroids, panels.normals)
            nmap = patch.condense_maps(full_n)
            del full_p, full_n
            system = dict(n=nmap, j=source_matrix(panels))
            construction = dict(map_source_pair_bound=map_pairs, map_array_byte_bound=map_bytes, orders=[])
            for order in (4, 8):
                system[f'F{order}'], diagnostics = source_flux(patch, panels, order,
                    maximum_pairs=args.maximum_pairs, maximum_bytes=args.maximum_bytes)
                construction['orders'].append(dict(order=order, **diagnostics))
        factor = lu_factor(system['j'])
        sigma = -lu_solve(factor, system['n']@psi)
        row = dict(level=level, body_panels=len(panels), construction=construction,
            centroid_normal_residual_max=float(abs(system['n']@psi+system['j']@sigma).max()), orders=[], probes=[])
        report['body_levels'].append(row)
        saved[f'level{level}_sigma'] = sigma
        for order in (4, 8):
            weak = snapshot.flux(psi)+system[f'F{order}']@sigma
            rate = snapshot.mass_solve(weak)
            quotient = np.vdot(psi, weak)/np.vdot(psi, mass@psi)
            row['orders'].append(dict(order=order, rayleigh_real=float(quotient.real), rayleigh_imag=float(quotient.imag),
                eta_rate_mass_norm=float(np.sqrt(np.vdot(rate, mass@rate).real/area))))
            if level:
                difference = rate-saved[f'level0_q{order}_eta_rate']
                row['orders'][-1]['coarse_fine_action_difference_mass_norm'] = float(np.sqrt(np.vdot(difference, mass@difference).real/area))
            saved[f'level{level}_q{order}_eta_rate'] = rate
        qdifference = saved[f'level{level}_q8_eta_rate']-saved[f'level{level}_q4_eta_rate']
        row['q4_q8_action_difference_mass_norm'] = float(np.sqrt(np.vdot(qdifference, mass@qdifference).real/area))
        for offset_id, fraction in enumerate((1., .5)):
            sp, sn = source_field(panels, sigma, targets[offset_id], normals)
            hp, hn = harmonic[offset_id]
            value, normal = hp+sp, hn+sn
            record = dict(offset_factor=fraction,
                potential_rms=float(np.sqrt(weights@abs(value)**2/weights.sum())),
                normal_rms=float(np.sqrt(weights@abs(normal)**2/weights.sum())), normal_max=float(abs(normal).max()),
                green_body_work_J=float(WATER_DENSITY*np.real(weights@(np.conj(value)*normal))))
            row['probes'].append(record)
            saved[f'level{level}_h{fraction}_potential'] = value
            saved[f'level{level}_h{fraction}_normal'] = normal
        row['seconds'] = perf_counter()-started
        print(json.dumps(row), flush=True)
        checkpoint(args, report, saved)
        if level == 1:
            fine = system
            for name, value in fine.items():
                saved['fine_'+name] = value
    run_spectrum(snapshot, fine['n'], fine['j'], fine['F4'], fine['F8'], psi, args, report,
        reference_rate=saved['level1_q8_eta_rate'], arrays=saved)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'odd', 'spectrum', 'body-probe', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--eta-rms', type=float, default=.001)
    parser.add_argument('--batch-columns', type=int, default=64)
    parser.add_argument('--maximum-pairs', type=int, default=500_000_000)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or not 0 < args.eta_rms <= .01 or min(args.batch_columns, args.maximum_pairs, args.maximum_bytes) <= 0):
        parser.error('fresh archive and positive bounded controls required')
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
        for name in ('water', 'phase', 'mode', 'odd', 'spectrum', 'body_probe')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
