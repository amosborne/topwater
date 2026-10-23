"""Raw2x2physical energy Gram on one owned odd mode, not a solver repair."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels, single_layer_trace
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import RegularizedSnapshot, owned
from .compact_p1_trace import CompactP1Trace
from .compact_p1_maps import body_maps
from .panel_gradient_reference import checked_influence
from .owned_odd_body_trace import owned_odd_trace


def complex_record(value):
    return dict(real=float(np.real(value)), imag=float(np.imag(value)))


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    odd, refinement = owned(args.odd), owned(args.refinement)
    parents = {name: sha256(getattr(args, name).read_bytes()).hexdigest()
               for name in ('water', 'phase', 'mode', 'odd', 'spectrum', 'body_probe')}
    if refinement['parents'] != parents or odd['parents'] != {k: parents[k] for k in ('water', 'phase', 'mode', 'body_probe')}:
        raise ValueError('Gram field/geometry ownership differs')
    owned(args.spectrum); owned(args.body_probe)
    with np.load(args.refinement.with_suffix('.npz')) as z:
        psi, sigma = z['frozen_physical_psi'].copy(), z['level0_sigma'].copy()
    with np.load(args.odd.with_suffix('.npz')) as z:
        F4, F8, jd = z['q4_F'].copy(), z['q8_F'].copy(), z['jd'].copy()
    patch, mass = snapshot.patch, snapshot.patch.mass
    if sha256((mass@psi).tobytes()).hexdigest() != refinement['physical_psi_moment_sha256']:
        raise ValueError('frozen complex-amplitude physical state changed')
    np.testing.assert_allclose(sigma, -lu_solve(lu_factor(jd), snapshot.coarse['n']@psi), rtol=1e-11, atol=1e-13)
    panels = TrianglePanels(snapshot.triangles, -1)
    source_pairs = sum(len(panels)*q*q*(2*len(patch.panels)+3*len(panels)-1) for q in args.orders)
    source_pairs += 6*len(patch.panels)+10*len(panels)-2  # two actual-point parity controls
    array_bound = 8*(2048*len(patch.panels)+64*len(panels)+32*len(panels)*max(args.orders)**2)+32*1024**2
    if source_pairs > args.maximum_pairs or array_bound > args.maximum_bytes:
        raise RuntimeError('owned Gram query/array preflight exceeds budget')
    full = patch.expand(psi).reshape(len(patch.y), len(patch.x))
    real = CompactP1Trace(patch.x, patch.y, full.real)
    imag = CompactP1Trace(patch.x, patch.y, full.imag)

    def harmonic(points, normals):
        pr, ur = real.bulk(points)
        pi, ui = imag.bulk(points)
        return pr+1j*pi, np.sum((ur+1j*ui)*normals, axis=1)

    # Independent actual-point potential/normal parity against the existing
    # all-nodal map owner before evaluating full body quadrature.
    ids = np.array([0, len(panels)//2])
    hp, hn = harmonic(panels.centroids[ids], panels.normals[ids])
    pmap, nmap = body_maps(patch, panels.centroids[ids], panels.normals[ids])
    np.testing.assert_allclose(hp, pmap@full.ravel(), atol=1e-12, rtol=3e-10)
    np.testing.assert_allclose(hn, nmap@full.ravel(), atol=1e-12, rtol=3e-10)
    rp, rn, _ = owned_odd_trace(panels, sigma, panels.centroids[ids], ids)
    trace = single_layer_trace(panels, ids)
    image = checked_influence(panels.centroids[ids]*[1., 1., -1.], panels)
    expected_r = (trace.single_layer-image.single_layer)@sigma
    expected_n = np.einsum('mni,mi->mn', trace.gradient_single_layer
        -image.gradient_single_layer*[1., 1., -1.], panels.normals[ids])@sigma
    np.testing.assert_allclose(rp, expected_r, atol=1e-12, rtol=3e-10)
    np.testing.assert_allclose(rn, expected_n, atol=1e-12, rtol=3e-10)
    plane_hh = np.vdot(psi, snapshot.flux(psi))
    plane_hr4, plane_hr8 = np.vdot(psi, F4@sigma), np.vdot(psi, F8@sigma)
    report.update(physical_psi_moment_sha256=refinement['physical_psi_moment_sha256'],
        complex_eta_amplitude_rms_m=refinement['physical_eta_rms_m'],
        normalization='Spatial physical-mass RMS of complex amplitude, not time RMS',
        gram_units='m5/s2; density/2 times Gram is complex-amplitude energy envelope, no time averaging',
        density_kg_m3=WATER_DENSITY, plane_quadrature_owners='Kphysical contactq4; Fsource q8 primary/q4 sensitivity',
        body_panels=len(panels), source_pair_preflight=source_pairs, array_byte_bound=array_bound,
        actual_harmonic_map_parity_max=float(max(abs(hp-pmap@full.ravel()).max(), abs(hn-nmap@full.ravel()).max())),
        actual_owned_trace_centroid_parity_max=float(max(abs(rp-expected_r).max(), abs(rn-expected_n).max())),
        plane_HH=complex_record(plane_hh), plane_HR_q4=complex_record(plane_hr4), plane_HR_q8=complex_record(plane_hr8),
        body_orders=[], qualified=False)
    arrays = dict(psi=psi, sigma=sigma)
    previous = None
    for order in args.orders:
        started = perf_counter()
        points, qweights = triangle_samples(snapshot.triangles, order)
        points = points.reshape(-1, 3)
        owners = np.repeat(np.arange(len(panels)), len(qweights))
        weights = (panels.areas[:, None]*qweights).ravel()
        normals = panels.normals[owners]
        hp, hn = harmonic(points, normals)
        rp, rn, diagnostics = owned_odd_trace(panels, sigma, points, owners)
        HH = plane_hh+weights@(np.conj(hp)*hn)
        RR = weights@(np.conj(rp)*rn)
        HR4 = plane_hr4+weights@(np.conj(hp)*rn)
        HR8 = plane_hr8+weights@(np.conj(hp)*rn)
        RH = weights@(np.conj(rp)*hn)
        gram = np.array([[HH, HR8], [RH, RR]])
        total_green = plane_hh+plane_hr8+weights@(np.conj(hp+rp)*(hn+rn))
        defect = HR8-np.conj(RH)
        row = dict(order=order, points=len(points), seconds=perf_counter()-started,
            HH=complex_record(HH), RR=complex_record(RR), HR_q8=complex_record(HR8), HR_q4=complex_record(HR4),
            RH=complex_record(RH), determinant=complex_record(np.linalg.det(gram)),
            cross_conjugacy_defect=complex_record(defect), cross_conjugacy_absolute=float(abs(defect)),
            cross_conjugacy_relative=float(abs(defect)/(abs(HR8)+abs(RH))) if abs(HR8)+abs(RH) else 0.,
            self_imaginary_absolute=[float(abs(HH.imag)), float(abs(RR.imag))],
            same_quadrature_total_green_sum_error=float(abs(total_green-gram.sum())),
            positive_energy_status='unresolved; raw quadrature/adjoint evidence only, no averaging',
            scalar_normal_trace_diagnostics=diagnostics)
        if previous is not None:
            row['previous_order_raw_gram_difference_max'] = float(abs(gram-previous).max())
        previous = gram
        report['body_orders'].append(row)
        arrays[f'q{order}_gram'] = gram
        arrays[f'q{order}_H'] = hp
        arrays[f'q{order}_NH'] = hn
        arrays[f'q{order}_R'] = rp
        arrays[f'q{order}_JR'] = rn
        print(json.dumps(row), flush=True)
        checkpoint(args, report, arrays)
    report.update(complete=True, process_peak_rss=peak_rss_bytes(), limits=[
        '2x2 Gram on one field pair is not a full discretization stability or passivity proof',
        'Kphysical plane quadrature is fixed; Fq4/q8 controls plane source sensitivity separately from body quadrature',
        'Both complex cross directions assembled independently; no averaging/eigenvalue repair',
        'Total Green bookkeeping uses the same quadrature and is not an independent accuracy oracle',
        'No moving-domain, nonlinear6DOF, open-water, transient or runtime qualification'])
    checkpoint(args, report, arrays)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'odd', 'spectrum', 'body-probe', 'refinement', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--orders', type=int, nargs='+', default=[2, 4, 8])
    parser.add_argument('--maximum-pairs', type=int, default=500_000_000)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    args = parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists() or args.orders != [2, 4, 8]
            or min(args.maximum_pairs, args.maximum_bytes) <= 0):
        parser.error('fresh archive, q2/4/8 and positive bounds required')
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
        for name in ('water', 'phase', 'mode', 'odd', 'spectrum', 'body_probe', 'refinement')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
