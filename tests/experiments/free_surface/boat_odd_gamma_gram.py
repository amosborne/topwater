"""One plane-only fixed-field Gram control; cached body16 is never reevaluated."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from numpy.polynomial.legendre import leggauss
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import RegularizedSnapshot, owned
from .plane_p1_field import PlaneP1Field
from .plane_single_layer import plane_single_layer
from .panel_gradient_reference import checked_influence
from .cut_hat_geometry import cut_hat


def number(value):
    return dict(real=float(value.real), imag=float(value.imag))


def parse(value):
    return complex(value['real'], value['imag'])


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def run(args, report):
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    body = owned(args.body16)
    for name in ('water', 'phase', 'mode'):
        if body['parents'][name] != sha256(getattr(args, name).read_bytes()).hexdigest():
            raise ValueError('fixed Gram field/snapshot ownership differs')
    for name in ('plane_p1_field.py', 'plane_single_layer.py', 'cut_hat_geometry.py'):
        source = Path('tests/experiments/free_surface')/name
        if sha256(source.read_bytes()).hexdigest() != snapshot.water['producer_sha256'].get(str(source)):
            raise ValueError('recorded physical contact source changed: '+name)
    with np.load(args.body16.with_suffix('.npz')) as z:
        psi, sigma, original = z['psi'].copy(), z['sigma'].copy(), z['q16_gram'].copy()
    patch = snapshot.patch
    if sha256((patch.mass@psi).tobytes()).hexdigest() != body['physical_psi_moment_sha256']:
        raise ValueError('fixed physical psi moments changed')
    full = patch.expand(psi).reshape(len(patch.y), len(patch.x))
    field = PlaneP1Field(patch.x, patch.y, np.stack((full.real, full.imag), axis=-1))
    np.testing.assert_array_equal(field.faces, patch.faces)
    np.testing.assert_array_equal(field.panels.triangles, patch.panels.triangles)
    affected = snapshot.flux.rows
    mapping = np.full(len(patch.nodes), -1, dtype=int)
    mapping[patch.active_nodes[affected]] = np.arange(len(affected))
    parents = patch.faces[patch.parent_id]
    area_ids = np.flatnonzero(np.any(mapping[parents] >= 0, axis=1))
    hx, hy = float(patch.x[1]-patch.x[0]), float(patch.y[1]-patch.y[0])
    footprint = _polygon(args.phase)
    lines = []
    for row, active in enumerate(affected):
        geometry = cut_hat(patch.nodes[patch.active_nodes[active]], hx, hy,
                           footprint)
        node_row, node_column = divmod(patch.active_nodes[active], len(patch.x))
        field.source_parents(node_row, node_column, geometry)
        line = geometry['waterline']
        for segment, normal, values in zip(line['segments'], line['hole_outward_normals'], line['values'], strict=True):
            lines.append((row, segment, normal, values))
    ns, nb = len(field.panels), len(sigma)
    preflight = {f'K{q}': (len(area_ids)*q*q+len(lines)*q)*ns for q in (4, 8)}
    preflight.update({f'F{q}': len(patch.cut_triangles)*q*q*nb for q in (8, 16)})
    if sum(preflight.values()) > args.maximum_pairs:
        raise RuntimeError('plane-only fixed-field source-pair preflight exceeds budget')
    array_bound = 32*1024**2+32*64*ns+8*128*nb*16
    if array_bound > args.maximum_bytes:
        raise MemoryError('plane-only query/work array preflight exceeds budget')
    report.update(physical_psi_moment_sha256=body['physical_psi_moment_sha256'],
        body_quadrature_held=16, source_pair_preflight=preflight, source_pairs_so_far=0,
        array_byte_bound=array_bound,
        storage_scope='Incremental query/work arrays; excludes owned input snapshots, scalar JIT and process runtime',
        contact_rows=len(affected), physical_area_pieces=len(area_ids),
        waterline_segments=len(lines), complete=False, qualified=False, stages=[])
    arrays = dict(psi=psi, sigma=sigma, affected_rows=affected)
    pairs = 0
    counters = dict(fallback_pairs=0, exact_zero_distance_products=0)
    source_gradient = field.gradient.reshape(ns, 4)

    def reserve(targets, sources):
        nonlocal pairs
        work = targets*sources
        if pairs+work > args.maximum_pairs:
            raise RuntimeError('plane-only source-pair budget exhausted')
        pairs += work

    def vector_potential(points):
        result = np.empty((len(points), 2), dtype=complex)
        for first in range(0, len(points), 64):
            target = points[first:first+64]
            reserve(len(target), ns)
            scalar, diagnostics = plane_single_layer(target, field.panels)
            for name in counters:
                counters[name] += diagnostics[name]
            values = (scalar@source_gradient).reshape(-1, 2, 2)
            result[first:first+len(target)] = values[:, :, 0]+1j*values[:, :, 1]
        return result

    triangle = patch.panels.triangles[:, :, :2]
    gradient = np.linalg.solve(triangle[:, 1:]-triangle[:, :1],
        np.broadcast_to(np.eye(3)[1:]-np.eye(3)[:1], (len(triangle), 2, 3)))
    saved_contact = snapshot.flux.contact@psi
    plane_HH4 = np.vdot(psi, snapshot.flux(psi))
    np.testing.assert_allclose(plane_HH4, parse(body['plane_HH']), rtol=1e-11, atol=1e-13)
    plane_hh = {}
    for order in (4, 8):
        started = perf_counter()
        result = np.zeros(len(affected), dtype=complex)
        bary, weights = triangle_samples(np.eye(3)[None], order)
        bary = bary[0]
        for cut in area_ids:
            parent = patch.parent_id[cut]
            xy = bary@patch.cut_triangles[cut]
            integrated = np.einsum('q,qi->i', patch.cut_areas[cut]*weights, vector_potential(xy))
            for local, row in enumerate(mapping[patch.faces[parent]]):
                if row >= 0:
                    result[row] += 2*gradient[parent, :, local]@integrated
        nodes, weights = leggauss(order)
        t, weights = (nodes+1)/2, weights/2
        for row, segment, normal, nodal in lines:
            points = segment[:1]+t[:, None]*(segment[1:]-segment[:1])
            measure = np.linalg.norm(segment[1]-segment[0])*weights*((1-t)*nodal[0]+t*nodal[1])
            result[row] += 2*np.einsum('q,i,qi->', measure, normal, vector_potential(points))
        if order == 4:
            np.testing.assert_allclose(result, saved_contact, rtol=1e-10, atol=1e-12)
        delta = np.vdot(psi[affected], result-saved_contact)
        plane_hh[order] = plane_HH4+delta
        row = dict(component='Kcontact', order=order, seconds=perf_counter()-started,
            maximum_action_difference_from_saved_q4=float(abs(result-saved_contact).max()),
            HH_plane=number(plane_hh[order]), HH_plane_change=number(delta))
        report['stages'].append(row)
        report['source_pairs_so_far'] = pairs
        arrays[f'contact_q{order}_action'] = result
        print(json.dumps(row), flush=True)
        checkpoint(args, report, arrays)
    panels = TrianglePanels(snapshot.triangles, -1)
    plane_hr = {}
    for order in (8, 16):
        started = perf_counter()
        bary, weights = triangle_samples(np.eye(3)[None], order)
        bary = bary[0]
        work = 0j
        for first in range(0, len(patch.cut_triangles), 16):
            last = min(first+16, len(patch.cut_triangles))
            xy = np.einsum('qj,tji->tqi', bary, patch.cut_triangles[first:last])
            values = np.einsum('qj,tji->tqi', bary, patch.basis_values[first:last])
            values = np.einsum('tqi,ti->tq', values, full.ravel()[patch.faces[patch.parent_id[first:last]]]).ravel()
            measure = (patch.cut_areas[first:last, None]*weights).ravel()
            target = np.column_stack((xy.reshape(-1, 2), np.zeros(len(values))))
            for start in range(0, len(target), 128):
                stop = min(start+128, len(target))
                reserve(stop-start, nb)
                flux = 2*checked_influence(target[start:stop], panels).gradient_single_layer[:, :, 2]@sigma
                work += np.sum(measure[start:stop]*np.conj(values[start:stop])*flux)
        if order == 8:
            np.testing.assert_allclose(work, parse(body['plane_HR_q8']), rtol=1e-11, atol=1e-13)
        plane_hr[order] = work
        row = dict(component='Fwork', order=order, seconds=perf_counter()-started, HR_plane=number(work))
        report['stages'].append(row)
        report['source_pairs_so_far'] = pairs
        print(json.dumps(row), flush=True)
        checkpoint(args, report, arrays)
    fixed_body = original.copy()
    fixed_body[0, 0] -= parse(body['plane_HH'])
    fixed_body[0, 1] -= parse(body['plane_HR_q8'])
    new = fixed_body.copy()
    new[0, 0] += plane_hh[8]
    new[0, 1] += plane_hr[16]
    arrays.update(cached_body_gram=fixed_body, original_gram=original, updated_gram=new)
    report.update(complete=True, raw_HH=number(new[0, 0]), raw_RR=number(new[1, 1]),
        raw_HR=number(new[0, 1]), raw_RH=number(new[1, 0]), determinant=number(np.linalg.det(new)),
        cross_conjugacy_defect=number(new[0, 1]-np.conj(new[1, 0])),
        cross_conjugacy_relative=float(abs(new[0, 1]-np.conj(new[1, 0]))/(abs(new[0, 1])+abs(new[1, 0]))),
        scalar_source_diagnostics=counters, process_peak_rss=peak_rss_bytes(),
        limits=['Only Gamma contact q8 and source-work q16 change; cached body16 never reevaluated',
            'All saved affected rows included, all positive pieces retained; other FFT rows unchanged',
            'Collapsed q4/q8-work parity are numerical ownership controls, not independent continuum error bounds',
            'Raw complex blocks are not averaged, repaired or used as a solver',
            'No full operator/F matrix, extra body quadrature, spectrum, transient or physics qualification'])
    checkpoint(args, report, arrays)


def _polygon(phase):
    with np.load(phase.with_suffix('.npz')) as z:
        return z['polygon'].copy()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water', 'phase', 'mode', 'body16', 'output'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--maximum-pairs', type=int, default=500_000_000)
    parser.add_argument('--maximum-bytes', type=int, default=512*1024**2)
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.maximum_pairs, args.maximum_bytes) <= 0:
        parser.error('fresh archive and positive bounds required')
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
        for name in ('water', 'phase', 'mode', 'body16')})
    started = perf_counter()
    try:
        run(args, report)
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        report['total_seconds'] = perf_counter()-started
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False))
