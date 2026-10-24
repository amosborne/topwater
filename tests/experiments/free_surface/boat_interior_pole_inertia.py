"""One fixed 68-interior-pole Ritz infinity family on the owned physical hull."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import cholesky, cho_solve, eigvalsh, solve_triangular
from scipy.spatial import ConvexHull
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import RegularizedSnapshot, owned


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def ratio(a, b):
    return float(a/b) if b else (0. if a == 0 else None)


def certified_poles(triangles):
    """Exact signs of ORIGINAL binary64 data; Qhull only discovers facets."""
    vertices = np.unique(triangles.reshape(-1, 3), axis=0)
    hull = ConvexHull(vertices)
    center = vertices.mean(axis=0)
    poles = .75*triangles.mean(axis=1)+.25*center
    rational = [tuple(Fraction(float(v)) for v in point) for point in vertices]
    targets = [tuple(Fraction(float(v)) for v in point) for point in np.vstack((center, poles))]
    clearance = []
    for facet in hull.simplices:
        a, b, c = [rational[i] for i in facet]
        ab, ac = [b[i]-a[i] for i in range(3)], [c[i]-a[i] for i in range(3)]
        normal = (ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0])
        if not any(normal):
            raise ValueError('degenerate discovered closed-solid facet')
        offset = sum(normal[i]*a[i] for i in range(3))
        signs = [sum(normal[i]*v[i] for i in range(3))-offset for v in rational]
        if any(s > 0 for s in signs) and any(s < 0 for s in signs):
            raise ValueError('Qhull facet is not an exact original-float supporting plane')
        if any(s > 0 for s in signs):
            normal, offset = tuple(-v for v in normal), -offset
        if not any(signs):
            raise ValueError('closed wet solid has no three-dimensional interior')
        distances = [offset-sum(normal[i]*v[i] for i in range(3)) for v in targets]
        if any(d <= 0 for d in distances):
            raise ValueError('fixed inset center/pole is not strictly inside exact wet-solid halfspaces')
        scale = np.linalg.norm([float(v) for v in normal])
        clearance.append([float(d)/scale for d in distances])
    if any(Fraction(float(z)) >= 0 for z in poles[:, 2]):
        raise ValueError('all fixed poles must lie strictly below the plane')
    image = poles*np.array([1., 1., -1.])
    if np.any(image[:, 2] <= 0):
        raise ValueError('reflected poles must lie strictly above the plane')
    return poles, dict(alpha=.25, center=center.tolist(), unique_wet_vertices=len(vertices),
        closed_support_facets=len(hull.simplices), exact_signs='Fraction.from_float original coordinates; strict, no tolerance',
        center_minimum_clearance_m=float(np.min(np.asarray(clearance)[:, 0])),
        pole_minimum_clearance_m=float(np.min(np.asarray(clearance)[:, 1:])),
        pole_maximum_z_m=float(poles[:, 2].max()), image_minimum_z_m=float(image[:, 2].min()))


def point_maps(points, normals, poles):
    """Ordinary target derivatives of positive 1/r minus its flat image."""
    potential = np.zeros((len(points), len(poles)))
    normal = np.zeros_like(potential)
    for sign, source in ((1., poles), (-1., poles*[1., 1., -1.])):
        delta = points[:, None]-source
        radius = np.linalg.norm(delta, axis=2)
        if not np.isfinite(radius).all() or np.any(radius <= 0):
            raise ValueError('finite ordinary point-kernel targets distinct from poles required')
        inverse = 1/(4*np.pi*radius)
        potential += sign*inverse
        normal -= sign*np.einsum('mni,mi->mn', delta, normals)*(inverse/radius**2)
    if not np.isfinite(potential).all() or not np.isfinite(normal).all():
        raise FloatingPointError('nonfinite ordinary pole trace')
    return potential, normal


def run(args, report):
    reference = owned(args.reference)
    snapshot = RegularizedSnapshot(args.water, args.phase, args.mode)
    for name in ('water', 'phase', 'mode'):
        if reference['parents'][name] != sha256(getattr(args, name).read_bytes()).hexdigest():
            raise ValueError('cached six-port references belong to different physical parents')
    with np.load(args.reference.with_suffix('.npz')) as z:
        cached = {name: z[name].copy() for name in ('ritz_added_mass', 'direct_added_mass',
            'capytaine_added_mass', 'dry_mass', 'triangles', 'world_cg')}
    np.testing.assert_array_equal(cached['triangles'], snapshot.triangles)
    np.testing.assert_array_equal(cached['world_cg'], snapshot.inputs.initial_state.position)
    panels = TrianglePanels(snapshot.triangles, -1)
    poles, certificate = certified_poles(snapshot.triangles)
    b = len(poles)
    if b != 68:
        raise ValueError('this verification family requires the owned 68 physical-face centroids')
    pair_bound = 2*b*b*(16**2+32**2+64**2+4**2)
    if pair_bound > args.maximum_pairs:
        raise RuntimeError('predeclared pole-family source-pair budget exceeded')
    # Four all-order target coordinate/weight tensors plus streamed kernels,
    # no Q-by-body field map retained for q16/32/64.
    array_bound = 8*(b*64**2*16+1024*b*8+8*b*b)+16*1024**2
    if array_bound > args.maximum_bytes:
        raise MemoryError('pole-family array preflight exceeds budget')
    factor_dry = cholesky(cached['dry_mass'], lower=True)

    def whiten(value):
        left = solve_triangular(factor_dry, value, lower=True)
        return solve_triangular(factor_dry, left.T, lower=True).T

    center = cached['world_cg']
    qpoints, qw = triangle_samples(panels.triangles, 4)
    qpoints = qpoints.reshape(-1, 3)
    qowners = np.repeat(np.arange(b), len(qw))
    qnormal = panels.normals[qowners]
    qT = np.column_stack((qnormal, np.cross(qpoints-center, qnormal)))
    qweights = (panels.areas[:, None]*qw).ravel()
    R4, J4 = point_maps(qpoints, qnormal, poles)
    pairs = 2*len(qpoints)*b
    report.update(complete=False, qualified=False, placement=certificate, source_pair_bound=pair_bound,
        array_workspace_bound_bytes=array_bound, body_orders=[], density_kg_m3=WATER_DENSITY,
        source_pairs_so_far=pairs,
        fixed_control=dict(order=4, points=len(qpoints), point_sha256=sha256(qpoints.tobytes()).hexdigest(),
            scope='Same OLD physical Gauss4 points for every pole-family energy quadrature'),
        archived_refinement=reference['archived_body_refinement'])
    arrays = dict(poles=poles, triangles=snapshot.triangles, world_cg=center,
        fixed_q4_points=qpoints, fixed_q4_weights=qweights, fixed_q4_T=qT,
        dry_mass=cached['dry_mass'], **{k:v for k,v in cached.items() if k.endswith('added_mass')})
    checkpoint(args, report, arrays)
    previous = None
    for order in (16, 32, 64):
        started = perf_counter()
        points, weights = triangle_samples(panels.triangles, order)
        points = points.reshape(-1, 3)
        owners = np.repeat(np.arange(b), len(weights))
        weights = (panels.areas[:, None]*weights).ravel()
        RR, BR = np.zeros((b, b)), np.zeros((6, b))
        for first in range(0, len(points), 1024):
            last = min(first+1024, len(points))
            pairs += 2*(last-first)*b
            if pairs > args.maximum_pairs:
                raise RuntimeError('pole-family source-pair counter exceeded predeclared budget')
            normal = panels.normals[owners[first:last]]
            R, J = point_maps(points[first:last], normal, poles)
            T = np.column_stack((normal, np.cross(points[first:last]-center, normal)))
            w = weights[first:last, None]
            RR += R.T@(w*J)
            BR += T.T@(w*R)
        E = .5*(RR+RR.T)  # Definition of this scalar symmetric-Green energy.
        row = dict(order=order, points=len(points), raw_adjoint_relative=ratio(np.linalg.norm(RR-RR.T), 2*np.linalg.norm(RR)),
            assembly_seconds=perf_counter()-started, positive_energy=False)
        report['source_pairs_so_far'] = pairs
        arrays.update({f'q{order}_raw_RR':RR, f'q{order}_BR':BR})
        report['body_orders'].append(row)
        checkpoint(args, report, arrays)
        try:
            diagonal = E.diagonal()
            if np.any(diagonal <= 0):
                raise ValueError('point-source energy has a nonpositive diagonal')
            scale = np.sqrt(diagonal)
            scaled = E/scale[:, None]/scale[None, :]
            eigenvalues = eigvalsh(scaled)
            rank_allowance = b*np.finfo(float).eps*max(abs(eigenvalues))
            row.update(scaled_energy_eigenvalues=eigenvalues.tolist(),
                arithmetic_rank_threshold=float(rank_allowance),
                resolved_rank=int((eigenvalues > rank_allowance).sum()))
            if eigenvalues[0] <= 0 or row['resolved_rank'] != b:
                raise ValueError('point-source energy is negative or numerically rank deficient; no repair')
            source_factor = cholesky(scaled, lower=True)
            sigma = cho_solve((source_factor, True), BR.T/scale[:, None])/scale[:, None]
            if not np.isfinite(sigma).all():
                raise ValueError('nonfinite source strengths; no regularization')
        except (ValueError, np.linalg.LinAlgError) as error:
            row['rejected'] = dict(type=type(error).__name__, message=str(error))
            row['seconds'] = perf_counter()-started
            checkpoint(args, report, arrays)
            print(json.dumps(row), flush=True)
            continue
        A = WATER_DENSITY*BR@sigma
        potential, normal = R4@sigma, J4@sigma
        leak = np.sqrt((qweights[:, None]*(normal-qT)**2).sum(axis=0)/qweights.sum())
        rigid = np.sqrt((qweights[:, None]*qT**2).sum(axis=0)/qweights.sum())
        control_mass = WATER_DENSITY*qT.T@(qweights[:, None]*potential)
        row.update(positive_energy=True, diagonal=A.diagonal().tolist(),
            fixed_q4_normal_relative_rms=[ratio(a, b) for a,b in zip(leak, rigid, strict=True)],
            fixed_q4_normal_max=abs(normal-qT).max(axis=0).tolist(),
            fixed_q4_integrated_pressure_added_mass=control_mass.tolist(),
            fixed_q4_pressure_port_difference_max=float(abs(control_mass-A).max()),
            algebraic_stationarity_max=float(abs(E@sigma-BR.T).max()),
            source_strength_l2=np.linalg.norm(sigma, axis=0).tolist(),
            source_strength_units='point strengths; not panel densities, no panel-area source normalization',
            scaled_energy_condition=float(eigenvalues[-1]/eigenvalues[0]))
        for name in ('ritz_added_mass', 'direct_added_mass', 'capytaine_added_mass'):
            comparison = A-cached[name]
            white = whiten(comparison)
            row[name+'_comparison'] = dict(diagonal_relative_difference=[ratio(a,b) for a,b in
                zip(comparison.diagonal(),cached[name].diagonal(),strict=True)],
                dry_mass_whitened_difference=white.tolist(),
                dry_mass_whitened_relative_operator_norm=ratio(np.linalg.norm(white,2),np.linalg.norm(whiten(cached[name]),2)))
        if previous is not None:
            row['previous_order_added_mass_difference_max'] = float(abs(A-previous).max())
        previous = A
        row['seconds'] = perf_counter()-started
        arrays.update({f'q{order}_added_mass':A, f'q{order}_sigma':sigma,
            f'q{order}_fixed_normal':normal, f'q{order}_fixed_potential':potential})
        checkpoint(args, report, arrays)
        print(json.dumps(row), flush=True)
    report.update(complete=True, process_peak_rss=peak_rss_bytes(),
        limits=['One predeclared alpha=.25 interior-pole family; no placement/refinement sweep',
            'Symmetric Green integrand defines energy; raw reciprocity retained and never repaired as an evolution operator',
            'Fixed Gauss4 field/load controls are sampled diagnostics, not certified continuum error norms',
            'Finite-q body energy can have placement/rank/conditioning error despite algebraic stationarity',
            'No free-surface wave operator, new transient, nonlinear6DOF, moving-domain or performance qualification'])
    checkpoint(args, report, arrays)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('water','phase','mode','reference','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--maximum-pairs',type=int,default=50_000_000)
    parser.add_argument('--maximum-bytes',type=int,default=128*1024**2)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.maximum_pairs,args.maximum_bytes)<=0:
        parser.error('fresh output/source archive and positive budgets required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest()
        for name in ('water','phase','mode','reference')})
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
