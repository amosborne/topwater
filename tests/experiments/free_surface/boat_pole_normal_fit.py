"""One weighted normal least-squares discriminator on the SAME saved poles."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import svd
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import owned
from .boat_interior_pole_inertia import point_maps


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def ratio(a, b):
    return float(a/b) if b else (0. if a == 0 else None)


def run(args, report):
    family, reference = owned(args.family), owned(args.reference)
    if family['parents']['reference'] != sha256(args.reference.read_bytes()).hexdigest():
        raise ValueError('physical six-port reference ownership differs')
    for name in ('water', 'phase', 'mode'):
        if family['parents'][name] != reference['parents'][name]:
            raise ValueError('fixed pole/reference physical parents differ')
    for name in ('tests/experiments/free_surface/boat_interior_pole_inertia.py',
                 'src/topwater/rankine.py', 'src/topwater/numerics.py', 'src/topwater/constants.py'):
        relative = Path(name)
        if sha256(relative.read_bytes()).hexdigest() != family['producer_sha256'].get(str(relative)):
            raise ValueError('saved pole geometry/kernel owner changed: '+name)
    with np.load(args.family.with_suffix('.npz')) as z:
        data = {name:z[name].copy() for name in ('poles','triangles','world_cg','fixed_q4_points',
            'fixed_q4_weights','fixed_q4_T','q64_sigma','q64_raw_RR','q64_BR','q64_fixed_normal')}
    with np.load(args.reference.with_suffix('.npz')) as z:
        capy, direct = z['capytaine_added_mass'].copy(), z['direct_added_mass'].copy()
        np.testing.assert_array_equal(z['triangles'], data['triangles'])
        np.testing.assert_array_equal(z['world_cg'], data['world_cg'])
    poles, center = data['poles'], data['world_cg']
    panels = TrianglePanels(data['triangles'], -1)
    b = len(poles)
    if b != len(panels):
        raise ValueError('saved pole-per-input-facet ownership differs')
    pair_bound = 2*b*b*(4**2+8**2)
    array_bound = 8*(17*b*b*(4**2+8**2)+12*b*b+32*b*(4**2+8**2))+16*1024**2
    if pair_bound > args.maximum_pairs or array_bound > args.maximum_bytes:
        raise RuntimeError('normal-fit source/workspace preflight exceeds budget')
    points4, weights4 = triangle_samples(panels.triangles, 4)
    points4 = points4.reshape(-1, 3)
    owners4 = np.repeat(np.arange(b), len(weights4))
    weights4 = (panels.areas[:,None]*weights4).ravel()
    normals4 = panels.normals[owners4]
    T4 = np.column_stack((normals4, np.cross(points4-center, normals4)))
    for actual, name in ((points4,'fixed_q4_points'),(weights4,'fixed_q4_weights'),(T4,'fixed_q4_T')):
        np.testing.assert_array_equal(actual, data[name])
    R4, J4 = point_maps(points4, normals4, poles)
    np.testing.assert_allclose(J4@data['q64_sigma'], data['q64_fixed_normal'], rtol=1e-12, atol=1e-12)
    weighted = np.sqrt(weights4)[:,None]*J4
    target = np.sqrt(weights4)[:,None]*T4
    U, singular, Vh = svd(weighted, full_matrices=False, check_finite=True)
    threshold = b*np.finfo(float).eps*singular[0]
    rank = int((singular > threshold).sum())
    arrays = dict(poles=poles, triangles=data['triangles'], world_cg=center,
        singular_values=singular, raw_RR64=data['q64_raw_RR'], BR64=data['q64_BR'])
    report.update(complete=False, qualified=False, source_pair_bound=pair_bound,
        array_workspace_bound_bytes=array_bound, pole_count=b, rank=rank,
        rank_threshold=float(threshold), singular_values=singular.tolist(),
        singular_rcond=ratio(singular[-1],singular[0]),
        svd_scope='Area-weighted OLD Gauss4 normals; all source singular directions must be resolved, no ridge/truncation')
    checkpoint(args,report,arrays)
    if rank != b:
        raise ValueError('normal-fit source matrix lacks resolved full rank; no rank manipulation')
    sigma_ls = Vh.T@((U.T@target)/singular[:,None])
    if not np.isfinite(sigma_ls).all():
        raise FloatingPointError('nonfinite unregularized least-squares strengths')
    points8, weights8 = triangle_samples(panels.triangles,8)
    points8 = points8.reshape(-1,3)
    owners8 = np.repeat(np.arange(b),len(weights8))
    weights8 = (panels.areas[:,None]*weights8).ravel()
    normals8 = panels.normals[owners8]
    T8 = np.column_stack((normals8,np.cross(points8-center,normals8)))
    R8,J8 = point_maps(points8,normals8,poles)
    E = .5*(data['q64_raw_RR']+data['q64_raw_RR'].T)
    BR = data['q64_BR']
    report['fields'] = {}
    for name,sigma in (('least_squares',sigma_ls),('ritz',data['q64_sigma'])):
        row = {}
        for order,R,J,T,w in ((4,R4,J4,T4,weights4),(8,R8,J8,T8,weights8)):
            potential,normal = R@sigma,J@sigma
            residual = normal-T
            rms = np.sqrt((w[:,None]*residual**2).sum(axis=0)/w.sum())
            base = np.sqrt((w[:,None]*T**2).sum(axis=0)/w.sum())
            pressure = WATER_DENSITY*T.T@(w[:,None]*potential)
            row[f'q{order}'] = dict(normal_relative_rms=[ratio(a,b) for a,b in zip(rms,base,strict=True)],
                normal_rms=rms.tolist(), normal_max=abs(residual).max(axis=0).tolist(),
                raw_pressure_port_matrix=pressure.tolist(),
                raw_pressure_port_diagonal_relative_to_capy=[ratio(a-b,b) for a,b in zip(pressure.diagonal(),capy.diagonal(),strict=True)],
                raw_pressure_port_diagonal_relative_to_direct=[ratio(a-b,b) for a,b in zip(pressure.diagonal(),direct.diagonal(),strict=True)],
                raw_pressure_port_asymmetry_max=float(abs(pressure-pressure.T).max()))
            arrays.update({f'{name}_q{order}_potential':potential,f'{name}_q{order}_normal':normal,
                           f'{name}_q{order}_pressure_port':pressure})
        energy = .5*np.sum(sigma*(E@sigma),axis=0)
        forcing = np.diag(BR@sigma)
        stationarity = E@sigma-BR.T
        row.update(saved_q64_energy_per_density=energy.tolist(), saved_q64_forcing_per_density=forcing.tolist(),
            saved_q64_mixed_functional_per_density=(energy-forcing).tolist(),
            saved_q64_stationarity_max=abs(stationarity).max(axis=0).tolist(),
            saved_q64_stationarity_relative=[ratio(np.linalg.norm(stationarity[:,j]),
                np.linalg.norm(E@sigma[:,j])+np.linalg.norm(BR.T[:,j])) for j in range(6)],
            source_strength_l2=np.linalg.norm(sigma,axis=0).tolist())
        report['fields'][name] = row
        arrays[name+'_sigma'] = sigma
    arrays.update(capytaine_reference=capy,direct_reference=direct)
    report.update(complete=True, source_pairs=pair_bound, process_peak_rss=peak_rss_bytes(),
        limits=['One SAME saved pole set; no placement or energy-quadrature sweep',
            'q4 fit and q8 heldout controls can both miss narrow peaks; neither is a continuum certificate',
            'Raw pressure-port matrices are nonsymmetric sampled diagnostics, not the Ritz added-mass energy',
            'Saved q64 scalar energy is not reevaluated or repaired; forcing/stationarity are distinct',
            'No wave operator, transient, moving-domain, nonlinear6DOF or speed qualification'])
    checkpoint(args,report,arrays)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('family','reference','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--maximum-pairs',type=int,default=1_000_000)
    parser.add_argument('--maximum-bytes',type=int,default=128*1024**2)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.maximum_pairs,args.maximum_bytes)<=0:
        parser.error('fresh source/output and positive budgets required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest()
        for name in ('family','reference')})
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
