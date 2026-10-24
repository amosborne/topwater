"""L2 trace-space discriminator using a cached finite coarse source control."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
from types import SimpleNamespace
import numpy as np
from scipy.linalg import cholesky, cho_solve, eigvalsh
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from topwater.surface import body_load
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import owned
from .boat_spatial_work import subdivide
from .owned_odd_body_trace import owned_odd_maps


def checkpoint(args,report,arrays):
    np.savez(args.output.with_suffix('.npz'),**arrays)
    report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False))


def ratios(a,b):
    return [float(x/y) if y else (0. if x==0 else None) for x,y in zip(a,b,strict=True)]


def run(args,report):
    gal,family,reference,odd=(owned(getattr(args,name)) for name in ('galerkin','family','reference','odd'))
    if (gal['parents']['family']!=sha256(args.family.read_bytes()).hexdigest()
            or gal['parents']['reference']!=sha256(args.reference.read_bytes()).hexdigest()
            or family['parents']['reference']!=sha256(args.reference.read_bytes()).hexdigest()
            or reference['parents']['odd']!=sha256(args.odd.read_bytes()).hexdigest()):
        raise ValueError('projection physical parents differ')
    if gal['refinement_level']!=1 or not gal['body_orders'][-1]['accepted'] or gal['body_orders'][-1]['order']!=8:
        raise ValueError('owned single-refined q8 Galerkin result required')
    for name in ('src/topwater/rankine.py','src/topwater/_rankine_kernel.py','src/topwater/_rankine_extended.py',
                 'src/topwater/numerics.py','src/topwater/constants.py',
                 'tests/experiments/free_surface/owned_odd_body_trace.py',
                 'tests/experiments/free_surface/panel_gradient_reference.py',
                 'tests/experiments/free_surface/boat_spatial_work.py'):
        relative=Path(name)
        if sha256(relative.read_bytes()).hexdigest()!=reference['producer_sha256'].get(str(relative)):
            raise ValueError('coarse source-control owner changed: '+name)
    with np.load(args.reference.with_suffix('.npz')) as z:
        coarse_tri,center,sigma,capy=(z[name].copy() for name in
            ('triangles','world_cg','sigma_collocation','capytaine_added_mass'))
    with np.load(args.family.with_suffix('.npz')) as z:
        np.testing.assert_array_equal(coarse_tri,z['triangles'])
        np.testing.assert_array_equal(center,z['world_cg'])
    with np.load(args.galerkin.with_suffix('.npz')) as z:
        fine,vertices,faces,lift,response,B,gal_mass,S=(z[name].copy() for name in
            ('triangles','vertices','faces','odd_lift','q8_odd_response','odd_lower_load','q8_added_mass','q8_half_energy'))
        np.testing.assert_array_equal(center,z['world_cg'])
    np.testing.assert_array_equal(fine,subdivide(coarse_tri))
    np.testing.assert_array_equal(vertices[faces[:len(fine)]],fine)
    coarse,fine_panels=TrianglePanels(coarse_tri,-1),TrianglePanels(fine,-1)
    nb,nf,n=len(coarse),len(fine),lift.shape[1]
    if S.shape!=(n,n) or not np.isfinite(S).all():
        raise ValueError('owned fine q8 half-energy dimensions/values differ')
    known_factor=cholesky(S,lower=True)
    parent=np.tile(np.arange(nb),4)
    with np.load(args.odd.with_suffix('.npz')) as z:
        jd=z['jd'].copy()
    modes=np.column_stack((coarse.normals,np.cross(coarse.centroids-center,coarse.normals)))
    residual=jd@sigma-modes
    source_allowance=128*np.finfo(float).eps*nb*(abs(jd)@abs(sigma)+abs(modes)).max()
    if abs(residual).max()>source_allowance:
        raise ValueError('cached source strengths do not satisfy their owned coarse collocation matrix')
    pair_bound=2*nb*nf*(4**2+8**2)
    memory=8*(nf*8**2*(2*nb+n+96)+2*len(vertices)**2+8*n*n+1024*nb*8)+32*1024**2
    if pair_bound>args.maximum_pairs or memory>args.maximum_bytes:
        raise RuntimeError('projection fixed q4/q8 query/workspace preflight exceeds budget')
    # A lower-body P1 consistent mass, excluding seam DOFs through the SAME
    # saved odd lift. This is a trace fitting mass, never physical water mass.
    full_mass=np.zeros((len(vertices),len(vertices)))
    local=(np.eye(3)+np.ones((3,3)))/12
    for face,area in zip(faces[:nf],fine_panels.areas,strict=True):
        full_mass[np.ix_(face,face)]+=area*local
    mass=lift.T@full_mass@lift
    diagonal=mass.diagonal()
    if np.any(diagonal<=0):
        raise ValueError('trace projection has a zero-support unknown; no fitting-mass invention')
    scale=np.sqrt(diagonal)
    scaled=mass/scale[:,None]/scale[None,:]
    eigenvalues=eigvalsh(scaled)
    rank_threshold=n*np.finfo(float).eps*max(abs(eigenvalues))
    if eigenvalues[0]<=rank_threshold:
        raise ValueError('trace projection mass lacks resolved full rank; no jitter')
    factor=cholesky(scaled,lower=True)
    independent_B=lift.T@body_load(SimpleNamespace(vertices=vertices,faces=faces),fine_panels,np.arange(nf),center)
    np.testing.assert_allclose(independent_B,B,rtol=1e-12,atol=1e-13)
    lower_lift=lift[faces[:nf]]
    if np.any((lower_lift!=0)&(lower_lift!=1)):
        raise ValueError('lower-body odd-lift basis must be nonnegative with zero seam trace')
    child=np.array([[[1,0,0],[.5,.5,0],[.5,0,.5]],
        [[.5,.5,0],[0,1,0],[0,.5,.5]],
        [[.5,0,.5],[0,.5,.5],[0,0,1]],
        [[.5,.5,0],[0,.5,.5],[.5,0,.5]]])
    child=np.repeat(child,nb,axis=0)
    report.update(complete=False,qualified=False,coarse_source_panels=nb,physical_target_faces=nf,
        trace_unknowns=n,source_pair_bound=pair_bound,array_workspace_bound_bytes=memory,
        memory_scope='Array/workspace bound excludes Python/JIT/library baseline, not a process-RSS cap',
        source_collocation_residual_max=float(abs(residual).max()),
        source_collocation_arithmetic_allowance=float(source_allowance),
        scaled_projection_mass_eigenvalues=eigenvalues.tolist(),
        projection_mass_rank_threshold=float(rank_threshold),projection_mass_rank=n,
        source_owner_vs_fine_normal_max=float(abs(coarse.normals[parent]-fine_panels.normals).max()),
        reference='Cached coarse original-panel P0 odd-source reconstruction, NOT converged Capytaine truth',orders=[])
    arrays=dict(source_triangles=coarse_tri,physical_triangles=fine,world_cg=center,
        source_sigma=sigma,projection_mass=mass,galerkin_trace_coefficients=response,
        exact_weak_load=B,cached_capytaine_centroid_matrix=capy,known_half_energy=S)
    checkpoint(args,report,arrays)
    for order in (4,8):
        started=perf_counter()
        bary,weights=triangle_samples(np.eye(3)[None],order);bary=bary[0]
        points=np.einsum('qj,fji->fqi',bary,fine).reshape(-1,3)
        original_bary=np.einsum('qj,fjk->fqk',bary,child).reshape(-1,3)
        owners=np.repeat(parent,len(weights))
        R,J,diagnostics=owned_odd_maps(coarse,points,owners,original_bary)
        reference_phi,reference_normal=R@sigma,J@sigma
        w=(fine_panels.areas[:,None]*weights).ravel()
        N=np.einsum('qj,fjc->fqc',bary,lower_lift).reshape(-1,n)
        rhs=N.T@(w[:,None]*reference_phi)
        arithmetic_factor=128*np.finfo(float).eps*(len(w)+n)
        mass_parity=float(abs(N.T@(w[:,None]*N)-mass).max())
        mass_allowance=arithmetic_factor*float(abs(mass).max())
        if mass_parity>mass_allowance:
            raise ValueError('positive quadrature basis differs from the exact trace fitting mass')
        coefficients=cho_solve((factor,True),rhs/scale[:,None])/scale[:,None]
        phi_gal,phi_projection=N@response,N@coefficients
        known_rhs=N.T@(w[:,None]*reference_normal)
        known_coefficients=cho_solve((known_factor,True),known_rhs)
        phi_known=N@known_coefficients
        known_residual=S@known_coefficients-known_rhs
        total=phi_gal-reference_phi
        projection=phi_projection-reference_phi
        gal_projection=phi_gal-phi_projection
        reference_norm2=(w[:,None]*reference_phi**2).sum(axis=0)
        norm2=[(w[:,None]*value**2).sum(axis=0) for value in (total,projection,gal_projection)]
        orthogonality=N.T@(w[:,None]*projection)
        pythagorean=norm2[0]-norm2[1]-norm2[2]
        orth_allowance=arithmetic_factor*float((abs(N).T@(w[:,None]*(abs(phi_projection)+abs(reference_phi)))).max())
        pyth_allowance=arithmetic_factor*float((norm2[0]+norm2[1]+norm2[2]+reference_norm2).max())
        if abs(orthogonality).max()>orth_allowance or abs(pythagorean).max()>pyth_allowance:
            raise ValueError('quadrature L2 projection ownership/orthogonality gate failed')
        source_normal=coarse.normals[owners]
        rigid_source=np.column_stack((source_normal,np.cross(points-center,source_normal)))
        normal_error=reference_normal-rigid_source
        normal_rms=np.sqrt((w[:,None]*normal_error**2).sum(axis=0)/w.sum())
        normal_scale=np.sqrt((w[:,None]*rigid_source**2).sum(axis=0)/w.sum())
        physical_normal=np.repeat(fine_panels.normals,len(weights),axis=0)
        T=np.column_stack((physical_normal,np.cross(points-center,physical_normal)))
        ports=[WATER_DENSITY*T.T@(w[:,None]*phi) for phi in (phi_gal,phi_projection,reference_phi)]
        known_port=WATER_DENSITY*T.T@(w[:,None]*phi_known)
        gal_proj_port,proj_ref_port,ref_capy_port=ports[0]-ports[1],ports[1]-ports[2],ports[2]-capy
        port_parity=max(float(abs(ports[0]-gal_mass).max()),
                        float(abs(ports[1]-WATER_DENSITY*B.T@coefficients).max()),
                        float(abs(known_port-WATER_DENSITY*B.T@known_coefficients).max()))
        port_allowance=arithmetic_factor*float((WATER_DENSITY*abs(T).T@(w[:,None]*(abs(phi_gal)+abs(phi_projection)+abs(phi_known)))).max())
        if port_parity>port_allowance:
            raise ValueError('polynomial trace load reconstruction gate failed')
        row=dict(order=order,points=len(points),seconds=perf_counter()-started,source_diagnostics=diagnostics,
            relative_trace_L2=dict(galerkin_minus_reference=ratios(np.sqrt(norm2[0]),np.sqrt(reference_norm2)),
                projection_minus_reference=ratios(np.sqrt(norm2[1]),np.sqrt(reference_norm2)),
                galerkin_minus_projection=ratios(np.sqrt(norm2[2]),np.sqrt(reference_norm2))),
            projection_orthogonality_max=float(abs(orthogonality).max()),
            pythagorean_squared_norm_defect=pythagorean.tolist(),
            arithmetic_ownership_gates=dict(mass_parity=mass_parity,mass_allowance=mass_allowance,
                orthogonality_allowance=orth_allowance,pythagorean_allowance=pyth_allowance,
                polynomial_port_allowance=port_allowance),
            reference_normal_relative_rms=ratios(normal_rms,normal_scale),
            reference_normal_max=abs(normal_error).max(axis=0).tolist(),
            known_harmonic_field_consistency=dict(
                scope='Same reconstructed harmonic coarse field with its OWN normal flux, not a rigid-motion/added-mass validation',
                relative_trace_L2_minus_reference=ratios(np.sqrt((w[:,None]*(phi_known-reference_phi)**2).sum(axis=0)),np.sqrt(reference_norm2)),
                relative_trace_L2_minus_projection=ratios(np.sqrt((w[:,None]*(phi_known-phi_projection)**2).sum(axis=0)),np.sqrt(reference_norm2)),
                energy_solve_residual_max=float(abs(known_residual).max()),
                energy_solve_relative_residual=[float(np.linalg.norm(known_residual[:,j])/(np.linalg.norm(S@known_coefficients[:,j])+np.linalg.norm(known_rhs[:,j])))
                    if np.linalg.norm(S@known_coefficients[:,j])+np.linalg.norm(known_rhs[:,j]) else None for j in range(6)],
                diagonal_port_minus_reference=(known_port-ports[2]).diagonal().tolist(),
                diagonal_port_minus_projection=(known_port-ports[1]).diagonal().tolist()),
            polynomial_galerkin_load_reconstruction_max=float(abs(ports[0]-gal_mass).max()),
            polynomial_projection_load_reconstruction_max=float(abs(ports[1]-WATER_DENSITY*B.T@coefficients).max()),
            load_decomposition_max_closure=float(abs(gal_proj_port+proj_ref_port+ref_capy_port-(ports[0]-capy)).max()),
            diagonal_load_contributions=dict(galerkin_minus_projection=gal_proj_port.diagonal().tolist(),
                projection_minus_reference=proj_ref_port.diagonal().tolist(),
                reference_integrated_minus_cached_capy_centroid=ref_capy_port.diagonal().tolist()))
        arrays.update({f'q{order}_projection_coefficients':coefficients,
            f'q{order}_reference_phi':reference_phi,f'q{order}_reference_normal':reference_normal,
            f'q{order}_galerkin_minus_projection_port':gal_proj_port,
            f'q{order}_projection_minus_reference_port':proj_ref_port,
            f'q{order}_reference_minus_capy_port':ref_capy_port,
            f'q{order}_galerkin_port':ports[0],f'q{order}_projection_port':ports[1],f'q{order}_reference_port':ports[2]})
        arrays.update({f'q{order}_known_flux_rhs':known_rhs,f'q{order}_known_trace_coefficients':known_coefficients,
            f'q{order}_known_port':known_port,f'q{order}_known_minus_reference_port':known_port-ports[2],
            f'q{order}_known_minus_projection_port':known_port-ports[1]})
        report['orders'].append(row)
        checkpoint(args,report,arrays)
        print(json.dumps(row),flush=True)
    rhs_change=arrays['q8_known_flux_rhs']-arrays['q4_known_flux_rhs']
    report.update(complete=True,process_peak_rss=peak_rss_bytes(),
        known_rhs_q4_q8_relative_change=[float(np.linalg.norm(rhs_change[:,j])/np.linalg.norm(arrays['q8_known_flux_rhs'][:,j]))
            if np.linalg.norm(arrays['q8_known_flux_rhs'][:,j]) else None for j in range(6)],limits=[
        'The reference is finite coarse source-collocation data; integrated loads expose its difference from cached centroid Capytaine values',
        'L2 projection minimizes only its sampled trace norm, not six loads or the Ritz energy norm',
        'Galerkin minus L2 projection is not a pure operator error; projection minus reference is not an unavoidable load floor',
        'Orthogonality/Pythagoras apply to the declared quadrature projection, not the continuum problem',
        'Original source normals own reference Neumann checks; fine physical normals own load ports',
        'Known-flux consistency uses only the SAME loaded fine q8 energy and source field; it is not rigid accuracy and q4/q8 flux sampling remains approximate',
        'No new body solve, energy quadrature, mesh refinement, placement, moving domain, nonlinear6DOF or trajectory qualification'])
    checkpoint(args,report,arrays)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('galerkin','family','reference','odd','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--maximum-pairs',type=int,default=4_000_000)
    parser.add_argument('--maximum-bytes',type=int,default=128*1024**2)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.maximum_pairs,args.maximum_bytes)<=0:
        parser.error('fresh output/source and positive budgets required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest()
        for name in ('galerkin','family','reference','odd')})
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
