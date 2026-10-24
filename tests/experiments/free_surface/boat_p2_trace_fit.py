"""Equal-DOF continuous P2 trace fitting on original facets; no BEM solve."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import cholesky, cho_solve, eigvalsh
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import owned


def shapes(bary):
    """Vertex then (01,12,20) edge shapes on an affine input triangle."""
    return np.concatenate((bary*(2*bary-1),
        4*np.stack((bary[...,0]*bary[...,1],bary[...,1]*bary[...,2],bary[...,2]*bary[...,0]),axis=-1)),axis=-1)


def checkpoint(args,report,arrays):
    np.savez(args.output.with_suffix('.npz'),**arrays)
    report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False))


def ratios(a,b):
    return [float(x/y) if y else (0. if x==0 else None) for x,y in zip(a,b,strict=True)]


def run(args,report):
    projection,gal,family=(owned(getattr(args,name)) for name in ('projection','galerkin','family'))
    if (projection['parents']['galerkin']!=sha256(args.galerkin.read_bytes()).hexdigest()
            or projection['parents']['family']!=sha256(args.family.read_bytes()).hexdigest()
            or gal['parents']['family']!=sha256(args.family.read_bytes()).hexdigest()):
        raise ValueError('equal-DOF fitting parents differ')
    for name in ('tests/laplace_reference.py','src/topwater/rankine.py','src/topwater/constants.py'):
        relative=Path(name)
        if sha256(relative.read_bytes()).hexdigest()!=projection['producer_sha256'].get(str(relative)):
            raise ValueError('owned sample/geometry reconstruction source differs: '+name)
    with np.load(args.family.with_suffix('.npz')) as z:
        original=z['triangles'].copy()
    with np.load(args.galerkin.with_suffix('.npz')) as z:
        fine,vertices,faces,lift,center=(z[name].copy() for name in
            ('triangles','vertices','faces','odd_lift','world_cg'))
    with np.load(args.projection.with_suffix('.npz')) as z:
        np.testing.assert_array_equal(fine,z['physical_triangles'])
        np.testing.assert_array_equal(original,z['source_triangles'])
        np.testing.assert_array_equal(center,z['world_cg'])
        data={f'q{q}_{name}':z[f'q{q}_{name}'].copy() for q in (4,8) for name in
            ('reference_phi','projection_port','reference_port')}
    nb,nf,n=len(original),len(fine),lift.shape[1]
    if nf!=4*nb:
        raise ValueError('owned one-midpoint-refinement targets required')
    np.testing.assert_array_equal(vertices[faces[:nf]],fine)
    # Exact coordinate identities with the ALREADY existing fine nodes. Shared
    # edges share IDs, and every z=0 vertex/midpoint has zero saved odd lift.
    lookup={tuple(point):i for i,point in enumerate(vertices)}
    edge=.5*(original+np.roll(original,-1,axis=1))
    nodes=np.concatenate((original,edge),axis=1)
    node_ids=np.array([[lookup[tuple(point)] for point in tri] for tri in nodes])
    node_lift=lift[node_ids]
    if np.any(node_lift[nodes[:,:,2]==0]) or np.any((node_lift!=0)&(node_lift!=1)):
        raise ValueError('P2 perimeter grounding or lower-body node identity differs')
    if np.count_nonzero(np.any(node_lift,axis=(0,1)))!=n:
        raise ValueError('P2 mapped nodes do not span the SAME saved free coefficient IDs')
    nodal_bary=np.concatenate((np.eye(3),.5*(np.eye(3)+np.roll(np.eye(3),-1,axis=0))))
    np.testing.assert_array_equal(shapes(nodal_bary),np.eye(6))
    child=np.array([[[1,0,0],[.5,.5,0],[.5,0,.5]],
        [[.5,.5,0],[0,1,0],[0,.5,.5]],
        [[.5,0,.5],[0,.5,.5],[0,0,1]],
        [[.5,.5,0],[0,.5,.5],[.5,0,.5]]])
    child=np.repeat(child,nb,axis=0)
    parent=np.tile(np.arange(nb),4)
    panels=TrianglePanels(fine,-1) # Geometry only; no influence evaluation.
    memory=8*(nf*8**2*(2*n+48)+12*n*n+nf*6*n)+32*1024**2
    if memory>args.maximum_bytes:
        raise MemoryError('P2 fitting array/workspace preflight exceeds budget')
    arrays=dict(original_triangles=original,physical_triangles=fine,world_cg=center,
        p2_node_ids=node_ids,p2_node_coordinates=nodes)
    report.update(complete=False,qualified=False,original_facets=nb,target_facets=nf,
        shared_free_trace_dofs=n,original_vertex_and_midpoint_nodes=len(np.unique(node_ids)),
        grounded_seam_nodes=int(np.sum(vertices[np.unique(node_ids),2]==0)),
        array_workspace_bound_bytes=memory,
        memory_scope='Array/workspace estimate excludes Python/library baseline; no fluid influence kernel runs',orders=[])
    checkpoint(args,report,arrays)
    mass=factor=scale=B=None
    previous=None
    for order in (4,8):
        started=perf_counter()
        bary,weight=triangle_samples(np.eye(3)[None],order);bary=bary[0]
        lam=np.einsum('qj,fjk->fqk',bary,child)
        N=np.einsum('fql,flc->fqc',shapes(lam),node_lift[parent]).reshape(-1,n)
        points=np.einsum('qj,fji->fqi',bary,fine).reshape(-1,3)
        original_points=np.einsum('fqj,fji->fqi',lam,original[parent]).reshape(-1,3)
        allowance=64*np.finfo(float).eps*np.einsum('fqj,fji->fqi',lam,abs(original[parent])).reshape(-1,3)
        if np.any(abs(points-original_points)>allowance):
            raise ValueError('dyadic original barycentric sample correspondence differs')
        w=(panels.areas[:,None]*weight).ravel()
        reference=data[f'q{order}_reference_phi']
        if reference.shape!=(len(w),6):
            raise ValueError('saved reference field/sample ordering differs')
        norm2=(w[:,None]*reference**2).sum(axis=0)
        if np.any(norm2<=0):
            raise ValueError('six-column reference normalization unresolved; no mode exclusion')
        normal=np.repeat(panels.normals,len(weight),axis=0)
        T=np.column_stack((normal,np.cross(points-center,normal)))
        current_mass=N.T@(w[:,None]*N)
        current_B=N.T@(w[:,None]*T)
        arithmetic=128*np.finfo(float).eps*(len(w)+n)
        if mass is None:
            mass=current_mass.copy();B=current_B.copy()
            diagonal=mass.diagonal()
            if np.any(diagonal<=0):
                raise ValueError('zero P2 fitting support; no mass invention')
            scale=np.sqrt(diagonal)
            scaled=mass/scale[:,None]/scale[None,:]
            eigen=eigvalsh(scaled)
            threshold=n*np.finfo(float).eps*max(abs(eigen))
            if eigen[0]<=threshold:
                raise ValueError('P2 fitting mass lacks resolved full rank; no jitter')
            factor=cholesky(scaled,lower=True)
            report.update(scaled_mass_eigenvalues=eigen.tolist(),mass_rank=n,
                rank_threshold=float(threshold),scaled_mass_condition=float(eigen[-1]/eigen[0]))
            arrays.update(p2_fitting_mass=mass,p2_polynomial_load=B)
        mass_error=float(abs(current_mass-mass).max());load_error=float(abs(current_B-B).max())
        if mass_error>arithmetic*abs(mass).max() or load_error>arithmetic*abs(B).max():
            raise ValueError('q4/q8 degree-four mass or degree-three load polynomial parity failed')
        rhs=N.T@(w[:,None]*reference)
        coefficients=cho_solve((factor,True),rhs/scale[:,None])/scale[:,None]
        fitted=N@coefficients
        error=fitted-reference
        error2=(w[:,None]*error**2).sum(axis=0)
        fitted2=(w[:,None]*fitted**2).sum(axis=0)
        floor=np.sqrt(error2/norm2)
        p1row=next(row for row in projection['orders'] if row['order']==order)
        p1floor=np.array(p1row['relative_trace_L2']['projection_minus_reference'])
        orth=N.T@(w[:,None]*error)
        pyth=norm2-fitted2-error2
        orth_allowance=arithmetic*float((abs(N).T@(w[:,None]*(abs(fitted)+abs(reference)))).max())
        pyth_allowance=arithmetic*float((norm2+fitted2+error2).max())
        if abs(orth).max()>orth_allowance or abs(pyth).max()>pyth_allowance:
            raise ValueError('P2 projection orthogonality/Pythagorean arithmetic gate failed')
        port=WATER_DENSITY*T.T@(w[:,None]*fitted)
        port_parity=float(abs(port-WATER_DENSITY*B.T@coefficients).max())
        port_allowance=arithmetic*float((WATER_DENSITY*abs(T).T@(w[:,None]*abs(fitted))).max())
        if port_parity>port_allowance:
            raise ValueError('P2 polynomial load-port reconstruction failed')
        reference_port=data[f'q{order}_reference_port']
        p1port=data[f'q{order}_projection_port']
        row=dict(order=order,seconds=perf_counter()-started,relative_P2_L2_floor=floor.tolist(),
            relative_P1_L2_floor=p1floor.tolist(),P2_over_P1_floor=ratios(floor,p1floor),
            mass_q4_parity_max=mass_error,polynomial_load_q4_parity_max=load_error,
            projection_orthogonality_max=float(abs(orth).max()),pythagorean_defect=pyth.tolist(),
            polynomial_port_parity_max=port_parity,
            P2_minus_reference_diagonal_load=(port-reference_port).diagonal().tolist(),
            P1_minus_reference_diagonal_load=(p1port-reference_port).diagonal().tolist(),
            P2_minus_P1_diagonal_load=(port-p1port).diagonal().tolist())
        if previous is not None:
            row.update(q4_q8_relative_floor_change=ratios(floor-previous[0],previous[0]),
                q4_q8_port_diagonal_change=(port-previous[1]).diagonal().tolist(),
                q4_q8_fitted_field_relative_L2_change=ratios(np.sqrt((w[:,None]*(N@(coefficients-previous[2]))**2).sum(axis=0)),np.sqrt(norm2)))
        previous=(floor,port,coefficients)
        arrays.update({f'q{order}_p2_coefficients':coefficients,f'q{order}_p2_port':port,
            f'q{order}_p2_minus_reference_port':port-reference_port,f'q{order}_p2_minus_p1_port':port-p1port})
        report['orders'].append(row);checkpoint(args,report,arrays);print(json.dumps(row),flush=True)
    report.update(complete=True,process_peak_rss=peak_rss_bytes(),limits=[
        'Cached coarse P0 source field is not converged physical truth and may contain original-triangulation edge artifacts',
        'P2 original-facet and P1 refined-facet spaces have equal free DOFs but are different spaces; no nesting or Galerkin accuracy implication',
        'Zero seam values are exact vertex/midpoint constraints, no tolerance welding or new geometry',
        'Fitting mass is not water mass, and no six-load fitting or BEM energy/kernel is used',
        'Equal DOF does not mean equal operator/quadrature cost or stable auxiliary-flux pairing',
        'A fitting improvement only motivates a bounded P2 kernel test; no moving-domain, nonlinear6DOF, motion or speed qualification'])
    checkpoint(args,report,arrays)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('projection','galerkin','family','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--maximum-bytes',type=int,default=128*1024**2)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or args.maximum_bytes<=0:
        parser.error('fresh source/output and positive array bound required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True);hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest()
        for name in ('projection','galerkin','family')});started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
