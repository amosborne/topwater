"""Contained 3-D local BEM/exterior FEM graph-domain accuracy/cost falsifier.

No boat or evolution: independent harmonic field on an actual curved graph.
The exterior supplies vertical-interface AND free-surface flux responses.
Shared potential nodes and summed weak flux are the coupling, with no blending.
"""
from itertools import permutations
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import solve
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits
from topwater import HullMesh
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples
from tests.experiments.free_surface.operator import assemble_energy
from tests.benchmark_rankine import peak_rss_bytes


def field(points):
    k=.7; kx=k*np.cos(.63); ky=k*np.sin(.63)
    phase=kx*points[...,0]+ky*points[...,1]+.37
    exponential=np.exp(k*points[...,2])
    value=exponential*np.cos(phase)
    gradient=exponential[...,None]*np.stack((-kx*np.sin(phase),
        -ky*np.sin(phase),k*np.cos(phase)),axis=-1)
    return value,gradient


def grid(n,amplitude,shift,*,exterior_scale=1.):
    if (n!=6 and (n<12 or n%12)) or not 0<=shift<n//3:
        raise ValueError('use n=6 or a positive multiple of 12 and an interior pocket')
    if not np.isfinite(exterior_scale) or exterior_scale<1:
        raise ValueError('exterior_scale must be finite and at least one')
    nx=ny=n; nz=n//3
    x=np.linspace(-2*np.pi,2*np.pi,nx+1)
    y=np.linspace(-np.pi,np.pi,ny+1)
    if exterior_scale!=1:
        # Move only exterior coordinates: local pocket and body resolution stay
        # fixed during outer-domain sensitivity checks. This is not an accuracy
        # claim for the resulting coarser exterior cells.
        for coordinates,count in ((x,nx),(y,ny)):
            first,last=count//3,2*count//3
            coordinates[:first]=coordinates[first]+exterior_scale*(coordinates[:first]-coordinates[first])
            coordinates[last+1:]=coordinates[last]+exterior_scale*(coordinates[last+1:]-coordinates[last])
    sigma=np.linspace(0,1,nz+1)
    xx,yy,ss=np.meshgrid(x,y,sigma,indexing='ij')
    eta=amplitude*(np.cos(.4*xx+.7*yy+.2)+.4*np.sin(.8*xx-.3*yy))
    points=np.column_stack((xx.ravel(),yy.ravel(),(-4+(4+eta)*ss).ravel()))
    numbering=np.arange(len(points)).reshape(xx.shape)
    # Each monotone path from one cube corner to its opposite defines a tet.
    steps=np.zeros((6,4,3),dtype=int)
    steps[:,1:]=np.cumsum(np.eye(3,dtype=int)[list(permutations(range(3)))],axis=1)
    offsets=steps@np.array([(ny+1)*(nz+1),nz+1,1])
    tets=(numbering[:-1,:-1,:-1].ravel()[:,None,None]+offsets).reshape(-1,4)
    i,j=np.indices((nx,ny))
    pocket=(nx//3+shift<=i)&(i<2*nx//3+shift)&(ny//3<=j)&(j<2*ny//3)
    near=np.repeat(pocket.ravel(),nz*6)
    return points,tets,near,numbering[:,:,-1].ravel()


def boundary(points,tets):
    """Exterior triangles of a tetrahedral union, each independently oriented."""
    faces=tets[:,[[1,2,3],[0,2,3],[0,1,3],[0,1,2]]].reshape(-1,3)
    _,first,counts=np.unique(np.sort(faces,axis=1),axis=0,
                             return_index=True,return_counts=True)
    # Keep first-encounter order as well as the original tetrahedral vertices.
    exterior=np.sort(first[counts==1])
    faces=faces[exterior]
    p=points[faces]
    # Candidate faces and tets.ravel() share tet-major/opposite-corner order.
    inward=np.einsum('ij,ij->i',np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),
                     points[tets.ravel()[exterior]]-p[:,0])>0
    faces[inward]=faces[inward,::-1]
    return faces


def stiffness(points,tets):
    coordinates=points[tets]
    edges=coordinates[:,1:]-coordinates[:,:1]
    inverse=np.linalg.inv(edges)
    gradient=np.empty((len(tets),4,3))
    gradient[:,1:]=inverse.transpose(0,2,1)
    gradient[:,0]=-gradient[:,1:].sum(axis=1)
    volume=np.abs(np.linalg.det(edges))/6
    local=np.einsum('tik,tjk,t->tij',gradient,gradient,volume)
    return coo_matrix((local.ravel(),(np.repeat(tets,4,axis=1).ravel(),
        np.tile(tets,(1,4)).ravel())),shape=(len(points),len(points))).tocsr()


def surface_data(points,faces,fs):
    panels=TrianglePanels(points[faces])
    samples,weights=triangle_samples(panels.triangles,12)
    bary,_=triangle_samples(np.eye(3)[None],12)
    _,grad=field(samples)
    exact_q=np.einsum('fqi,fi->fq',grad,panels.normals)
    load=np.zeros(len(points)); nonfree=~np.all(np.isin(faces,fs),axis=1)
    for corner in range(3):
        values=np.einsum('f,q,q,fq->f',panels.areas,weights,bary[0,:,corner],exact_q)
        np.add.at(load,faces[nonfree,corner],values[nonfree])
    return panels,samples,weights,bary,exact_q,load,~nonfree,faces


def solve_case(points,tets,near,fs,kind,physical):
    start=perf_counter()
    if kind=='pocket':
        energy=stiffness(points,tets[~near])
        near_faces=boundary(points,tets[near])
        nodes,inverse=np.unique(near_faces,return_inverse=True)
        mesh=HullMesh(points[nodes],inverse.reshape(-1,3))
        local,diagnostics=assemble_energy(mesh.faces,TrianglePanels.from_mesh(mesh),3)
        energy+=coo_matrix((local.ravel(),(np.repeat(nodes,len(nodes)),
            np.tile(nodes,len(nodes)))),shape=energy.shape).tocsr()
        near_nodes=len(nodes);near_panels=len(near_faces)
    else:
        energy=stiffness(points,tets)
        near_nodes=near_panels=0;diagnostics={}
    assembly_seconds=perf_counter()-start
    active=np.flatnonzero(energy.diagonal()>0)
    unknown=np.setdiff1d(active,fs)
    exact,_=field(points)
    panels,samples,weights,bary,exact_q,load,top,physical_faces=physical
    rhs=load[unknown]-energy[unknown][:,fs]@exact[fs]
    factor_start=perf_counter();factor=splu(energy[unknown][:,unknown].tocsc())
    factor_seconds=perf_counter()-factor_start
    factor_bytes=sum(m.data.nbytes+m.indices.nbytes+m.indptr.nbytes
                     for m in (factor.L,factor.U))
    times=[]
    for _ in range(6):
        solve_start=perf_counter();answer=factor.solve(rhs);times.append(perf_counter()-solve_start)
    potential=np.zeros(len(points));potential[fs]=exact[fs];potential[unknown]=answer
    residual=energy@potential-load
    lookup=np.full(len(points),-1);lookup[fs]=np.arange(len(fs))
    top_faces=lookup[physical_faces[top]]
    mass=np.zeros((len(fs),len(fs)))
    for face,area in zip(top_faces,panels.areas[top]):
        mass[np.ix_(face,face)]+=area*(np.ones((3,3))+np.eye(3))/12
    q=solve(mass,residual[fs],assume_a='pos')
    q_sample=np.einsum('qc,fc->fq',bary[0],q[top_faces])
    error=np.sqrt(np.einsum('f,q,fq->',panels.areas[top],weights,
        (q_sample-exact_q[top])**2)/np.einsum('f,q,fq->',panels.areas[top],weights,exact_q[top]**2))
    # Independent exact volume energy: |grad(phi)|²=k² exp(2 k z).
    k=.7
    projected=panels.areas[top]*panels.normals[top,2]
    exact_energy=np.einsum('f,q,fq->',projected,weights,
        k/2*(np.exp(2*k*samples[top,:,2])-np.exp(-8*k)))
    actual_energy=potential@(energy@potential)
    algebra=float(np.linalg.norm(residual[unknown])/np.linalg.norm(rhs))
    report=dict(kind=kind,active_nodes=len(active),unknowns=len(unknown),
        near_nodes=near_nodes,near_panels=near_panels,matrix_nonzeros=energy.nnz,
        matrix_bytes=energy.data.nbytes+energy.indices.nbytes+energy.indptr.nbytes,
        factor_bytes=factor_bytes,surface_mass_bytes=mass.nbytes,
        assembly_seconds=assembly_seconds,factor_seconds=factor_seconds,
        solve_first_seconds=times[0],solve_warm_seconds=times[1:],
        continuous_surface_flux_relative_rms=float(error),
        energy_relative_error=float(abs(actual_energy/exact_energy-1)),
        original_algebraic_relative_residual=algebra,
        assembly_factor_first_solve_seconds=assembly_seconds+factor_seconds+times[0],
        solve_and_diagnostics_seconds=perf_counter()-start,**diagnostics)
    return report,q_sample


def placement_study():
    """Move the partition by pi/3 on the unchanged physical mesh, then refine."""
    start=perf_counter()
    rows=[]
    for n in (12,24):
        points,tets,near,fs=grid(n,.04,0)
        physical_faces=boundary(points,tets)
        physical=surface_data(points,physical_faces,fs)
        base=solve_case(points,tets,near,fs,'full_fem',physical)
        first=solve_case(points,tets,near,fs,'pocket',physical)
        shifted_points,shifted_tets,shifted_near,shifted_fs=grid(n,.04,n//12)
        np.testing.assert_array_equal(shifted_points,points)
        np.testing.assert_array_equal(shifted_tets,tets)
        np.testing.assert_array_equal(shifted_fs,fs)
        second=solve_case(points,tets,shifted_near,fs,'pocket',physical)
        panels,samples,weights,bary,exact_q,load,top,faces=physical
        denominator=np.einsum('f,q,fq->',panels.areas[top],weights,exact_q[top]**2)
        def norm(values):
            return float(np.sqrt(np.einsum('f,q,fq->',panels.areas[top],weights,values**2)/denominator))
        row=dict(n=n,full_fem=base[0],central=first[0],shifted=second[0],
            partition_flux_relative_difference=norm(first[1]-second[1]),
            central_minus_fem_relative_difference=norm(first[1]-base[1]),
            shifted_minus_fem_relative_difference=norm(second[1]-base[1]))
        rows.append(row)
        print(json.dumps(row),flush=True)
    keys=('partition_flux_relative_difference','central_minus_fem_relative_difference',
          'shifted_minus_fem_relative_difference')
    passed=all(max(row[key] for key in keys)<row['full_fem']['continuous_surface_flux_relative_rms'] for row in rows)
    passed=passed and all(rows[1][key]<rows[0][key] for key in keys)
    return dict(qualified=False,interface_implementation_check_passed=bool(passed),
        results=rows,peak_rss_bytes=peak_rss_bytes(),
        complete_experiment_seconds=perf_counter()-start,
        scope='Body-free graph-domain harmonic interface check; physical flux accuracy, wave evolution, open boundary, hull and realtime remain unqualified.')


if __name__=='__main__':
    with threadpool_limits(1):
        report=placement_study()
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-pocket.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['interface_implementation_check_passed']:
        raise RuntimeError('partition dependence did not converge below the full-FEM discretization error')
