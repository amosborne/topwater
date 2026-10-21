"""Independent canonical-hull Dirichlet energy/Neumann potential probe."""
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.linalg import solve
from threadpoolctl import threadpool_limits
from .operator import assemble_energy
from tests.reference_hull import mesh_arrays
from tests.laplace_reference import harmonic_field,point_source,triangle_samples
from topwater import HullMesh
from topwater.rankine import TrianglePanels


def run(refinement,bounded):
    vertices,faces=mesh_arrays(refinement)
    mesh=HullMesh(vertices,faces);p=TrianglePanels.from_mesh(mesh,1 if bounded else -1)
    start=perf_counter();energy,diagnostics=assemble_energy(mesh.faces,p,3)
    seconds=perf_counter()-start
    points,weights=triangle_samples(p.triangles,16)
    if bounded:
        def field(x):
            ph,g=harmonic_field(x);ps,gs=point_source(x,[3.5,1.,.7])
            return np.column_stack((ph,ps)),np.concatenate((g,gs),axis=1)
    else:
        def field(x):return point_source(x,[0,0,-.175])
    phi,_=field(vertices)
    reference_phi,gradient=field(points.reshape(-1,3))
    reference_phi=reference_phi.reshape(len(p),len(weights),-1)
    gradient=gradient.reshape(len(p),len(weights),-1,3)
    reference_q=np.einsum('fqrj,fj->fqr',gradient,p.normals)
    exact_energy=np.einsum('f,q,fqr,fqs->rs',p.areas,weights,reference_phi,reference_q)
    opposite=np.roll(p.triangles,-2,axis=1)-np.roll(p.triangles,-1,axis=1)
    shape_gradient=np.cross(p.geometric_normals[:,None,:],opposite)/(2*p.areas[:,None,None])
    shape=1/3+np.einsum('fqj,fcj->fqc',points-p.centroids[:,None],shape_gradient)
    load=np.zeros_like(phi);area=np.zeros(len(vertices))
    for corner in range(3):
        values=np.einsum('f,q,fq,fqr->fr',p.areas,weights,shape[:,:,corner],reference_q)
        np.add.at(load,faces[:,corner],values)
        np.add.at(area,faces[:,corner],p.areas/3)
    if bounded:
        matrix=np.block([[energy,area[:,None]],[area[None],np.zeros((1,1))]])
        solved=solve(matrix,np.vstack((load,np.zeros((1,phi.shape[1])))),assume_a='sym')[:-1]
        expected=phi-(area@phi)/area.sum()
    else:
        solved=solve(energy,load,assume_a='pos');expected=phi
    neumann_error=np.sqrt(area@((solved-expected)**2)/(area@(expected**2)))
    return dict(panels=len(p),bounded=bounded,seconds=seconds,**diagnostics,
        energy_relative_diagonal_error=(np.abs(np.diag(phi.T@energy@phi-exact_energy))/np.diag(exact_energy)).tolist(),
        neumann_nodal_potential_relative_rms=neumann_error.tolist(),
        net_reference_flux=load.sum(axis=0).tolist())


if __name__=='__main__':
    rows=[]
    with threadpool_limits(1):
        for bounded in (True,False):
            for refinement in (2,3,4):
                row=run(refinement,bounded);rows.append(row);print(json.dumps(row),flush=True)
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-hull-check.json').write_text(json.dumps(rows,indent=2)+'\n')
