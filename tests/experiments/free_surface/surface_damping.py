"""Consistent weak P1 surface damping used by radiation validation cases."""
import numpy as np
from scipy.sparse import coo_matrix
from tests.laplace_reference import triangle_samples


def matrices(points,faces,free_nodes,damping):
    """Return projected surface mass M and Q=integral(sigma*N_i*N_j)."""
    damping=np.asarray(damping)
    if damping.shape!=(len(free_nodes),) or not np.all(np.isfinite(damping)) or np.any(damping<0):
        raise ValueError('damping must be finite, nonnegative and indexed by free_nodes')
    lookup=np.full(len(points),-1,dtype=int)
    lookup[free_nodes]=np.arange(len(free_nodes))
    local=lookup[faces]
    if np.any(local<0):
        raise ValueError('surface faces reference non-free nodes')
    triangles=points[faces]
    normal=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0])
    area=normal[:,2]/2
    if np.any(area<=0):
        raise ValueError('surface has inverted or vertical projected triangles')
    shape,weights=triangle_samples(np.eye(3)[None],3)
    shape=shape[0]
    sigma=np.einsum('qi,fi->fq',shape,damping[local])
    elements=area[:,None,None]*(np.ones((3,3))+np.eye(3))/12
    weighted=area[:,None,None]*np.einsum('q,qi,qj,fq->fij',weights,shape,shape,sigma)
    indices=(np.repeat(local,3,axis=1).ravel(),np.tile(local,(1,3)).ravel())
    size=(len(free_nodes),len(free_nodes))
    return tuple(coo_matrix((values.ravel(),indices),shape=size).tocsr()
        for values in (elements,weighted))
