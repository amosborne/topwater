"""Nonlinear graph rates on the canoe's moving, body-fitted surface mesh."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from tests.laplace_reference import triangle_samples


def rates(points,faces,free_nodes,psi,flux,height_jacobians,body_jacobian,velocity,*,damping=None):
    """Return ALE height/potential rates and nodal mesh velocity.

    Each free node owns its height coordinate. Horizontal mesh motion may
    depend on several heights: J_height*height_rate + J_body*world_twist.
    The three sparse height maps include the redistributed wet hull as well
    as free nodes. The free-node z map must be identity, with no body z term.
    ``flux`` is the free-surface weak load after removing prescribed hull flux.
    """
    lookup=np.full(len(points),-1,dtype=int)
    lookup[free_nodes]=np.arange(len(free_nodes))
    local=lookup[faces]
    if np.any(local<0):
        raise ValueError('free-surface faces reference non-free nodes')
    triangles=points[faces]
    edges=triangles[:,1:,:2]-triangles[:,:1,:2]
    determinant=edges[:,0,0]*edges[:,1,1]-edges[:,0,1]*edges[:,1,0]
    if np.any(determinant<=0):
        raise ValueError('ALE surface has inverted or vertical projected triangles')
    area=determinant/2
    opposite=np.roll(triangles[:,:,:2],-2,axis=1)-np.roll(triangles[:,:,:2],-1,axis=1)
    gradient=np.stack((-opposite[:,:,1],opposite[:,:,0]),axis=-1)/determinant[:,None,None]
    slope=np.einsum('fci,fc->fi',gradient,triangles[:,:,2])
    potential_gradient=np.einsum('fci,fc->fi',gradient,psi[local])
    unit_mass=(np.ones((3,3))+np.eye(3))/12
    element_mass=area[:,None,None]*unit_mass
    shape,weights=triangle_samples(np.eye(3)[None],3)
    shape=shape[0]
    z=np.einsum('qc,fc->fq',shape,triangles[:,:,2])
    if damping is None:
        sigma=np.zeros_like(z)
    else:
        damping=np.asarray(damping)
        if damping.shape!=(len(free_nodes),) or not np.all(np.isfinite(damping)) or np.any(damping<0):
            raise ValueError('damping must be finite, nonnegative and indexed by free_nodes')
        sigma=np.einsum('qc,fc->fq',shape,damping[local])
    indices=(np.repeat(local,3,axis=1).ravel(),np.tile(local,(1,3)).ravel())
    size=(len(free_nodes),len(free_nodes))
    mass=coo_matrix((element_mass.ravel(),indices),shape=size).tocsr()
    transport=[coo_matrix(((element_mass*slope[:,axis,None,None]).ravel(),indices),
        shape=size).tocsr() for axis in (0,1)]
    body_motion=np.einsum('nij,j->ni',body_jacobian,velocity)
    right=flux.copy()
    moving_mass=mass.copy()
    for axis in (0,1):
        moving_mass-=transport[axis]@height_jacobians[axis][free_nodes]
        right+=transport[axis]@body_motion[free_nodes,axis]
    for i in range(3):
        np.add.at(right,local[:,i],-np.einsum('f,q,q,fq->f',area,weights,shape[:,i],sigma*z))
    height_rate=spsolve(moving_mass.tocsc(),right)
    mesh_velocity=body_motion+np.column_stack([jacobian@height_rate for jacobian in height_jacobians])
    graph_rate=np.einsum('qc,fc->fq',shape,height_rate[local])
    horizontal_motion=np.einsum('qc,fci->fqi',shape,mesh_velocity[faces,:2])
    # In the exterior sponge the surface rate includes -sigma*eta. Add it
    # back when recovering the physical fluid velocity for Bernoulli.
    eulerian_rate=graph_rate-np.einsum('fqi,fi->fq',horizontal_motion,slope)+sigma*z
    vertical_velocity=(eulerian_rate+
        np.einsum('fi,fi->f',slope,potential_gradient)[:,None])/(1+np.sum(slope*slope,axis=1))[:,None]
    bernoulli=-9.81*z-.5*np.sum(potential_gradient**2,axis=1)[:,None]
    bernoulli=bernoulli+.5*(1+np.sum(slope*slope,axis=1))[:,None]*vertical_velocity**2
    bernoulli+=np.einsum('fqi,fi->fq',horizontal_motion,potential_gradient)
    bernoulli-=sigma*np.einsum('qc,fc->fq',shape,psi[local])
    load=np.zeros(len(free_nodes))
    for corner in range(3):
        np.add.at(load,local[:,corner],
            np.einsum('f,q,q,fq->f',area,weights,shape[:,corner],bernoulli))
    potential_rate=spsolve(mass.tocsc(),load)
    if not np.all(np.isfinite(height_rate)) or not np.all(np.isfinite(potential_rate)):
        raise FloatingPointError('nonfinite moving free-surface rates')
    return height_rate,potential_rate,mesh_velocity
