"""Nonlinear graph rates on the canoe's moving, body-fitted surface mesh."""
import numpy as np
from scipy.linalg import solve
from tests.laplace_reference import triangle_samples


def rates(points,faces,free_nodes,psi,flux,eta_jacobian,body_jacobian,velocity):
    """Return ALE height/potential rates and nodal mesh velocity.

    Each free node owns its height coordinate. Its horizontal mesh motion is
    J_eta*height_rate + J_body*world_twist. Contact nodes slide on material
    hull edges; all other free nodes retain fixed horizontal coordinates.
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
    mass=np.zeros((len(free_nodes),len(free_nodes)))
    transport=np.zeros_like(mass)
    body_motion=np.einsum('nij,j->ni',body_jacobian,velocity)
    contact_slope=np.einsum('fi,fci->fc',slope,eta_jacobian[faces,:2])
    prescribed_transport=np.einsum('fi,fci->fc',slope,body_motion[faces,:2])
    right=flux.copy()
    for i in range(3):
        for j in range(3):
            np.add.at(mass,(local[:,i],local[:,j]),element_mass[:,i,j])
            np.add.at(transport,(local[:,i],local[:,j]),
                element_mass[:,i,j]*contact_slope[:,j])
            np.add.at(right,local[:,i],element_mass[:,i,j]*prescribed_transport[:,j])
    height_rate=solve(mass-transport,right)
    mesh_velocity=body_motion.copy()
    mesh_velocity[free_nodes]+=eta_jacobian[free_nodes]*height_rate[:,None]
    shape,weights=triangle_samples(np.eye(3)[None],3)
    shape=shape[0]
    graph_rate=np.einsum('qc,fc->fq',shape,height_rate[local])
    horizontal_motion=np.einsum('qc,fci->fqi',shape,mesh_velocity[faces,:2])
    eulerian_rate=graph_rate-np.einsum('fqi,fi->fq',horizontal_motion,slope)
    vertical_velocity=(eulerian_rate+
        np.einsum('fi,fi->f',slope,potential_gradient)[:,None])/(1+np.sum(slope*slope,axis=1))[:,None]
    z=np.einsum('qc,fc->fq',shape,triangles[:,:,2])
    bernoulli=-9.81*z-.5*np.sum(potential_gradient**2,axis=1)[:,None]
    bernoulli=bernoulli+.5*(1+np.sum(slope*slope,axis=1))[:,None]*vertical_velocity**2
    bernoulli+=np.einsum('fqi,fi->fq',horizontal_motion,potential_gradient)
    load=np.zeros(len(free_nodes))
    for corner in range(3):
        np.add.at(load,local[:,corner],
            np.einsum('f,q,q,fq->f',area,weights,shape[:,corner],bernoulli))
    potential_rate=solve(mass,load,assume_a='pos')
    if not np.all(np.isfinite(height_rate)) or not np.all(np.isfinite(potential_rate)):
        raise FloatingPointError('nonfinite moving free-surface rates')
    return height_rate,potential_rate,mesh_velocity
