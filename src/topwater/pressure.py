"""Acceleration-potential/body coupling on one explicit P0 mixed boundary."""
from dataclasses import dataclass
import numpy as np
from .frames import rotation_matrix
from .bem import assemble_laplace
from .boundary import BoundarySystem
from .rankine import TrianglePanels
from .fluid_boundary import FluidBoundary


@dataclass(frozen=True,init=False)
class PressureSystem:
    """Keep one immutable geometry, boundary assignment and factors together."""
    panels: object
    body_mask: np.ndarray
    system: BoundarySystem
    boundary: FluidBoundary

    def __init__(self,boundary,dirichlet,*,length_scale,controls=None):
        if not isinstance(boundary,FluidBoundary):raise ValueError('owned closed fluid/body geometry required')
        panels=TrianglePanels.from_mesh(boundary.mesh)
        body=boundary.tags=='body';mask=np.asarray(dirichlet)
        if (body.dtype!=np.bool_ or body.shape!=(len(panels.areas),) or not np.any(body)
                or mask.dtype!=np.bool_ or mask.shape!=body.shape or np.any(mask[body])):
            raise ValueError('one immutable Neumann body assignment required')
        body=body.copy();body.setflags(write=False)
        operators=assemble_laplace(panels,controls=controls)
        system=BoundarySystem(operators,dirichlet=mask,areas=panels.areas,
            length_scale=length_scale,bounded=True)
        object.__setattr__(self,'panels',panels);object.__setattr__(self,'body_mask',body)
        object.__setattr__(self,'system',system)
        object.__setattr__(self,'boundary',boundary)


def body_normal_acceleration(centers,normals,state,velocity,hessian):
    """Known Eulerian chi_n, excluding six unknown CG accelerations.

    Fluid-outward n; fluid gradient/Hessian are mandatory physical inputs.
    This owner does not recover derivatives or invent them when unavailable.
    """
    x=np.asarray(centers);n=np.asarray(normals);u=np.asarray(velocity);H=np.asarray(hessian)
    if (x.ndim!=2 or x.shape[1]!=3 or n.shape!=x.shape or u.shape!=x.shape
            or H.shape!=(len(x),3,3) or not all(np.isfinite(a).all() for a in (x,n,u,H))):
        raise ValueError('finite body centers, normals, fluid velocities and Hessians required')
    omega=rotation_matrix(state.orientation)@state.angular_velocity
    r=x-state.position;v=state.linear_velocity+np.cross(omega,r)
    return (np.einsum('ni,ni->n',n,np.cross(omega,np.cross(omega,r)))
        -np.einsum('ni,nij,nj->n',n,H,v)+np.einsum('ni,ni->n',v-u,np.cross(omega,n)))


def hydrostatic_load(triangles,normals,center,density,gravity):
    """Exact linear pressure force and quadratic moment on planar wet facets."""
    t=np.asarray(triangles);n=np.asarray(normals)
    area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    z=t[:,:,2];r=t-np.asarray(center)
    mean_rz=(r.sum(axis=1)*z.sum(axis=1)[:,None]+np.einsum('nki,nk->ni',r,z))/12
    force=-density*gravity*np.sum(area[:,None]*z.mean(axis=1)[:,None]*n,axis=0)
    moment=-density*gravity*np.sum(area[:,None]*np.cross(mean_rz,n),axis=0)
    return np.r_[force,moment]


@dataclass(frozen=True)
class PressureAcceleration:
    acceleration_world: np.ndarray
    added_mass: np.ndarray
    acceleration_potential: np.ndarray
    fluid_load: np.ndarray
    external_load: np.ndarray
    diagnostics: dict


def pressure_acceleration(owner,*,known_potential,
        known_normal_acceleration,body_speed_squared,density,gravity,applied_load=None):
    """Solve seven chi fields and six coupled rigid accelerations.

    Caller supplies chi Dirichlet data, known Neumann terms (including explicit
    outer/bottom conditions) and body |grad(phi)|^2. No derivative data is inferred.
    A matching factorization and fluid-outward panels are required. Weight and
    the world-frame gyroscopic load are included once, here.
    """
    if not isinstance(owner,PressureSystem):raise ValueError('one owned pressure geometry/factorization required')
    panels,system,mask=owner.panels,owner.system,owner.body_mask;n=len(panels.areas)
    body,state=owner.boundary.body,owner.boundary.state
    d=np.asarray(known_potential);q=np.asarray(known_normal_acceleration);speed=np.asarray(body_speed_squared)
    if (mask.shape!=(n,) or mask.dtype!=np.bool_ or not np.any(mask)
            or d.shape!=(n,) or q.shape!=(n,) or speed.shape!=(int(mask.sum()),)
            or not all(np.isfinite(a).all() for a in (d,q,speed)) or np.any(speed<0)
            or not np.isfinite([density,gravity]).all() or min(density,gravity)<=0
            or np.any(system.dirichlet[mask]) or not np.array_equal(system.areas,panels.areas)):
        raise ValueError('matching Neumann body panels and explicit finite pressure/kinetic data required')
    applied=np.zeros(6) if applied_load is None else np.asarray(applied_load,float)
    if applied.shape!=(6,) or not np.isfinite(applied).all():raise ValueError('finite six-component world-CG load required')
    centers=panels.centroids[mask];normal=panels.normals[mask];areas=panels.areas[mask]
    modes=np.c_[normal,np.cross(centers-state.position,normal)]
    phi=np.zeros((n,7));neumann=np.zeros((n,7));phi[:,0]=d;neumann[:,0]=q;neumann[mask,1:]=modes
    result=system.solve(potential=phi,normal_velocity=neumann)
    load=modes.T*areas
    added=density*load@result.potential[mask,1:]
    known=-density*load@(result.potential[mask,0]+.5*speed)
    known+=hydrostatic_load(panels.triangles[mask],normal,state.position,density,gravity)
    R=rotation_matrix(state.orientation);I=R@body.inertia_cg@R.T;omega=R@state.angular_velocity
    mass=np.zeros((6,6));mass[:3,:3]=body.mass*np.eye(3);mass[3:,3:]=I
    external=applied.copy();external[2]-=body.mass*gravity;external[3:]-=np.cross(omega,I@omega)
    total=mass+added
    if not np.isfinite(total).all() or np.linalg.cond(total)>1e12:raise np.linalg.LinAlgError('singular coupled mass matrix')
    acceleration=np.linalg.solve(total,known+external)
    chi=result.potential[:,0]+result.potential[:,1:]@acceleration
    fluid=-density*load@(chi[mask]+.5*speed)+hydrostatic_load(panels.triangles[mask],normal,state.position,density,gravity)
    balance=mass@acceleration-fluid-external
    scale=1+np.linalg.norm(mass@acceleration)+np.linalg.norm(fluid)+np.linalg.norm(external)
    residual=np.linalg.norm(balance)/scale
    if not np.isfinite(residual) or residual>1e-10:raise np.linalg.LinAlgError('pressure/body balance failed')
    return PressureAcceleration(acceleration,added,chi,fluid,external,
        dict(balance_relative=float(residual),raw_added_mass_reciprocity=float(np.max(abs(added-added.T))),
             boundary=result.diagnostics,scope='pressure block only; full derivative recovery/RHS not supplied'))
