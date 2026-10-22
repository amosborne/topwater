"""Current-geometry nonlinear six-DOF canoe experiment, under construction."""
from dataclasses import dataclass
import json
from types import SimpleNamespace
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.spatial.transform import Rotation
from topwater import SimulationInputs,RigidBody,InitialState
from topwater.hydrostatics import equilibrate_heave
from topwater.rankine import TrianglePanels
from tests.reference_hull import input_data
from .canoe import body_load
from .canoe_contact import CanoeContact
from .pocket import grid,boundary,stiffness
from .ale_surface import rates


@dataclass(frozen=True)
class CanoeDomain:
    body: RigidBody
    rest: InitialState
    contact: CanoeContact
    reference_points: np.ndarray
    exterior_tets: np.ndarray
    local_mapping: np.ndarray
    free_nodes: np.ndarray
    free_faces: np.ndarray
    local_free_indices: np.ndarray
    column_free_indices: np.ndarray
    volume_nodes: int
    column_index: np.ndarray
    wall_mapping: np.ndarray
    gpu: bool


def make_domain(n=12,gpu=True):
    data=SimulationInputs.from_json(json.dumps(input_data(2)))
    rest,_=equilibrate_heave(data.body,data.initial_state)
    points,tets,near,top_nodes=grid(n,0.,0)
    near &= np.all(points[tets][:,:,2]>=-1.,axis=1)
    local_boundary=boundary(points,tets[near])
    top=np.all(points[local_boundary,2]==0,axis=1)
    contact=CanoeContact.make(data.body,rest,points[local_boundary[~top]],
        points[np.unique(local_boundary[top])])
    lookup={tuple(point):i for i,point in enumerate(points)}
    additional=[];mapping=[]
    for point in contact.reference_points:
        key=tuple(point)
        if key not in lookup:
            lookup[key]=len(points)+len(additional)
            additional.append(point)
        mapping.append(lookup[key])
    mapping=np.asarray(mapping)
    exterior_boundary=boundary(points,tets[~near])
    exterior_top=exterior_boundary[np.all(points[exterior_boundary,2]==0,axis=1)]
    free_faces=np.vstack((exterior_top,mapping[contact.faces[contact.free_faces]]))
    free_nodes=np.unique(free_faces)
    free_lookup=np.full(len(points)+len(additional),-1,dtype=int)
    free_lookup[free_nodes]=np.arange(len(free_nodes))
    return CanoeDomain(data.body,rest,contact,np.vstack((points,additional)),tets[~near],
        mapping,free_nodes,free_faces,free_lookup[mapping[contact.free_nodes]],
        free_lookup[top_nodes],len(points),np.repeat(np.arange(len(top_nodes)),n//3+1),
        mapping[contact.wall_nodes],gpu)


def geometry(domain,position,orientation,eta):
    local=domain.contact.evaluate(position,orientation,eta[domain.local_free_indices])
    points=domain.reference_points.copy()
    # Interior deformation ends at the pocket bottom. Columns underneath the
    # canoe have no free node and no exterior elements above that level.
    column_eta=np.zeros(len(domain.column_free_indices))
    active=domain.column_free_indices>=0
    column_eta[active]=eta[domain.column_free_indices[active]]
    volume=points[:domain.volume_nodes]
    volume[:,2]+=np.maximum(volume[:,2]+1,0)*column_eta[domain.column_index]
    wall_points=points[domain.wall_mapping].copy()
    points[domain.local_mapping]=local.points
    points[domain.wall_mapping]=wall_points
    eta_jacobian=np.zeros_like(points);eta_jacobian[domain.free_nodes,2]=1
    eta_jacobian[domain.local_mapping]=local.eta_jacobian
    body_jacobian=np.zeros((len(points),3,6))
    body_jacobian[domain.local_mapping]=local.body_jacobian
    return points,eta_jacobian,body_jacobian


def field(domain,position,orientation,eta,psi,velocity):
    points,eta_jacobian,body_jacobian=geometry(domain,position,orientation,eta)
    contact=domain.contact;mapping=domain.local_mapping
    panels=TrianglePanels(points[mapping][contact.faces],1)
    if domain.gpu:
        from .cuda_operator import assemble_energy
    else:
        from .operator import assemble_energy
    local,_=assemble_energy(contact.faces,panels,3)
    energy=stiffness(points,domain.exterior_tets)
    energy+=coo_matrix((local.ravel(),(np.repeat(mapping,len(mapping)),
        np.tile(mapping,len(mapping)))),shape=energy.shape).tocsr()
    load=np.zeros((len(points),6))
    # body_load integrates fluid-outward normals, including waterline rows.
    local_mesh=SimpleNamespace(vertices=points[mapping],faces=contact.faces)
    load[mapping]=body_load(local_mesh,panels,contact.body_faces,position)
    active=np.flatnonzero(energy.diagonal()>0)
    unknown=np.setdiff1d(active,domain.free_nodes)
    factor=splu(energy[unknown][:,unknown].tocsc())
    potential=np.zeros(len(points));potential[domain.free_nodes]=psi
    right=load@velocity-energy@potential
    potential[unknown]=factor.solve(right[unknown])
    flux=(energy@potential-load@velocity)[domain.free_nodes]
    acceleration_potentials=np.zeros((len(points),6))
    acceleration_potentials[unknown]=factor.solve(load[unknown])
    coordinates=points[domain.exterior_tets]
    tet_volume=np.abs(np.linalg.det(coordinates[:,1:]-coordinates[:,:1]))/6
    local_volume=np.sum(panels.areas*np.einsum('fi,fi->f',panels.centroids,panels.normals))/3
    z=panels.triangles[:,:,2]
    local_first_moment=np.sum(panels.areas*panels.normals[:,2]*
        (z.sum(axis=1)**2+(z*z).sum(axis=1)))/24
    return dict(points=points,eta_jacobian=eta_jacobian,body_jacobian=body_jacobian,
        potential=potential,flux=flux,acceleration_potentials=acceleration_potentials,
        fluid_volume=float(local_volume+tet_volume.sum()),
        fluid_energy=float(500*potential@energy@potential+9810*(local_first_moment+
            tet_volume@coordinates[:,:,2].mean(axis=1))))


def surface(domain,solution,psi,velocity):
    height_rate,potential_rate,mesh_velocity=rates(solution['points'],domain.free_faces,
        domain.free_nodes,psi,solution['flux'],solution['eta_jacobian'],
        solution['body_jacobian'],velocity)
    column_rate=np.zeros(len(domain.column_free_indices))
    active=domain.column_free_indices>=0
    column_rate[active]=height_rate[domain.column_free_indices[active]]
    wall=domain.wall_mapping
    volume_wall=wall[wall<domain.volume_nodes]
    mesh_velocity[volume_wall,2]=np.maximum(domain.reference_points[volume_wall,2]+1,0)*\
        column_rate[domain.column_index[volume_wall]]
    return height_rate,potential_rate,mesh_velocity


def moved_orientation(orientation,omega,dt):
    return (Rotation.from_rotvec(dt*omega)*
        Rotation.from_quat(orientation,scalar_first=True)).as_quat(scalar_first=True)


def rhs(domain,state,derivative_step=1e-4,*,diagnostics=None):
    """Fresh nonlinear fluid/pressure/body stage; world twist is held in FD.

    Two additional current-geometry field solves differentiate the known ALE
    potential in the actual state-rate direction. The six acceleration columns
    supply its remaining linear dependence on body acceleration. This avoids
    estimating a fluid Hessian at the hull; it is not a frozen-operator model.
    """
    from .body_pressure import evaluate_body_pressure
    size=len(domain.free_nodes)
    position=state[:3]
    orientation=state[3:7]/np.linalg.norm(state[3:7])
    velocity=state[7:13]
    eta=state[13:13+size];psi=state[13+size:]
    current=field(domain,position,orientation,eta,psi,velocity)
    height_rate,potential_rate,mesh_velocity=surface(domain,current,psi,velocity)
    offsets=[]
    for dt in (derivative_step,-derivative_step):
        moved=field(domain,position+dt*velocity[:3],
            moved_orientation(orientation,velocity[3:],dt),
            eta+dt*height_rate,psi+dt*potential_rate,velocity)
        offsets.append(moved['potential'])
    known_potential_rate=(offsets[0]-offsets[1])/(2*derivative_step)
    mapping=domain.local_mapping
    rotation=Rotation.from_quat(orientation,scalar_first=True).as_matrix()
    pressure=evaluate_body_pressure(current['points'][mapping],domain.contact.faces,
        domain.contact.body_faces,position,domain.body.mass,
        rotation@domain.body.inertia_cg@rotation.T,velocity,
        current['potential'][mapping],known_potential_rate[mapping],mesh_velocity[mapping],
        current['acceleration_potentials'][mapping])
    quat_rate=.5*np.r_[-velocity[3:]@orientation[1:],
        orientation[0]*velocity[3:]+np.cross(velocity[3:],orientation[1:])]
    result=np.r_[velocity[:3],quat_rate,pressure.acceleration,height_rate,potential_rate]
    if not np.all(np.isfinite(result)):
        raise FloatingPointError('nonfinite nonlinear canoe stage')
    if diagnostics is not None:
        diagnostics.update(pressure=pressure.diagnostics,
            acceleration=pressure.acceleration.tolist(),
            hydrodynamic_load=pressure.hydrodynamic_load.tolist(),
            added_mass=pressure.added_mass.tolist(),
            fluid_volume=current['fluid_volume'],
            mechanical_energy=current['fluid_energy']+domain.body.mass*9.81*position[2]+
                .5*domain.body.mass*velocity[:3]@velocity[:3]+
                .5*velocity[3:]@(rotation@domain.body.inertia_cg@rotation.T)@velocity[3:],
            maximum_mesh_speed=float(np.linalg.norm(mesh_velocity,axis=1).max()))
    return result


def initial_state(domain,perturbed=True):
    size=len(domain.free_nodes)
    position=domain.rest.position.copy()
    orientation=domain.rest.orientation.copy()
    velocity=np.zeros(6);eta=np.zeros(size);psi=np.zeros(size)
    if perturbed:
        position+=np.array([.0001,-.0001,.0002])
        orientation=moved_orientation(orientation,np.array([.001,-.0007,.0004]),1.)
        velocity=np.array([.001,-.001,.002,.002,-.001,.001])
        x,y=domain.reference_points[domain.free_nodes,:2].T
        eta=.0001*np.cos(.5*x)*np.cos(y)
        psi=.0001*np.cos(.5*x)
    return np.r_[position,orientation,velocity,eta,psi]


def run(n=12,duration=.02,rtol=1e-4,gpu=True,derivative_step=1e-4,max_step=.005):
    from time import perf_counter
    from scipy.integrate import RK23
    if not np.isfinite(duration) or duration<=0 or not 0<rtol<1 or derivative_step<=0 or max_step<=0:
        raise ValueError('positive duration/derivative step and 0 < rtol < 1 required')
    domain=make_domain(n,gpu)
    state=initial_state(domain)
    calls=0;start=perf_counter();last={}
    def evaluate(time,value):
        nonlocal calls
        output=rhs(domain,value,derivative_step,diagnostics=last)
        calls+=1
        if calls%10==0:
            print(json.dumps(dict(progress=True,time=time,calls=calls)),flush=True)
        return output
    solver=RK23(evaluate,0.,state,duration,rtol=rtol,atol=1e-8,max_step=max_step)
    initial={};rhs(domain,state,derivative_step,diagnostics=initial)
    states=[state];times=[0.];energies=[initial['mechanical_energy']];volumes=[initial['fluid_volume']]
    while solver.status=='running':
        message=solver.step()
        if solver.status=='failed':
            raise RuntimeError(message)
        states.append(solver.y.copy());times.append(solver.t)
        energies.append(last['mechanical_energy']);volumes.append(last['fluid_volume'])
    states=np.asarray(states)
    diagnostics={};rhs(domain,states[-1],derivative_step,diagnostics=diagnostics)
    report=dict(qualified=False,case='nonlinear six-DOF canoe with moving contact',
        simulated_seconds=solver.t,n=n,gpu=gpu,rtol=rtol,derivative_step=derivative_step,max_step=max_step,
        accepted_steps=len(times)-1,rhs_evaluations=calls,
        maximum_position_change=np.max(np.abs(states[:,:3]-states[0,:3]),axis=0).tolist(),
        maximum_velocity_change=np.max(np.abs(states[:,7:13]-states[0,7:13]),axis=0).tolist(),
        quaternion_norm_error=float(np.max(np.abs(np.linalg.norm(states[:,3:7],axis=1)-1))),
        maximum_fluid_volume_change=float(np.max(np.abs(np.asarray(volumes)-volumes[0]))),
        maximum_mechanical_energy_change=float(np.max(np.abs(np.asarray(energies)-energies[0]))),
        final=diagnostics,
        scope='Current nonlinear geometry, free surface, pressure and six-DOF acceleration. '
            'Fixed contact topology and reflecting finite tank; no remeshing, outgoing '
            'waves or finite-frequency Capytaine qualification yet. '
            'No performance claim before those correctness requirements pass.',
        diagnostic_wall_seconds=perf_counter()-start)
    return report,dict(times=np.asarray(times),states=states,energy=energies,fluid_volume=volumes)


if __name__=='__main__':
    import argparse
    from pathlib import Path
    from threadpoolctl import threadpool_limits
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=12)
    parser.add_argument('--duration',type=float,default=.02)
    parser.add_argument('--rtol',type=float,default=1e-4)
    parser.add_argument('--derivative-step',type=float,default=1e-4)
    parser.add_argument('--max-step',type=float,default=.005)
    parser.add_argument('--cpu',action='store_true')
    args=parser.parse_args()
    with threadpool_limits(1):
        report,states=run(args.n,args.duration,args.rtol,not args.cpu,args.derivative_step,args.max_step)
    path=Path('artifacts');path.mkdir(exist_ok=True)
    suffix=f'n{args.n}-duration{args.duration:g}-rtol{args.rtol:g}-step{args.max_step:g}'
    np.savez(path/f'nonlinear-canoe-{suffix}.npz',**states)
    (path/f'nonlinear-canoe-{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
