"""Current-geometry variational six-DOF canoe validation runner.

Forces differentiate the same discrete action as the fluid solve. The slow
canonical_canoe finite-difference reference remains independent. No nonlinear
wave/body qualification or production-performance claim follows from this.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import RK23
from threadpoolctl import threadpool_limits

from .canonical_canoe import body_map,encode,decode
from .canonical_geometry import gradient as geometric_gradient
from .moving_canoe import make_domain,initial_state


def rhs(domain,canonical,*,diagnostics=None):
    from .cuda_shape_gradient import vertex_gradient
    size=len(canonical)//2
    q=canonical[:size]
    if np.linalg.norm(q[3:6])>=np.pi:
        raise ValueError('local rotation chart boundary reached; recharting required')
    current={}
    physical,hamiltonian,qdot=decode(domain,canonical,diagnostics=current)
    mapping=domain.local_mapping
    rim=domain.contact.rim_clearance(current['points'][mapping],physical[:3],physical[3:7])
    configuration_gradient=geometric_gradient(domain,q,qdot,current['potential'])
    local=vertex_gradient(current['points'][mapping],domain.contact.faces,
        current['potential'][mapping])
    world_body=np.einsum('nij,ni->j',current['body_jacobian'][mapping],local)
    configuration_gradient[:6]+=body_map(q[3:6]).T@world_body
    configuration_gradient[6:]+=sum(matrix[mapping].T@local[:,axis]
        for axis,matrix in enumerate(current['height_jacobians']))
    result=np.r_[qdot,-configuration_gradient]
    if not np.all(np.isfinite(result)):
        raise FloatingPointError('nonfinite variational canoe stage')
    if diagnostics is not None:
        diagnostics.update(physical_state=physical,hamiltonian=hamiltonian,
            fluid_volume=current['fluid_volume'],minimum_rim_clearance=rim)
    return result


def run(n=12,duration=1.,rtol=1e-4,max_step=.005,wet_subdivisions=2,initial_scale=10.):
    if not (np.isfinite(duration) and duration>0 and np.isfinite(rtol) and 0<rtol<1
            and np.isfinite(max_step) and max_step>0):
        raise ValueError('finite positive duration/max_step and 0 < rtol < 1 required')
    start=perf_counter()
    domain=make_domain(n,True,wet_subdivisions)
    canonical=encode(domain,initial_state(domain,scale=initial_scale))
    calls=0
    cache={}
    initial={}
    def evaluate(time,value):
        nonlocal calls
        details={}
        rate=rhs(domain,value,diagnostics=details)
        calls+=1
        cache.clear()
        cache.update(details,time=time,canonical=value.copy())
        if not initial:
            initial.update(cache)
        if calls%20==0:
            print(json.dumps(dict(progress=True,simulated_time=time,rhs_evaluations=calls)),flush=True)
        return rate
    solver=RK23(evaluate,0.,canonical,duration,rtol=rtol,atol=1e-8,max_step=max_step)
    times=[0.]
    states=[initial['physical_state']]
    canonical_states=[canonical.copy()]
    energies=[initial['hamiltonian']]
    volumes=[initial['fluid_volume']]
    clearances=[initial['minimum_rim_clearance']]
    while solver.status=='running':
        message=solver.step()
        if solver.status=='failed':
            raise RuntimeError(message)
        if cache['time']!=solver.t or not np.array_equal(cache['canonical'],solver.y):
            raise RuntimeError('accepted state has no matching variational diagnostics')
        times.append(solver.t)
        states.append(cache['physical_state'])
        canonical_states.append(solver.y.copy())
        energies.append(cache['hamiltonian'])
        volumes.append(cache['fluid_volume'])
        clearances.append(cache['minimum_rim_clearance'])
    states=np.asarray(states)
    energy_change=float(np.max(np.abs(np.asarray(energies)-energies[0])))
    report=dict(qualified=False,case='fully nonlinear variational six-DOF canoe',
        n=n,wet_subdivisions=wet_subdivisions,initial_scale=initial_scale,
        simulated_seconds=solver.t,rtol=rtol,max_step=max_step,
        accepted_steps=len(times)-1,rhs_evaluations=calls,
        initial_excitation_energy=float(energies[0]),
        maximum_energy_change_joules=energy_change,
        maximum_energy_change_over_excitation=energy_change/energies[0] if energies[0]>0 else None,
        maximum_fluid_volume_change=float(np.max(np.abs(np.asarray(volumes)-volumes[0]))),
        maximum_position_change=np.max(np.abs(states[:,:3]-states[0,:3]),axis=0).tolist(),
        maximum_velocity_change=np.max(np.abs(states[:,7:13]-states[0,7:13]),axis=0).tolist(),
        minimum_rim_clearance=float(min(clearances)),
        diagnostic_wall_seconds=perf_counter()-start,
        scope='Current nonlinear BEM/FEM geometry and exact discrete kinetic-energy gradient; '
            'coupled body/surface canonical evolution with adaptive RK23. Undamped reflecting '
            'tank and fixed physical-facet contact topology. Local rotation chart limited to '
            'less than pi. No remeshing, general failure-event coverage, outgoing-wave, '
            'Capytaine family or performance qualification.')
    return report,dict(times=np.asarray(times),states=states,
        canonical_states=np.asarray(canonical_states),hamiltonian=np.asarray(energies),
        fluid_volume=np.asarray(volumes),minimum_rim_clearance=np.asarray(clearances))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=12)
    parser.add_argument('--duration',type=float,default=1.)
    parser.add_argument('--rtol',type=float,default=1e-4)
    parser.add_argument('--max-step',type=float,default=.005)
    parser.add_argument('--wet-subdivisions',type=int,default=2)
    parser.add_argument('--initial-scale',type=float,default=10.)
    args=vars(parser.parse_args())
    with threadpool_limits(1):
        report,saved=run(**args)
    suffix=(f"n{args['n']}-wet{args['wet_subdivisions']}-initial{args['initial_scale']:g}"
        f"-duration{args['duration']:g}-rtol{args['rtol']:g}-step{args['max_step']:g}")
    Path('artifacts').mkdir(exist_ok=True)
    np.savez(f'artifacts/variational-canoe-{suffix}.npz',**saved)
    Path(f'artifacts/variational-canoe-{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
