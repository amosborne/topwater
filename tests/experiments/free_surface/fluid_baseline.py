"""Adaptive nonlinear near/far fluid run. Moving canoe coupling is unfinished."""
from time import perf_counter
_process_start=perf_counter()
import argparse
from dataclasses import replace
import json
from pathlib import Path
import numpy as np
from scipy.integrate import RK23
from threadpoolctl import threadpool_limits
from tests.benchmark_rankine import peak_rss_bytes
from .pocket_evolution import make_tank,rhs


def run(n=24,periods=.25,rtol=1e-3,gpu=True):
    if periods<=0 or not np.isfinite(periods) or not 0<rtol<1:
        raise ValueError('periods must be positive and 0 < rtol < 1')
    start=perf_counter();data=make_tank(n,gpu=gpu)
    surface=data.surface;fs=surface.surface_nodes;mass=surface.projected_mass
    x=surface.mesh.vertices[fs,0]+2*np.pi
    amplitude=.02;k=.25;g=9.81
    omega=np.sqrt(g*k*np.tanh(4*k));period=2*np.pi/omega
    size=len(fs);eta0=amplitude*np.cos(k*x);psi0=np.zeros(size)
    preparation_seconds=perf_counter()-start
    # A nonzero-potential, deformed-surface comparison checks the entire stage,
    # rather than the trivially zero fluid velocity of the initial turning point.
    phase_eta=eta0/np.sqrt(2)
    phase_psi=-amplitude*g/omega*np.cos(k*x)/np.sqrt(2)
    stage=perf_counter();verification={}
    if gpu:
        cpu=rhs(replace(data,gpu=False),phase_eta,phase_psi)
        candidate=rhs(data,phase_eta,phase_psi)
        for name,index in (('elevation_rate',0),('potential_rate',1)):
            difference=candidate[index]-cpu[index]
            relative=float(np.sqrt(difference@mass@difference/(cpu[index]@mass@cpu[index])))
            verification[name+'_relative_difference']=relative
        if max(verification.values())>1e-10:
            raise RuntimeError('GPU and CPU complete fluid stages disagree')
        gauge=rhs(data,phase_eta,np.ones(size))
        verification['gauge_elevation_rate_max']=float(np.max(np.abs(gauge[0])))
        verification['gauge_nonlinear_potential_rate_max']=float(np.max(np.abs(gauge[1]+g*phase_eta)))
        verification['deformed_stage_volume_rate']=candidate[2]
        if max(abs(verification[key]) for key in (
                'gauge_elevation_rate_max','gauge_nonlinear_potential_rate_max',
                'deformed_stage_volume_rate'))>1e-10:
            raise RuntimeError('GPU fluid stage fails gauge or volume invariance')
    verification_seconds=perf_counter()-stage
    timings=[]
    def evaluate(time,state):
        timing={}
        value=rhs(data,state[:size],state[size:],timing=timing)
        timings.append(timing)
        if len(timings)%20==0:
            print(json.dumps(dict(progress=True,simulated_time=time,fluid_solves=len(timings),
                integration_wall_seconds=perf_counter()-integration_start)),flush=True)
        return np.r_[value[0],value[1]]
    integration_start=perf_counter()
    solver=RK23(evaluate,0.,np.r_[eta0,psi0],periods*period,
        rtol=rtol,atol=amplitude*rtol*.01,max_step=period/32)
    states=[solver.y.copy()];times=[solver.t];attempts=0
    initial_energy=.5*g*eta0@mass@eta0
    energies=[initial_energy]
    while solver.status=='running':
        before=solver.nfev
        message=solver.step()
        if solver.status=='failed':
            raise RuntimeError(message)
        # Each RK23 attempt uses its three new stages; FSAL reuses the start.
        evaluations=solver.nfev-before
        if evaluations%solver.n_stages:
            raise RuntimeError('unexpected RK23 stage accounting')
        attempts+=evaluations//solver.n_stages
        states.append(solver.y.copy());times.append(solver.t)
        # RK23's last derivative evaluation is at the accepted endpoint (FSAL).
        energies.append(timings[-1]['total_energy'])
    final_timing={}
    final=rhs(data,solver.y[:size],solver.y[size:],timing=final_timing)
    timings.append(final_timing)
    integration_seconds=perf_counter()-integration_start
    states=np.asarray(states);times=np.asarray(times)
    volume=states[:,:size]@(np.ones(size)@mass)
    report=dict(qualified=False,phase='body-free closed-tank fluid baseline',n=n,gpu=gpu,
        amplitude=amplitude,wavenumber=k,rtol=rtol,simulated_seconds=float(solver.t),
        preparation_seconds=preparation_seconds,verification_seconds=verification_seconds,
        integration_seconds=integration_seconds,
        integration_wall_per_simulated_second=integration_seconds/solver.t,
        solver_rhs_evaluations=solver.nfev,total_integration_fluid_solves=len(timings),
        accepted_steps=len(times)-1,rejected_steps=attempts-(len(times)-1),
        minimum_step=float(np.diff(times).min()),maximum_step=float(np.diff(times).max()),
        final_energy_relative_change=float(final[3]/initial_energy-1),
        maximum_accepted_energy_relative_change=float(np.max(np.abs(np.asarray(energies)/initial_energy-1))),
        maximum_volume_change=float(np.max(np.abs(volume-volume[0]))),
        verification=verification,stage_timings=timings,peak_rss_bytes=peak_rss_bytes(),
        scope='Every RK stage rebuilds current geometry and the coupled BEM/FEM operator. '
            'Final energy diagnostic is included in integration time/solve count. '
            'No moving hull, pressure/body integration or outgoing-wave boundary yet. '
            'This does not qualify the requested complete canoe simulation.')
    return report,dict(times=times,elevation=states[:,:size],potential=states[:,size:],energy=energies)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=24)
    parser.add_argument('--periods',type=float,default=.25)
    parser.add_argument('--rtol',type=float,default=1e-3)
    parser.add_argument('--cpu',action='store_true')
    args=parser.parse_args()
    with threadpool_limits(1):
        report,states=run(args.n,args.periods,args.rtol,not args.cpu)
    Path('artifacts').mkdir(exist_ok=True)
    backend='cpu' if args.cpu else 'gpu'
    suffix=f'{backend}-n{args.n}-rtol{args.rtol:g}-periods{args.periods:g}'
    np.savez(f'artifacts/fluid-baseline-{suffix}.npz',**states)
    report['process_seconds_through_state_output']=perf_counter()-_process_start
    Path(f'artifacts/fluid-baseline-{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='stage_timings'}),flush=True)
