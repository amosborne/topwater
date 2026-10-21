"""Short actual nonlinear BEM evolution; no production/open-water claim."""
import json
import argparse
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from .evolution import make_tank,rhs


def run(resolution,tolerance):
    data=make_tank(resolution);mesh,fs=data.mesh,data.surface_nodes;n=len(fs);mass=data.projected_mass
    amplitude=.02;g=9.81;omega=np.sqrt(g*np.tanh(4.));period=2*np.pi/omega
    x=mesh.vertices[fs,0]+2*np.pi;eta0=amplitude*np.cos(x);psi0=np.zeros(n)
    initial_energy=.5*g*eta0@mass@eta0
    count=0;start=perf_counter()
    def evaluate(time,state):
        nonlocal count
        result=rhs(data,state[:n],state[n:]);count+=1
        if count%20==0:print(json.dumps(dict(progress=True,resolution=resolution,calls=count,time=time,seconds=perf_counter()-start)),flush=True)
        return np.r_[result[0],result[1]]
    result=solve_ivp(evaluate,(0,period/4),np.r_[eta0,psi0],method='RK23',rtol=tolerance,atol=amplitude*tolerance*.01,max_step=period/32)
    assert result.success,result.message
    final=rhs(data,result.y[:n,-1],result.y[n:,-1])
    volume=np.ones(n)@mass@result.y[:n]
    return dict(resolution=resolution,rtol=tolerance,amplitude=amplitude,
        elapsed_seconds=perf_counter()-start,rhs_evaluations=result.nfev,accepted_steps=len(result.t)-1,
        final_energy_relative_change=float(final[3]/initial_energy-1),
        maximum_volume_change=float(np.max(np.abs(volume-volume[0]))),
        final_elevation_rms=float(np.sqrt(result.y[:n,-1]@mass@result.y[:n,-1]/mass.sum())))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rtol',type=float,default=1e-3)
    arguments=parser.parse_args()
    rows=[]
    with threadpool_limits(1):
        for resolution in ((8,4,4),(16,8,8)):
            row=run(resolution,arguments.rtol);rows.append(row);print(json.dumps(row),flush=True)
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-evolution.json').write_text(json.dumps(rows,indent=2)+'\n')
