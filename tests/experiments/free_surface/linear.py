"""Actual BEM frozen-flat linear wave evolution; not nonlinear qualification."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import solve
from scipy.integrate import solve_ivp
from threadpoolctl import threadpool_limits
from .domain import tank
from .operator import assemble_energy
from topwater.rankine import TrianglePanels
from tests.laplace_reference import triangle_samples


def run(resolution,tolerance):
    mesh,free,node_map=tank(*resolution)
    panels=TrianglePanels.from_mesh(mesh)
    energy,_=assemble_energy(mesh.faces,panels,3)
    fs=node_map.ravel();other=np.setdiff1d(np.arange(len(mesh.vertices)),fs)
    s=energy[np.ix_(fs,fs)]-energy[np.ix_(fs,other)]@solve(energy[np.ix_(other,other)],energy[np.ix_(other,fs)],assume_a='pos')
    lookup=np.full(len(mesh.vertices),-1,dtype=int);lookup[fs]=np.arange(len(fs))
    mass=np.zeros_like(s)
    for face,area in zip(lookup[mesh.faces[free]],panels.areas[free]):
        mass[np.ix_(face,face)]+=area*(np.ones((3,3))+np.eye(3))/12
    # Only this frozen-geometry linear experiment may reuse the full DtN map.
    dtn=solve(mass,s,assume_a='pos');n=len(fs);g=9.81
    amplitude=.001
    eta0=amplitude*np.cos(mesh.vertices[fs,0]+2*np.pi)
    omega=np.sqrt(g*np.tanh(4.));period=2*np.pi/omega
    times=np.linspace(0,5*period,1001)
    def rhs(time,state):return np.r_[dtn@state[n:],-g*state[:n]]
    result=solve_ivp(rhs,(times[0],times[-1]),np.r_[eta0,np.zeros(n)],method='RK23',rtol=tolerance,atol=tolerance*1e-3,t_eval=times)
    assert result.success,result.message
    eta,psi=result.y[:n],result.y[n:]
    expected=eta0[:,None]*np.cos(omega*times)
    nodal_reference_error=np.sqrt(np.einsum('it,ij,jt->',eta-expected,mass,eta-expected)/np.einsum('it,ij,jt->',expected,mass,expected))
    # Independent spatial quadrature compares against the continuous wave,
    # rather than its P1 nodal interpolant (which has substantial coarse error).
    points,weights=triangle_samples(mesh.vertices[mesh.faces[free]],16)
    shape,_=triangle_samples(np.eye(3)[None],16)
    exact_shape=amplitude*np.cos(points[:,:,0]+2*np.pi)
    b=np.zeros(n);c=float(np.einsum('f,q,fq->',panels.areas[free],weights,exact_shape**2))
    for corner in range(3):
        load=np.einsum('f,q,q,fq->f',panels.areas[free],weights,shape[0,:,corner],exact_shape)
        np.add.at(b,lookup[mesh.faces[free,corner]],load)
    cosine=np.cos(omega*times)
    numerical_square=np.sum(eta*(mass@eta),axis=0)
    error_square=numerical_square-2*cosine*(b@eta)+cosine**2*c
    continuous_error=float(np.sqrt(np.trapezoid(error_square,times)/np.trapezoid(cosine**2*c,times)))
    e=.5*(g*np.sum(eta*(mass@eta),axis=0)+np.sum(psi*(s@psi),axis=0))
    volume=np.ones(n)@mass@eta
    return dict(resolution=resolution,rtol=tolerance,rhs_evaluations=result.nfev,periods=5,
        continuous_elevation_time_space_relative_rms=continuous_error,
        p1_reference_elevation_time_space_relative_rms=float(nodal_reference_error),
        maximum_relative_energy_drift=float(np.max(np.abs(e/e[0]-1))),
        maximum_volume_change=float(np.max(np.abs(volume-volume[0]))))


if __name__=='__main__':
    rows=[]
    with threadpool_limits(1):
        for resolution,tolerance in (((8,4,4),1e-7),((16,8,8),1e-7),((16,8,8),1e-9)):
            row=run(resolution,tolerance);rows.append(row);print(json.dumps(row),flush=True)
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-linear.json').write_text(json.dumps(rows,indent=2)+'\n')
