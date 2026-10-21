"""Independent directional energy variation of actual moving-graph BEM."""
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from .evolution import make_tank,rhs


def run(resolution):
    data=make_tank(resolution);mesh,fs=data.mesh,data.surface_nodes;mass=data.projected_mass
    x=mesh.vertices[fs,0]+2*np.pi;g=9.81;a=.02;theta=np.pi/4
    omega=np.sqrt(g*np.tanh(4.));eta=a*np.cos(x)*np.cos(theta)
    psi=-a*g/omega*np.cos(x)*np.sin(theta);direction=np.cos(2*x)
    base=rhs(data,eta,psi)
    predicted=-direction@mass@(base[1]+g*eta)
    rows=[]
    for epsilon in (1e-4,5e-5):
        energies=[]
        for sign in (1,-1):
            moved=eta+sign*epsilon*direction
            value=rhs(data,moved,psi)[3]-.5*g*moved@mass@moved
            energies.append(value)
        derivative=(energies[0]-energies[1])/(2*epsilon)
        rows.append(dict(epsilon=epsilon,finite_difference=derivative,dynamic_prediction=float(predicted),
            relative_discrepancy=float(abs(derivative-predicted)/abs(derivative))))
    return dict(resolution=resolution,results=rows)


if __name__=='__main__':
    rows=[]
    with threadpool_limits(1):
        for resolution in ((8,4,4),(16,8,8)):
            row=run(resolution);rows.append(row);print(json.dumps(row),flush=True)
    Path('artifacts').mkdir(exist_ok=True)
    Path('artifacts/free-surface-shape-energy.json').write_text(json.dumps(rows,indent=2)+'\n')
