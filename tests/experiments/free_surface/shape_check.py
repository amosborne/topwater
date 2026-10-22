"""Independent directional energy variation of actual moving-graph BEM."""
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from .evolution import make_tank,rhs


def check(data,evaluate,k=1.):
    mesh,fs=data.mesh,data.surface_nodes;mass=data.projected_mass
    x=mesh.vertices[fs,0]+2*np.pi;g=9.81;a=.02;theta=np.pi/4
    omega=np.sqrt(g*k*np.tanh(4*k));eta=a*np.cos(k*x)*np.cos(theta)
    psi=-a*g/omega*np.cos(k*x)*np.sin(theta);direction=np.cos(2*k*x)
    base=evaluate(eta,psi)
    predicted=-direction@mass@(base[1]+g*eta)
    rows=[]
    for epsilon in (1e-4,5e-5):
        energies=[]
        for sign in (1,-1):
            moved=eta+sign*epsilon*direction
            value=evaluate(moved,psi)[3]-.5*g*moved@mass@moved
            energies.append(value)
        derivative=(energies[0]-energies[1])/(2*epsilon)
        rows.append(dict(epsilon=epsilon,finite_difference=derivative,dynamic_prediction=float(predicted),
            relative_discrepancy=float(abs(derivative-predicted)/abs(derivative))))
    return dict(wavenumber=k,results=rows)


def run(resolution):
    data=make_tank(resolution)
    return dict(resolution=resolution,**check(data,lambda eta,psi:rhs(data,eta,psi)))


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pocket',action='store_true')
    args=parser.parse_args()
    rows=[]
    with threadpool_limits(1):
        if args.pocket:
            from .pocket_evolution import make_tank as make_pocket,rhs as pocket_rhs
            for n in (12,24):
                data=make_pocket(n,gpu=False)
                for k in (.25,1.):
                    row=dict(n=n,**check(data.surface,lambda eta,psi:pocket_rhs(data,eta,psi),k))
                    rows.append(row);print(json.dumps(row),flush=True)
        else:
            for resolution in ((8,4,4),(16,8,8)):
                row=run(resolution);rows.append(row);print(json.dumps(row),flush=True)
    Path('artifacts').mkdir(exist_ok=True)
    suffix='pocket' if args.pocket else 'bem'
    Path(f'artifacts/free-surface-shape-energy-{suffix}.json').write_text(json.dumps(rows,indent=2)+'\n')
