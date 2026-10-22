"""Exact rest linearization of the nonlinear variational canoe action.

All six radiation columns are solved together. This validates the linear
limit against Capytaine; it does not replace nonlinear trajectory validation.
The finite tank and sponge require independent spatial/domain convergence.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from .moving_canoe import make_domain,operator_system
from .surface_damping import matrices
from .canoe_radiation import sponge,MODES


def solve(domain,frequencies,*,sponge_strength=4.,sponge_buffer=.5):
    frequencies=np.asarray(frequencies,dtype=float)
    if frequencies.ndim!=1 or not len(frequencies) or not np.all(np.isfinite(frequencies)) or np.any(frequencies<=0):
        raise ValueError('nonempty positive finite frequency sequence required')
    size=len(domain.free_nodes)
    points,_,_,_,energy,load=operator_system(domain,domain.rest.position,
        domain.rest.orientation,np.zeros(size))
    damping,widths=sponge(points[domain.free_nodes],sponge_strength,sponge_buffer)
    mass,weighted=matrices(points,domain.free_faces,domain.free_nodes,damping)
    mass_factor=splu(mass.tocsc())
    active=np.flatnonzero(energy.diagonal()>0)
    free=domain.free_nodes
    indices=(np.repeat(free,size),np.tile(free,size))
    reports=[];potentials=[];heights=[]
    for omega in frequencies:
        # exp(-i omega t): D eta = Schur psi - C V; D psi = -g M eta.
        # Eliminating eta adds D M^-1 D/g to the full fluid energy matrix.
        derivative=weighted.toarray().astype(complex)-1j*omega*mass.toarray()
        robin=(weighted@mass_factor.solve(weighted.toarray())-omega**2*mass.toarray()
            -2j*omega*weighted.toarray())/9.81
        system=energy.astype(complex)+coo_matrix((robin.ravel(),indices),shape=energy.shape).tocsr()
        factor=splu(system[active][:,active].tocsc())
        potential=np.zeros_like(load,dtype=complex)
        potential[active]=factor.solve(load[active].astype(complex))
        response=load.T@potential
        # Columns currently correspond to UNIT VELOCITY; convert heights too.
        trace=potential[free]
        eta=-(mass_factor.solve((derivative@trace).real)
              +1j*mass_factor.solve((derivative@trace).imag))/9.81
        residual=system@potential-load
        dissipation=(2*omega/9.81)*(trace.conj().T@weighted@trace)
        reports.append(dict(omega=float(omega),added_mass=(1000*response.real).tolist(),
            radiation_damping=(1000*omega*response.imag).tolist(),
            relative_system_residual=float(np.linalg.norm(residual[active])/np.linalg.norm(load[active])),
            reciprocity_relative=float(np.linalg.norm(response-response.T)/np.linalg.norm(response)),
            passivity_identity_relative=float(np.linalg.norm(response.imag-dissipation)/max(np.linalg.norm(response.imag),1e-30)),
            minimum_damping_eigenvalue=float(np.linalg.eigvalsh(500*omega*(response.imag+response.imag.T)).min())))
        potentials.append(potential);heights.append(eta)
    report=dict(qualified=False,case='rest linearization of nonlinear variational canoe',
        mode_order=list(MODES),sponge_strength=sponge_strength,sponge_buffer=sponge_buffer,
        sponge_widths=widths.tolist(),water_depth=4.,cases=reports,
        scope='Full six-mode linear limit only. Finite absorbing tank; spatial, domain, '
            'sponge and Capytaine reference convergence required. No nonlinear evolution or timing qualification.')
    return report,dict(potential_per_velocity=np.asarray(potentials),
        height_per_velocity=np.asarray(heights),damping=damping)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=12)
    parser.add_argument('--wet-subdivisions',type=int,default=2)
    parser.add_argument('--exterior-scale',type=float,default=1.)
    parser.add_argument('--omega',type=float,nargs='+',default=[4.])
    parser.add_argument('--sponge-strength',type=float,default=4.)
    parser.add_argument('--sponge-buffer',type=float,default=.5)
    args=parser.parse_args()
    with threadpool_limits(1):
        domain=make_domain(args.n,True,args.wet_subdivisions,args.exterior_scale)
        report,saved=solve(domain,args.omega,sponge_strength=args.sponge_strength,sponge_buffer=args.sponge_buffer)
    report.update(n=args.n,wet_subdivisions=args.wet_subdivisions,exterior_scale=args.exterior_scale)
    suffix=f'n{args.n}-wet{args.wet_subdivisions}-extent{args.exterior_scale:g}-sponge{args.sponge_strength:g}'
    Path('artifacts').mkdir(exist_ok=True)
    np.savez(f'artifacts/linear-canoe-radiation-{suffix}.npz',**saved)
    Path(f'artifacts/linear-canoe-radiation-{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
