"""Exact rest linearization of the nonlinear variational canoe action.

All six radiation columns are solved together. This validates the linear
limit against Capytaine; it does not replace nonlinear trajectory validation.
The finite tank and sponge require independent spatial/domain convergence.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from scipy.sparse import bmat,coo_matrix
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
    active=np.flatnonzero(energy.diagonal()>0)
    free=domain.free_nodes
    free_active=np.searchsorted(active,free)
    if np.any(free_active>=len(active)) or not np.array_equal(active[free_active],free):
        raise ValueError('free surface references an inactive fluid potential node')
    right=np.vstack((load[active],np.zeros((size,6)))).astype(complex)
    reports=[];potentials=[];heights=[]
    for omega in frequencies:
        # exp(-i omega t): D eta = Schur psi - C V; D psi = -g M eta.
        # Retain eta to avoid densifying the entire far free surface through
        # M^-1. Eliminating this block recovers E + P.T*D*M^-1*D*P/g.
        derivative=(weighted-1j*omega*mass).tocoo()
        coupling=coo_matrix((-derivative.data,(free_active[derivative.row],derivative.col)),
            shape=(len(active),size)).tocsc()
        system=bmat([[energy[active][:,active].astype(complex),coupling],
            [coupling.T,-9.81*mass]],format='csc')
        solved=splu(system).solve(right)
        potential=np.zeros_like(load,dtype=complex)
        potential[active]=solved[:len(active)]
        response=load.T@potential
        # Both potential and height columns correspond to unit body velocity.
        trace=potential[free]
        eta=solved[len(active):]
        residual=system@solved-right
        dissipation=(2*omega/9.81)*(trace.conj().T@weighted@trace)
        reports.append(dict(omega=float(omega),added_mass=(1000*response.real).tolist(),
            radiation_damping=(1000*omega*response.imag).tolist(),
            relative_system_residual=float(np.linalg.norm(residual)/np.linalg.norm(right)),
            reciprocity_relative=float(np.linalg.norm(response-response.T)/np.linalg.norm(response)),
            passivity_identity_relative=float(np.linalg.norm(response.imag-dissipation)/max(np.linalg.norm(response.imag),1e-30)),
            minimum_damping_eigenvalue=float(np.linalg.eigvalsh(500*omega*(response.imag+response.imag.T)).min())))
        potentials.append(potential);heights.append(eta)
    wet=points[domain.local_mapping][domain.contact.faces[domain.contact.body_faces]][:,::-1]
    report=dict(qualified=False,case='rest linearization of nonlinear variational canoe',
        mode_order=list(MODES),sponge_strength=sponge_strength,sponge_buffer=sponge_buffer,
        sponge_widths=widths.tolist(),water_depth=4.,cases=reports,
        wetted_triangles_sha256=sha256(np.asarray(wet,dtype='<f8').tobytes()).hexdigest(),
        wetted_panels=len(wet),rotation_center=domain.rest.position.tolist(),
        scope='Full six-mode linear limit only. Finite absorbing tank; spatial, domain, '
            'sponge and Capytaine reference convergence required. Q=0 is checked against '
            'the undamped nonlinear runner; Q>0 adds a validation absorber. '
            'No nonlinear outgoing-wave evolution or timing qualification.')
    return report,dict(potential_per_velocity=np.asarray(potentials),
        height_per_velocity=np.asarray(heights),damping=damping,wetted_triangles=wet)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=12)
    parser.add_argument('--wet-subdivisions',type=int,default=2)
    parser.add_argument('--exterior-scale',type=float,default=1.)
    parser.add_argument('--omega',type=float,nargs='+',default=[4.])
    parser.add_argument('--sponge-strength',type=float,default=4.)
    parser.add_argument('--sponge-buffer',type=float,default=.5)
    parser.add_argument('--capytaine',action='store_true')
    args=parser.parse_args()
    with threadpool_limits(1):
        domain=make_domain(args.n,True,args.wet_subdivisions,args.exterior_scale)
        report,saved=solve(domain,args.omega,sponge_strength=args.sponge_strength,sponge_buffer=args.sponge_buffer)
    if args.capytaine:
        from tests.capytaine_reference import frequency_response
        from .validate_radiation import coefficient_gates
        reference=frequency_response(saved['wetted_triangles'],domain.rest.position,
            domain.body.mass,args.omega,[],water_depth=4.)
        for i,case in enumerate(report['cases']):
            case['capytaine']=dict(added_mass=reference['added_mass'][i].tolist(),
                radiation_damping=reference['radiation_damping'][i].tolist())
            d=np.array([1.,1.,1.,3.5,3.5,3.5])
            scale=1000*3.5**3*d[:,None]*d[None,:]
            case.update(added_mass_nondimensional=(np.asarray(case['added_mass'])/scale).tolist(),
                reference_added_mass_nondimensional=(reference['added_mass'][i]/scale).tolist(),
                damping_nondimensional=(np.asarray(case['radiation_damping'])/(args.omega[i]*scale)).tolist(),
                reference_damping_nondimensional=(reference['radiation_damping'][i]/(args.omega[i]*scale)).tolist())
        report['capytaine_green_function_settings']=reference['green_function_settings']
        report['coefficient_gates']=coefficient_gates(report['cases'])
    report.update(n=args.n,wet_subdivisions=args.wet_subdivisions,exterior_scale=args.exterior_scale)
    frequencies='-'.join(f'{value:g}' for value in args.omega)
    suffix=(f'n{args.n}-wet{args.wet_subdivisions}-extent{args.exterior_scale:g}'
        f'-sponge{args.sponge_strength:g}-buffer{args.sponge_buffer:g}-omega{frequencies}'
        f'-capytaine{int(args.capytaine)}')
    Path('artifacts').mkdir(exist_ok=True)
    np.savez(f'artifacts/linear-canoe-radiation-{suffix}.npz',**saved)
    Path(f'artifacts/linear-canoe-radiation-{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
