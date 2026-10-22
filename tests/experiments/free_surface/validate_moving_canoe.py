"""Equilibrium, pressure derivative and Capytaine checks of the coupled canoe."""
import argparse
import json
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from .moving_canoe import make_domain,initial_state,rhs,field


def run(n=12,gpu=True,capytaine=False,material_subdivisions=2):
    domain=make_domain(n,gpu,material_subdivisions)
    rest=initial_state(domain,False)
    equilibrium=rhs(domain,rest)
    state=initial_state(domain)
    derivatives=[]
    for step in (2e-4,1e-4,5e-5):
        value=rhs(domain,state,step)
        derivatives.append(dict(step=step,acceleration=value[7:13].tolist()))
    accelerations=np.array([row['acceleration'] for row in derivatives])
    changes=np.linalg.norm(np.diff(accelerations,axis=0),axis=1)
    relative=changes/np.linalg.norm(accelerations[-1])
    report=dict(qualified=False,n=n,gpu=gpu,material_subdivisions=material_subdivisions,
        rest_maximum_rate=float(np.max(np.abs(equilibrium))),
        directional_derivatives=derivatives,
        acceleration_relative_changes=relative.tolist(),
        scope='Checks the actual nonlinear six-DOF stage. Instantaneous Capytaine '
            'reaction is one limit, not finite-frequency or nonlinear validation.')
    if capytaine:
        from tests.capytaine_reference import instantaneous_inertia
        current=field(domain,rest[:3],rest[3:7],rest[13:13+len(domain.free_nodes)],
            rest[13+len(domain.free_nodes):],rest[7:13])
        contact=domain.contact;mapping=domain.local_mapping
        body_points=current['points'][mapping][contact.faces[contact.body_faces]]
        expected,damping=instantaneous_inertia(body_points[:,::-1],rest[:3],domain.body.mass)
        diagnostics={};rhs(domain,rest,diagnostics=diagnostics)
        actual=np.asarray(diagnostics['added_mass'])
        report.update(capytaine_instantaneous=dict(actual=actual.tolist(),
            reference=expected.tolist(),
            diagonal_relative_error=(np.diag(actual-expected)/np.diag(expected)).tolist(),
            reference_damping_maximum=float(np.max(np.abs(damping))),
            scope='Same actual wet hull; finite reflecting tank/depth versus '
                'infinite-depth Capytaine. Spatial and domain errors not separated.'))
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n',type=int,default=12)
    parser.add_argument('--cpu',action='store_true')
    parser.add_argument('--capytaine',action='store_true')
    parser.add_argument('--material-subdivisions',type=int,default=2)
    args=parser.parse_args()
    with threadpool_limits(1):
        report=run(args.n,not args.cpu,args.capytaine,args.material_subdivisions)
    Path('artifacts').mkdir(exist_ok=True)
    Path(f'artifacts/nonlinear-canoe-validation-n{args.n}-material{args.material_subdivisions}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
