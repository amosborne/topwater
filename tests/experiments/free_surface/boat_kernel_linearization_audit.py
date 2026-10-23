"""Bounded one-anchor audit on archived actual canonical boat states.

No trajectory, tuning or target-derived slopes. A chart failure is reported,
not refreshed. Input is a driver NPZ and its geometry-0 JSON/NPZ pair.
"""
import argparse
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from topwater import InitialState, sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import WATER_DENSITY
from .cuda_operator import cp, prepare_energy
from .cuda_kernel_linearization import RawKernelAnchor
from .cuda_shape_gradient import vertex_gradient
from .boat_chart_tangent import chart_tangent
from .boat_action_tangent import refine_tangent, action_terms
from .boat_world_field import current_field
from .incident_baseline_benchmark import coefficients


def run(states,geometry,output):
    states,geometry,output=map(Path,(states,geometry,output))
    if output.exists():
        raise FileExistsError(output)
    reference=dict(np.load(geometry.with_suffix('.npz'),allow_pickle=False))
    metadata=json.loads(geometry.read_text())
    saved=np.load(states,allow_pickle=False)
    times=saved['accepted_times'];sequence=saved['accepted_states']
    values=coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    report=dict(complete=False,qualified=False,scope='One raw-kernel anchor, archived actual states on one guarded chart; no motion qualification',
                states=str(states),geometry=str(geometry),quadrature_order=3,body_refinement=1,samples=[])

    def phase(name,fn,timings):
        cp.cuda.get_current_stream().synchronize();start=perf_counter()
        result=fn()
        cp.cuda.get_current_stream().synchronize();timings[name]=perf_counter()-start
        return result

    try:
        cp.zeros(1).get()
        cache=None
        for time,state in zip(times[:3],sequence[:3],strict=True):
            timings={};sample=dict(time=float(time),timings=timings)
            report['samples'].append(sample)
            pose=InitialState(state[:3],state[3:7]);n=len(reference['parent_nodes'])
            raw,tangent=phase('guarded_chart_tangent',lambda:chart_tangent(reference,pose=pose,
                parent_eta=state[13:13+n],parent_psi=state[13+n:]),timings)
            arrays,dx,ds=phase('refine_tangent',lambda:refine_tangent(raw,tangent['point_jacobian'],
                tangent['parent_map_jacobian'],1),timings)
            points,faces=arrays['points'],arrays['faces']
            sample.update(nodes=len(points),faces=len(faces))
            token=str(reference['node_chart_key'])
            if cache is None:
                cache=phase('anchor_build',lambda:RawKernelAnchor(points,faces,chart_token=token),timings)
                report['retained_pair_bytes']=cache.pair_bytes
            fresh=phase('fresh_energy',lambda:prepare_energy(faces,TrianglePanels(points[faces],1)),timings)
            approximate=phase('anchor_evaluate',lambda:cache.evaluate(points,faces,chart_token=token),timings)
            twist=state[7:13]
            actual=phase('fresh_current_field',lambda:current_field(arrays,fresh.energy,values,twist,time=float(time)),timings)
            candidate=phase('linearized_current_field',lambda:current_field(arrays,approximate.energy,values,twist,time=float(time)),timings)
            exact_g=phase('fresh_gradient',lambda:vertex_gradient(points,faces,actual['potential'],assembly=fresh,density=WATER_DENSITY),timings)
            approx_g=phase('linearized_gradient',lambda:approximate.vertex_gradient(candidate['potential'],density=WATER_DENSITY),timings)
            exact_terms=phase('fresh_action_terms',lambda:action_terms(arrays,actual,exact_g,values,twist,float(time),dx,ds),timings)
            approx_terms=phase('linearized_action_terms',lambda:action_terms(arrays,candidate,approx_g,values,twist,float(time),dx,ds),timings)
            sample.update(max_displacement=float(np.max(abs(points-cache.points))),
                relative_energy_error=float(np.linalg.norm(approximate.energy-fresh.energy)/np.linalg.norm(fresh.energy)),
                relative_potential_error=float(np.linalg.norm(candidate['potential']-actual['potential'])/np.linalg.norm(actual['potential'])),
                relative_gradient_error=float(np.linalg.norm(approx_g-exact_g)/np.linalg.norm(exact_g)),
                gradient_max_absolute_error=float(np.max(abs(approx_g-exact_g))),
                h_constant_norm_approximate=approximate.diagnostics['h_constant_norm'],
                h_constant_norm_fresh=fresh.diagnostics['h_constant_norm'],
                energy_constant_norm_approximate=approximate.diagnostics['energy_constant_norm'],
                energy_constant_norm_fresh=float(np.max(abs(fresh.energy@np.ones(len(points))))))
            for key in ('eta_rate','psi_rate','body_pose_load'):
                sample[key+'_fresh']=exact_terms[key].tolist()
                sample[key+'_approximate']=approx_terms[key].tolist()
                sample[key+'_max_absolute_error']=float(np.max(abs(approx_terms[key]-exact_terms[key])))
            sample['acceleration_potential_max_error']=float(np.max(abs(candidate['acceleration_potentials']-actual['acceleration_potentials'])))
            print(json.dumps(sample),flush=True)
        report['complete']=True
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error))
        raise
    finally:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--states',required=True);parser.add_argument('--geometry',required=True)
    parser.add_argument('--output',required=True)
    print(json.dumps(run(**vars(parser.parse_args())),indent=2))
