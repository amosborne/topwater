"""Bounded canonical q3 contracted-scalar gate against frozen old producers."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import perf_counter
import types
import numpy as np
from topwater import InitialState,sea_states
from topwater.rankine import TrianglePanels
from tests.reference_boat import WATER_DENSITY
from . import cuda_operator
from . import cuda_scalar_gradient
from .boat_chart_tangent import chart_tangent
from .boat_action_tangent import refine_tangent,action_terms
from .boat_world_field import current_field
from .incident_baseline_benchmark import coefficients


def run(states,geometry,frozen_report,output):
    states,geometry,frozen_report,output=map(Path,(states,geometry,frozen_report,output))
    if output.exists():raise FileExistsError(output)
    control=json.loads(frozen_report.read_text())
    archive=frozen_report.with_name(frozen_report.stem+'-source')
    for name,digest in control['source_sha256'].items():
        if sha256((archive/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('frozen producer archive hash changed: '+name)
    def frozen_module(name,file):
        module=types.ModuleType(name);module.__file__=str(file)
        module.__package__='tests.experiments.free_surface';sys.modules[name]=module
        exec(compile(file.read_text(),str(file),'exec'),module.__dict__)
        return module
    # Immutable old files are a control only; production keeps one shared owner.
    old_operator=frozen_module('tests.experiments.free_surface._frozen_scalar_operator',archive/'cuda_operator.py')
    old_shape=frozen_module('tests.experiments.free_surface._frozen_scalar_shape',archive/'cuda_shape_gradient.py')
    saved=np.load(states,allow_pickle=False);state=saved['accepted_states'][0];time=float(saved['accepted_times'][0])
    reference=dict(np.load(geometry.with_suffix('.npz'),allow_pickle=False))
    metadata=json.loads(geometry.read_text())
    n=len(reference['parent_nodes']);twist=state[7:13]
    raw,tangent=chart_tangent(reference,pose=InitialState(state[:3],state[3:7]/np.linalg.norm(state[3:7])),
        parent_eta=state[13:13+n],parent_psi=state[13+n:])
    arrays,dx,ds=refine_tangent(raw,tangent['point_jacobian'],tangent['parent_map_jacobian'],1)
    points,faces=arrays['points'],arrays['faces']
    values=coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    folder=Path(__file__).parent
    names=('cuda_pair.cuh','cuda_scalar_shape.cu','cuda_scalar_gradient.py','cuda_shape_gradient.py',
        'cuda_shape.cu','cuda_dual.cuh','cuda_wide_dual.cuh','cuda_operator.py','cuda_energy.cu',
        'boat_scalar_gradient_audit.py','boat_action_tangent.py','boat_chart_tangent.py','boat_world_field.py',
        'incident_baseline_benchmark.py')
    report=dict(complete=False,qualified=False,scope='One actual canonical 472-face q3 state; frozen pre-refactor source comparison, no physical or speed qualification',
        order=3,nodes=len(points),faces=len(faces),time=time,
        frozen_report=str(frozen_report),frozen_report_sha256=sha256(frozen_report.read_bytes()).hexdigest(),
        frozen_source_sha256=control['source_sha256'],states_sha256=sha256(states.read_bytes()).hexdigest(),
        geometry_sha256=sha256(geometry.with_suffix('.npz').read_bytes()).hexdigest(),
        geometry_json_sha256=sha256(geometry.read_bytes()).hexdigest(),
        source_sha256={name:sha256((folder/name).read_bytes()).hexdigest() for name in names},rows=[])
    cp=cuda_operator.cp
    def measured(fn):
        cp.cuda.get_current_stream().synchronize();start=perf_counter();result=fn()
        cp.cuda.get_current_stream().synchronize();return result,perf_counter()-start
    try:
        cp.zeros(1).get()
        old,old_seconds=measured(lambda:old_operator.prepare_energy(faces,TrianglePanels(points[faces],1),order=3))
        new,new_seconds=measured(lambda:cuda_operator.prepare_energy(faces,TrianglePanels(points[faces],1),order=3))
        report.update(old_energy_seconds=old_seconds,new_energy_seconds=new_seconds,
            energy_relative_difference=float(np.linalg.norm(new.energy-old.energy)/np.linalg.norm(old.energy)),
            energy_max_absolute_difference=float(np.max(abs(new.energy-old.energy))),
            mixed_trace_max_absolute_difference=float(np.max(abs(new.mixed_trace.get()-old.mixed_trace.get()))))
        field=current_field(arrays,old.energy,values,twist,time=time)
        phi=field['potential'];gradients=[]
        for label,width in (('old_raw',3),('contracted',1),('contracted',1),('old_raw',3)):
            function=old_shape.vertex_gradient if label=='old_raw' else cuda_scalar_gradient.vertex_gradient
            assembly=old if label=='old_raw' else new
            gradient,seconds=measured(lambda:function(points,faces,phi,assembly=assembly,
                density=WATER_DENSITY,order=3,batch_width=width))
            if not gradients:baseline=gradient
            direction=np.random.default_rng(71).normal(size=points.shape)
            chain=action_terms(arrays,field,gradient,values,twist,time,dx,ds)
            if not gradients:baseline_chain=chain
            kernel=old_shape._kernels[width] if label=='old_raw' else cuda_scalar_gradient._kernels[width]
            report['rows'].append(dict(label=label,batch_width=width,seconds=seconds,
                relative_gradient_difference=float(np.linalg.norm(gradient-baseline)/np.linalg.norm(baseline)),
                maximum_absolute_gradient_difference=float(np.max(abs(gradient-baseline))),
                random_direction_contraction_difference=float(np.sum((gradient-baseline)*direction)),
                parent_rate_max_difference=float(np.max(abs(chain['psi_rate']-baseline_chain['psi_rate']))),
                body_pose_load_difference=(chain['body_pose_load']-baseline_chain['body_pose_load']).tolist(),
                translation_residual=gradient.sum(axis=0).tolist(),
                kernel_attributes=kernel.attributes))
            gradients.append(gradient)
            print(json.dumps(report['rows'][-1]),flush=True)
        report['complete']=True
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('states','geometry','frozen_report','output'):parser.add_argument('--'+name.replace('_','-'),required=True)
    print(json.dumps(run(**vars(parser.parse_args())),indent=2))
