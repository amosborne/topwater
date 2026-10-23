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
from . import cuda_reverse_gradient
from .cuda_pair_gradient_probe import selected_pair_derivatives
from .boat_chart_tangent import chart_tangent
from .boat_action_tangent import refine_tangent,action_terms
from .boat_world_field import current_field
from .incident_baseline_benchmark import coefficients


def run(frozen_report,output,states=None,geometry=None,reverse=False,order=3,contact_snapshot=None):
    frozen_report,output=map(Path,(frozen_report,output))
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
    contact=None
    if contact_snapshot is not None:
        contact_snapshot=Path(contact_snapshot);contact=json.loads(contact_snapshot.read_text())
        if not contact['complete'] or contact.get('extra_free_refinement',0):
            raise ValueError('complete body-only refined contact snapshot required')
        root=contact_snapshot.resolve().parent.parent
        def resolve(path):
            path=Path(path);return path if path.is_absolute() else root/path
        for path,digest in contact['input_sha256'].items():
            if sha256(resolve(path).read_bytes()).hexdigest()!=digest:
                raise ValueError('saved contact input changed: '+path)
        geometry=resolve(contact['chart_identity']['path']).with_suffix('.json')
        states=contact_snapshot.with_suffix('.npz');contact_data=dict(np.load(states,allow_pickle=False))
        state=contact_data['state'];time=float(contact['time_seconds'])
        pilot_path=resolve(next(path for path in contact['input_sha256'] if path.endswith('.json') and '-geometry-' not in path))
        refinement=json.loads(pilot_path.read_text())['boundary_refinement']['levels']
        local_refinement=contact.get('local_body_refinement')
    else:
        if states is None or geometry is None:raise ValueError('states/geometry or contact snapshot required')
        states,geometry=Path(states),Path(geometry)
        saved=np.load(states,allow_pickle=False);state=saved['accepted_states'][0];time=float(saved['accepted_times'][0])
        refinement=1;local_refinement=None
    reference=dict(np.load(geometry.with_suffix('.npz'),allow_pickle=False))
    metadata=json.loads(geometry.read_text())
    n=len(reference['parent_nodes']);twist=state[7:13]
    raw,tangent=chart_tangent(reference,pose=InitialState(state[:3],state[3:7]/np.linalg.norm(state[3:7])),
        parent_eta=state[13:13+n],parent_psi=state[13+n:])
    arrays,dx,ds=refine_tangent(raw,tangent['point_jacobian'],tangent['parent_map_jacobian'],refinement,
        local_body_refinement=local_refinement)
    if contact is not None:
        if str(reference['node_chart_key'])!=contact['chart_identity']['node_chart_key']:
            raise ValueError('contact chart identity changed')
        for name in ('points','faces','free_nodes'):np.testing.assert_array_equal(arrays[name],contact_data[name])
    points,faces=arrays['points'],arrays['faces']
    values=coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    folder=Path(__file__).parent
    names=('cuda_pair.cuh','cuda_scalar_shape.cu','cuda_scalar_gradient.py','cuda_shape_gradient.py',
        'cuda_shape.cu','cuda_dual.cuh','cuda_wide_dual.cuh','cuda_operator.py','cuda_energy.cu',
        'boat_scalar_gradient_audit.py','boat_action_tangent.py','boat_chart_tangent.py','boat_world_field.py',
        'incident_baseline_benchmark.py','cuda_reverse_gradient.py','cuda_scalar_reverse.cu','cuda_pair_gradient_probe.py')
    report=dict(complete=False,qualified=False,scope='One actual canonical state; frozen pre-refactor source comparison, no physical or speed qualification',
        order=order,nodes=len(points),faces=len(faces),time=time,
        frozen_report=str(frozen_report),frozen_report_sha256=sha256(frozen_report.read_bytes()).hexdigest(),
        frozen_source_sha256=control['source_sha256'],states_sha256=sha256(states.read_bytes()).hexdigest(),
        geometry_sha256=sha256(geometry.with_suffix('.npz').read_bytes()).hexdigest(),
        geometry_json_sha256=sha256(geometry.read_bytes()).hexdigest(),
        source_sha256={name:sha256((folder/name).read_bytes()).hexdigest() for name in names},rows=[],reverse=reverse,
        body_refinement=refinement,local_body_refinement=local_refinement)
    if contact is not None:
        report['contact_snapshot']=str(contact_snapshot)
        report['contact_snapshot_sha256']=sha256(contact_snapshot.read_bytes()).hexdigest()
        report['contact_input_sha256']=contact['input_sha256']
    cp=cuda_operator.cp
    def measured(fn):
        cp.cuda.get_current_stream().synchronize();start=perf_counter();result=fn()
        cp.cuda.get_current_stream().synchronize();return result,perf_counter()-start
    try:
        cp.zeros(1).get()
        old,old_seconds=measured(lambda:old_operator.prepare_energy(faces,TrianglePanels(points[faces],1),order=order))
        new,new_seconds=measured(lambda:cuda_operator.prepare_energy(faces,TrianglePanels(points[faces],1),order=order))
        report.update(old_energy_seconds=old_seconds,new_energy_seconds=new_seconds,
            energy_relative_difference=float(np.linalg.norm(new.energy-old.energy)/np.linalg.norm(old.energy)),
            energy_max_absolute_difference=float(np.max(abs(new.energy-old.energy))),
            mixed_trace_max_absolute_difference=float(np.max(abs(new.mixed_trace.get()-old.mixed_trace.get()))))
        field=current_field(arrays,old.energy,values,twist,time=time)
        phi=field['potential'];gradients=[]
        candidate='reverse' if reverse else 'contracted'
        for label,width in (('old_raw',3),(candidate,1),(candidate,1),('old_raw',3)):
            function=old_shape.vertex_gradient if label=='old_raw' else (cuda_reverse_gradient.vertex_gradient if reverse else cuda_scalar_gradient.vertex_gradient)
            assembly=old if label=='old_raw' else new
            options=dict(assembly=assembly,density=WATER_DENSITY,order=order)
            if label!='reverse':options['batch_width']=width
            gradient,seconds=measured(lambda:function(points,faces,phi,**options))
            if not gradients:baseline=gradient
            direction=np.random.default_rng(71).normal(size=points.shape)
            chain=action_terms(arrays,field,gradient,values,twist,time,dx,ds)
            if not gradients:baseline_chain=chain
            kernel=old_shape._kernels[width] if label=='old_raw' else (cuda_reverse_gradient._kernel if reverse else cuda_scalar_gradient._kernels[width])
            report['rows'].append(dict(label=label,batch_width=None if label=='reverse' else width,seconds=seconds,
                relative_gradient_difference=float(np.linalg.norm(gradient-baseline)/np.linalg.norm(baseline)),
                maximum_absolute_gradient_difference=float(np.max(abs(gradient-baseline))),
                random_direction_contraction_difference=float(np.sum((gradient-baseline)*direction)),
                parent_rate_max_difference=float(np.max(abs(chain['psi_rate']-baseline_chain['psi_rate']))),
                body_pose_load_difference=(chain['body_pose_load']-baseline_chain['body_pose_load']).tolist(),
                translation_residual=gradient.sum(axis=0).tolist(),
                kernel_attributes=kernel.attributes))
            gradients.append(gradient)
            print(json.dumps(report['rows'][-1]),flush=True)
        if reverse:
            body,free,wall=(arrays[key] for key in ('body_faces','free_faces','wall_faces'))
            tri=points[faces]
            normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
            normal/=np.linalg.norm(normal,axis=1)[:,None]
            def shared_edge(first,second,coplanar_only=False):
                for a in first:
                    for b in second:
                        if (a!=b and len(np.intersect1d(faces[a],faces[b]))==2
                                and (not coplanar_only or normal[a]@normal[b]>1-1e-12)):
                            return int(a),int(b)
                raise ValueError('canonical pair diagnostic requires a shared edge')
            contact=shared_edge(body,free);coplanar=shared_edge(wall,wall,True)
            area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
            edge=np.roll(tri,-1,axis=1)-tri
            quality=area/np.max(np.sum(edge*edge,axis=2),axis=1)
            slender=int(free[np.argmin(quality[free])])
            center=tri.mean(axis=1)
            nearby=int(body[np.argmin(np.linalg.norm(center[body]-center[slender],axis=1))])
            selected=[(int(body[0]),int(body[0])),(slender,slender),contact,contact[::-1],coplanar,
                      coplanar[::-1],(nearby,slender),(slender,nearby)]
            labels=['body-self','slender-free-self','contact-edge','contact-edge-reverse',
                    'wall-edge','wall-edge-reverse','near-slender','near-slender-reverse']
            forward_pairs,reverse_pairs=selected_pair_derivatives(points,faces,phi,np.array(selected),
                assembly=new,density=WATER_DENSITY,order=order)
            report['selected_pairs']=[dict(label=label,faces=list(pair),
                face_quality=quality[list(pair)].tolist(),
                forward=forward_pairs[k].tolist(),reverse=reverse_pairs[k].tolist(),
                max_absolute_difference=float(np.max(abs(reverse_pairs[k]-forward_pairs[k]))),
                relative_difference=float(np.linalg.norm(reverse_pairs[k]-forward_pairs[k])/max(np.linalg.norm(forward_pairs[k]),1e-30)))
                for k,(label,pair) in enumerate(zip(labels,selected,strict=True))]
        report['complete']=True
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('states','geometry','frozen_report','output'):parser.add_argument('--'+name.replace('_','-'),required=name in ('frozen_report','output'))
    parser.add_argument('--contact-snapshot')
    parser.add_argument('--order',type=int,choices=[3,6,12],default=3)
    parser.add_argument('--reverse',action='store_true')
    print(json.dumps(run(**vars(parser.parse_args())),indent=2))
