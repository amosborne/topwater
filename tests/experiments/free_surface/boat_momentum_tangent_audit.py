"""Bounded exact-current field-rate audit against four fresh offset solves."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
import numpy as np
from scipy.spatial.transform import Rotation
from topwater import InitialState,sea_states
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests.reference_boat import reference_body,WATER_DENSITY,GRAVITY
from .cuda_operator import cp,prepare_energy
from .cuda_shape_gradient import vertex_gradient
from .boat_chart_tangent import chart_tangent
from .boat_geometry_template import evaluate_template,InvalidTemplate
from .boat_action_tangent import refine_tangent,action_terms
from .boat_world_field import current_field,refine_boundary
from .boat_momentum_reference import incident_load_rate
from .boat_momentum_tangent import field_momentum_rate
from .boat_body_action import momentum
from .boat_ale_pressure import differentiate
from .body_pressure import evaluate_body_pressure
from .incident_baseline_benchmark import coefficients
from .boat_incident_tangent import incident_rate as analytic_incident_rate


def run(output,states=None,geometry=None,sample_indices=(0,2),order=3,analytic_incident=False,contact_snapshot=None):
    output=Path(output)
    if output.exists():
        raise FileExistsError(output)
    contact=None
    if contact_snapshot is not None:
        contact_snapshot=Path(contact_snapshot)
        contact=json.loads(contact_snapshot.read_text())
        if not contact['complete'] or contact.get('extra_free_refinement',0):
            raise ValueError('complete body-only refined contact snapshot required')
        root=contact_snapshot.resolve().parent.parent
        def resolve(path):
            path=Path(path)
            return path if path.is_absolute() else root/path
        geometry=resolve(contact['chart_identity']['path']).with_suffix('.json')
        states=contact_snapshot.with_suffix('.npz')
        contact_saved=dict(np.load(states,allow_pickle=False))
        saved=dict(accepted_times=np.array([contact['time_seconds']]),accepted_states=contact_saved['state'][None])
        sample_indices=(0,)
        pilot_path=resolve(next(path for path in contact['input_sha256'] if path.endswith('.json') and '-geometry-' not in path))
        pilot=json.loads(pilot_path.read_text())
        refinement=pilot['boundary_refinement']['levels']
        local_refinement=contact.get('local_body_refinement')
    else:
        if states is None or geometry is None:
            raise ValueError('states and geometry or a contact snapshot are required')
        states,geometry=Path(states),Path(geometry)
        saved=np.load(states,allow_pickle=False);refinement=1;local_refinement=None
    reference=dict(np.load(geometry.with_suffix('.npz'),allow_pickle=False))
    if contact is not None and str(reference['node_chart_key'])!=contact['chart_identity']['node_chart_key']:
        raise ValueError('contact snapshot chart identity changed')
    metadata=json.loads(geometry.read_text())
    values=coefficients(sea_states.loads_realization(json.dumps(metadata['actual_sea_realization'])))
    folder=Path(__file__).parent
    files=('boat_momentum_tangent.py', 'boat_momentum_reference.py','boat_momentum_tangent_audit.py','cuda_directional.py',
           'cuda_directional.cu','cuda_pair.cuh','boat_chart_tangent.py','boat_action_tangent.py',
           'boat_world_field.py','boat_geometry_template.py','boat_ale_pressure.py','body_pressure.py',
           'incident_baseline_benchmark.py','boat_body_action.py','cuda_shape_gradient.py', 'cuda_panel_geometry.py','cuda_shape.cu',
           'cuda_dual.cuh','cuda_wide_dual.cuh','cuda_operator.py','cuda_energy.cu','boat_incident_tangent.py')
    report=dict(complete=False,qualified=False,scope='Exact-current directional field derivative; archived canonical states, no trajectory',
        states=str(states),geometry=str(geometry),states_sha256=sha256(states.read_bytes()).hexdigest(),
        geometry_sha256=sha256(geometry.with_suffix('.npz').read_bytes()).hexdigest(),
        geometry_json_sha256=sha256(geometry.read_bytes()).hexdigest(),order=order,
        analytic_incident=analytic_incident,body_refinement=refinement,
        source_sha256={name:sha256((folder/name).read_bytes()).hexdigest() for name in files},samples=[])
    body=reference_body()
    if contact is not None:
        report['contact_snapshot']=str(contact_snapshot)
        report['contact_snapshot_sha256']=sha256(contact_snapshot.read_bytes()).hexdigest()
        report['local_body_refinement']=local_refinement

    def phase(name,fn,timings):
        cp.cuda.get_current_stream().synchronize();start=perf_counter();result=fn()
        cp.cuda.get_current_stream().synchronize();timings[name]=perf_counter()-start
        return result

    try:
        cp.zeros(1).get()
        for index in sample_indices:
            time=float(saved['accepted_times'][index]);state=saved['accepted_states'][index]
            timings={};sample=dict(index=int(index),time=time,timings=timings,controls=[])
            report['samples'].append(sample)
            pose=InitialState(state[:3],state[3:7]/np.linalg.norm(state[3:7]));n=len(reference['parent_nodes'])
            twist,eta,psi=state[7:13],state[13:13+n],state[13+n:]
            raw,tangent=phase('chart_tangent',lambda:chart_tangent(reference,pose=pose,parent_eta=eta,parent_psi=psi),timings)
            arrays,dx,ds=phase('refine_tangent',lambda:refine_tangent(raw,tangent['point_jacobian'],tangent['parent_map_jacobian'],
                refinement,local_body_refinement=local_refinement),timings)
            if contact is not None:
                for name in ('points','faces','free_nodes'):
                    np.testing.assert_array_equal(arrays[name],contact_saved[name])
            points,faces=arrays['points'],arrays['faces'];sample.update(nodes=len(points),faces=len(faces))
            assembly=phase('center_energy',lambda:prepare_energy(faces,TrianglePanels(points[faces],1),order=order),timings)
            current=phase('center_field',lambda:current_field(arrays,assembly.energy,values,twist,time=time),timings)
            gradient=phase('center_gradient',lambda:vertex_gradient(points,faces,current['potential'],assembly=assembly,density=WATER_DENSITY,order=order),timings)
            terms=phase('action_terms',lambda:action_terms(arrays,current,gradient,values,twist,time,dx,ds),timings)
            eta_rate,psi_rate=terms['eta_rate'],terms['psi_rate']
            before=perf_counter()
            for halving in range(13):
                h=1e-5/2**halving;offsets={}
                try:
                    for dt in (h,-h,h/2,-h/2):
                        q=(Rotation.from_rotvec(dt*twist[3:])*Rotation.from_quat(pose.orientation,scalar_first=True)).as_quat(scalar_first=True)
                        moved,_=evaluate_template(reference,pose=InitialState(pose.position+dt*twist[:3],q),
                            parent_eta=eta+dt*eta_rate,parent_psi=psi+dt*psi_rate)
                        for _ in range(refinement):
                            moved=refine_boundary(moved,region='body')
                        for _ in range((local_refinement or {}).get('levels',0)):
                            node=local_refinement['contact_node']
                            selected=moved['body_faces'][np.any(moved['faces'][moved['body_faces']]==node,axis=1)]
                            moved=refine_boundary(moved,region='body',marked_faces=selected)
                        offsets[dt]=moved
                        for key in ('faces','free_nodes','free_faces','body_faces','wall_faces','bottom_faces',
                                    'body_source_group','free_source_parent','free_parent_indices','node_chart_key'):
                            if not np.array_equal(offsets[dt][key],arrays[key]):
                                raise ValueError('exact field-rate oracle left its current semantic chart/basis')
                        if (offsets[dt]['points'].shape!=points.shape
                                or offsets[dt]['free_parent_weights'].shape!=arrays['free_parent_weights'].shape
                                or not np.isfinite(offsets[dt]['points']).all()
                                or not np.isfinite(offsets[dt]['free_parent_weights']).all()):
                            raise ValueError('invalid exact field-rate oracle geometry/map')
                except InvalidTemplate:
                    continue
                break
            else:
                raise InvalidTemplate('no common exact field-rate oracle stencil')
            timings['four_offset_geometries']=perf_counter()-before;sample['step_seconds']=h
            incident={}
            before=perf_counter()
            for step in (h,h/2):
                incident[step]=incident_load_rate(arrays,offsets[step],offsets[-step],values,twist,time,step)
            timings['four_offset_boundary_loads']=perf_counter()-before
            sample['incident_h_h2_max_difference']=float(np.max(abs(incident[h]-incident[h/2])))
            chosen_incident=incident[h/2]
            if analytic_incident:
                moving=np.einsum('nid,d->ni',dx,np.r_[twist,eta_rate])
                chosen_incident=phase('analytic_incident_first_call',lambda:analytic_incident_rate(arrays,moving,values,time),timings)
                warmed=phase('analytic_incident_warm_call',lambda:analytic_incident_rate(arrays,moving,values,time),timings)
                np.testing.assert_array_equal(chosen_incident,warmed)
                sample['analytic_incident_h_h2_max_errors']=[float(np.max(abs(chosen_incident-incident[step]))) for step in (h,h/2)]
            candidate=phase('directional_field_momentum_rate',lambda:field_momentum_rate(arrays,assembly,current,values,twist,time,
                eta_rate,psi_rate,dx,ds,terms['object_derivatives']['body_load_jacobian'],chosen_incident,body_pose_load=terms['body_pose_load']),timings)
            sample['candidate_diagnostics']=candidate['diagnostics']
            fields={}
            before=perf_counter()
            for dt,moved in offsets.items():
                fresh=prepare_energy(moved['faces'],TrianglePanels(moved['points'][moved['faces']],1),order=order)
                fields[dt]=current_field(moved,fresh.energy,values,twist,time=time+dt,acceleration_columns=False)
            cp.cuda.get_current_stream().synchronize();timings['four_offset_energies_fields']=perf_counter()-before
            rotation=rotation_matrix(pose.orientation);inertia=rotation@body.inertia_cg@rotation.T
            rigid=np.zeros((6,6));rigid[:3,:3]=body.mass*np.eye(3);rigid[3:,3:]=inertia
            added=WATER_DENSITY*current['body_load'].T@current['acceleration_potentials']
            external=np.r_[0.,0.,-body.mass*GRAVITY,-np.cross(twist[3:],inertia@twist[3:])]
            def pressure(rate,velocity):
                return evaluate_body_pressure(points,faces,arrays['body_faces'],pose.position,body.mass,inertia,
                    twist,current['potential'],rate,velocity,current['acceleration_potentials'],density=WATER_DENSITY,gravity=GRAVITY)
            candidate_pressure=pressure(candidate['potential_rate_ale'],candidate['mesh_velocity'])
            sample['candidate_known_load']=candidate['known_load'].tolist()
            sample['candidate_pressure_known_load']=candidate_pressure.known_pressure_load.tolist()
            phi_rates=[];load_controls=[]
            for step in (h,h/2):
                rate,velocity,ale=differentiate(arrays,offsets[step],offsets[-step],fields[step]['potential'],
                    fields[-step]['potential'],step,pose.position,twist)
                transport=(momentum(fields[step])-momentum(fields[-step]))/(2*step)
                load=terms['body_pose_load']+transport+candidate['angular_transport']
                oracle_pressure=pressure(rate,velocity)
                phi_rates.append(rate);load_controls.append(load)
                sample['controls'].append(dict(step_seconds=step,
                    potential_rate_max_error=float(np.max(abs(candidate['potential_rate_ale']-rate))),
                    potential_rate_relative_error=float(np.linalg.norm(candidate['potential_rate_ale']-rate)/np.linalg.norm(rate)),
                    mesh_velocity_max_error=float(np.max(abs(candidate['mesh_velocity']-velocity))),
                    momentum_rate_error=(candidate['momentum_rate']-transport).tolist(),
                    known_load=load.tolist(),known_load_error=(candidate['known_load']-load).tolist(),
                    pressure_known_load_error=(candidate_pressure.known_pressure_load-oracle_pressure.known_pressure_load).tolist(),
                    acceleration_error=(np.linalg.solve(rigid+added,candidate['known_load']+external)-np.linalg.solve(rigid+added,load+external)).tolist(),ale=ale))
            sample['oracle_h_h2_potential_rate_max_difference']=float(np.max(abs(phi_rates[0]-phi_rates[1])))
            sample['oracle_h_h2_known_load_difference']=(load_controls[0]-load_controls[1]).tolist()
            print(json.dumps(sample),flush=True)
        report['complete']=True
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('states','geometry','output'):parser.add_argument('--'+name,required=name=='output')
    parser.add_argument('--contact-snapshot')
    parser.add_argument('--analytic-incident',action='store_true')
    parser.add_argument('--sample-indices',nargs='+',type=int,default=[0,2])
    parser.add_argument('--order',type=int,default=3)
    print(json.dumps(run(**vars(parser.parse_args())),indent=2))
