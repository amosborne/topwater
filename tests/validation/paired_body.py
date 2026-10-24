"""Controlled paired Ainf-only study; full finite-wave reference remains missing."""
import argparse
from hashlib import sha256
from importlib import import_module
from io import BytesIO
import json
from pathlib import Path
import shutil
from time import perf_counter
from zipfile import ZipFile
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from topwater.inertia import instantaneous_added_mass
from topwater.numerics import LaplaceControls,DenseOperatorPlan
from .common import prepare,save,read_owned,digest,current_sources,midpoint_children,study_preflight,mesh_quality,entrywise
from .reference_mesh import half_surface,farthest_first,verify_surface


def controls(path):
    c=json.loads(Path(path).read_text())
    required=dict(format_version=1,refinements=[0,1,2,3],cpu_threads=1,
        roundoff_coordinate_absolute_m=1e-12,density_kg_m3=1025.,gravity_m_s2=9.81,
        maximum_reference_workspace_bytes=12884901888,maximum_topwater_operator_bytes=536870912,
        maximum_elapsed_seconds=180.,reference_relative=.001,reference_absolute_normalized=1e-6,
        candidate_relative=.01,candidate_absolute_normalized=1e-5,trace_scaled_tolerance=1e-10,
        full_reference_scope=False,symmetry_half='world y>=0 then orientation-preserving reflection')
    if any(c.get(k)!=v for k,v in required.items()):raise ValueError('unsupported prospective paired-body contract')
    return c


def symmetric_family(triangles,refinements=(0,1,2,3),tolerance=1e-12):
    half=half_surface(triangles).copy();snap=abs(half[:,:,1])<=tolerance
    maximum=float(np.max(abs(half[:,:,1][snap]),initial=0.));half[snap,1]=0.
    retained=[];points=half.reshape(-1,3);mapped=[]
    for point in points:
        matches=[i for i,p in enumerate(retained) if np.max(abs(point-p))<=tolerance]
        if matches:
            chosen=retained[matches[0]];maximum=max(maximum,float(np.max(abs(point-chosen))))
        else:retained.append(point.copy());chosen=retained[-1]
        mapped.append(chosen)
    half=np.array(mapped).reshape(half.shape)
    family=[];half_levels=[]
    for level in range(max(refinements)+1):
        if level in refinements:
            reflected=half[:,[0,2,1]]*np.array([1.,-1.,1.])
            full=farthest_first(np.r_[half,reflected]);verify_surface(triangles,full)
            family.append(full);half_levels.append(half.copy())
        if level<max(refinements):half=midpoint_children(half)
    return family,half_levels,maximum


def owned_parent(bundle,c):
    if digest(bundle)!=c['parent_bundle_sha256']:raise ValueError('parent bundle differs')
    stem=c['parent_stem'];base='validation/'+stem
    with ZipFile(bundle) as z:
        raw=z.read(base+'.json');r=json.loads(raw)
        if (r['kind']!='pressure-discrepancy-diagnostic' or not r['complete'] or r['passed']
                or r['full_reference_qualified'] or r['completed_capytaine_calls']!=6
                or r['contract_sha256']!=digest('validation/pressure-discrepancy-controls.json')):
            raise ValueError('parent scope/controls/coverage mismatch')
        if sha256(z.read(base+'.source/contract.json')).hexdigest()!=r['contract_sha256']:
            raise ValueError('parent archived contract differs')
        for name,h in r['producer_sha256'].items():
            if sha256(z.read(base+'.source/'+name)).hexdigest()!=h:raise ValueError('corrupt parent producer '+name)
        current_sources(r,tuple(r['producer_sha256']))
        payload=z.read(base+'.npz')
        if sha256(payload).hexdigest()!=r['arrays_sha256']:raise ValueError('corrupt parent arrays')
        physical=z.read(base+'.source/calm-inputs.json')
        if sha256(physical).hexdigest()!=r['calm_inputs_sha256']:raise ValueError('parent physical input differs')
        with np.load(BytesIO(payload)) as a:data={name:a[name].copy() for name in a.files}
    return r,data,physical,sha256(raw).hexdigest()


def freeze(args,c,r,a):
    parent,data,physical,parent_hash=owned_parent(args.parent_bundle,c)
    archive=args.output.with_suffix('.source');shutil.copyfile(args.parent_bundle,archive/'parent.zip')
    shutil.copyfile(c['protocol'],archive/'protocol.md');(archive/'inputs.json').write_bytes(physical)
    r.update(input_sha256=digest(archive/'inputs.json'),parent_manifest_sha256=parent_hash,
        parent_bundle_sha256=digest(args.parent_bundle),protocol_sha256=digest(c['protocol']),scope=__doc__,
        full_reference_qualified=False)
    save(args.output,r,a)
    mean=SimulationInputs.from_json(physical.decode())
    if mean.initial_state.position[1]!=0 or len(mean.waves.amplitudes):raise ValueError('verified calm symmetric reference required')
    original=hydrostatics(mean.body,mean.initial_state,density=c['density_kg_m3'],gravity=c['gravity_m_s2'])
    verify_surface(original.wetted_triangles,data['body_triangles'])
    family,halves,weld=symmetric_family(data['body_triangles'],c['refinements'],c['roundoff_coordinate_absolute_m'])
    R=rotation_matrix(mean.initial_state.orientation);mass=np.zeros((6,6));mass[:3,:3]=mean.body.mass*np.eye(3);mass[3:,3:]=R@mean.body.inertia_cg@R.T
    a.update(world_cg=mean.initial_state.position,dry_mass=mass,original_wetted_triangles=original.wetted_triangles,
        parent_body_triangles=data['body_triangles'],previous_fine_Ainf=data['fine_reference_Ainf'],
        previous_fine_uncertainty=data['fine_reference_uncertainty'])
    rows=[]
    for level,(triangles,half) in enumerate(zip(family,halves,strict=True)):
        proof=verify_surface(original.wetted_triangles,triangles)
        panels=TrianglePanels(triangles,-1);edges=panels.edge_lengths
        rows.append(dict(level=level,panels=len(triangles),proof=proof,
            maximum_edge_m=float(edges.max()),minimum_area_m2=float(panels.areas.min()),
            maximum_aspect_indicator=float((edges.max(axis=1)**2/(2*panels.areas)).max())))
        a[f'triangles_{level}']=triangles;a[f'half_{level}']=half
    counts=[row['panels'] for row in rows]
    reference=study_preflight(counts[0],c['refinements'],c['maximum_reference_workspace_bytes'],panel_counts=counts)
    top=[DenseOperatorPlan(n,LaplaceControls(max_operator_bytes=c['maximum_topwater_operator_bytes'])).operator_bytes for n in counts]
    if counts!=[counts[0]*4**n for n in c['refinements']]:raise ValueError('unexpected midpoint family counts')
    r.update(complete=True,passed=True,levels=rows,reference_resource_plan=reference,
        topwater_operator_bytes=top,maximum_roundoff_weld_m=weld,
        qualification_scope='geometry/control freeze only; no hydrodynamic solve')


def scaled_check(actual,expected,scale,tolerance):
    x=np.asarray(actual);y=np.asarray(expected)
    if x.shape!=y.shape or not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('missing/nonfinite paired trace')
    ratio=abs(x-y)/(tolerance*(np.asarray(scale)+abs(y)))
    if np.any(ratio>1):raise ValueError('paired source/load consistency failed')
    return float(np.max(ratio))


def entry_diagnostics(candidate,reference,scale,c):
    C=np.asarray(candidate)/scale;A=np.asarray(reference)/scale
    previous=abs(A[-2]-A[-3]);latest=abs(A[-1]-A[-2]);U=np.maximum(previous,latest)
    ref_allow=c['reference_absolute_normalized']+c['reference_relative']*abs(A[-1])
    candidate_allow=c['candidate_absolute_normalized']+c['candidate_relative']*abs(A[-1])
    mesh_allow=c['candidate_absolute_normalized']+c['candidate_relative']*abs(C[-1])
    ratios=dict(reference_previous=previous/ref_allow,reference_latest=latest/ref_allow,
        reference_trend=np.maximum(latest-previous,0)/c['reference_absolute_normalized'],
        candidate_mesh=abs(C[-1]-C[-2])/mesh_allow,
        candidate_comparison=(abs(C[-1]-A[-1])+U)/candidate_allow)
    arrays={name+'_ratio':value for name,value in ratios.items()}
    arrays.update({name+'_failed':value>1 for name,value in ratios.items()})
    return arrays,{name:np.argwhere(value>1).tolist() for name,value in ratios.items()}


def study(args,c,r,a):
    frozen,data=read_owned(args.frozen,digest(args.controls),'paired-body-freeze',passed=True)
    current_sources(frozen,tuple(frozen['producer_sha256']))
    if (len(frozen['levels'])!=4 or [x['level'] for x in frozen['levels']]!=c['refinements']
            or [x['panels'] for x in frozen['levels']]!=[frozen['levels'][0]['panels']*4**n for n in c['refinements']]):
        raise ValueError('all four declared ordered refinement levels required')
    counts=[len(data[f'triangles_{n}']) for n in c['refinements']]
    if counts!=[x['panels'] for x in frozen['levels']]:raise ValueError('frozen actual panel counts differ')
    study_preflight(counts[0],c['refinements'],c['maximum_reference_workspace_bytes'],panel_counts=counts)
    for count in counts:DenseOperatorPlan(count,LaplaceControls(max_operator_bytes=c['maximum_topwater_operator_bytes']))
    src=args.frozen.with_suffix('.source');archive=args.output.with_suffix('.source')
    for suffix in ('.json','.npz'):shutil.copyfile(args.frozen.with_suffix(suffix),archive/('freeze'+suffix))
    for source,target in ((src/'parent.zip','parent.zip'),(src/'inputs.json','inputs.json'),(src/'protocol.md','protocol.md')):
        shutil.copyfile(source,archive/target)
    if (digest(archive/'parent.zip')!=c['parent_bundle_sha256'] or digest(archive/'inputs.json')!=frozen['input_sha256']
            or digest(archive/'protocol.md')!=frozen['protocol_sha256']):raise ValueError('frozen actual inputs differ')
    mean=SimulationInputs.from_json((archive/'inputs.json').read_text())
    r.update(freeze_manifest_sha256=digest(args.frozen),input_sha256=frozen['input_sha256'],
        parent_bundle_sha256=c['parent_bundle_sha256'],protocol_sha256=frozen['protocol_sha256'],
        completed_capytaine_calls=0,completed_topwater_levels=0,levels=[],trace_checks=[],
        full_reference_qualified=False,scope=__doc__)
    a.update(data);save(args.output,r,a);started=perf_counter()
    def budget():
        if perf_counter()-started>c['maximum_elapsed_seconds']:raise TimeoutError('paired study elapsed budget exhausted')
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles,MODES
        if cpt.__version__!='3.0.0' or any(p['num_threads']!=1 for p in threadpool_info()):raise ValueError('pinned Capytaine and one native thread required')
        r['capytaine_version']=cpt.__version__;r['native_thread_pools']=threadpool_info()
        r.update(time_convention='Re(amplitude*exp(-i*omega*t))',
            trace_normalization='Capytaine phi/sources/bc divided by -i*omega before symbolic infinity cancellation; unit CG acceleration potentials',
            auxiliary_potential_units=['m','m','m','m^2','m^2','m^2'],
            trace_normals='Topwater fluid-outward; Capytaine body-outward; physical potential unchanged')
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        r['green_settings']=green.exportable_settings
        top=[];reference=[];rho=c['density_kg_m3'];scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:])
        for level,row in enumerate(frozen['levels']):
            budget();triangles=data[f'triangles_{level}'];panels=TrianglePanels(triangles,-1)
            N=np.c_[panels.normals,np.cross(panels.centroids-data['world_cg'],panels.normals)]
            begin=perf_counter();t=instantaneous_added_mass(panels,data['world_cg'],density=rho,
                controls=LaplaceControls(max_operator_bytes=c['maximum_topwater_operator_bytes']))
            budget();top.append(t.added_mass.copy());a[f'topwater_aux_{level}']=t.acceleration_potential
            ratio=scaled_check(rho*N.T@(panels.areas[:,None]*t.acceleration_potential),t.added_mass,scale,c['trace_scaled_tolerance'])
            record=dict(level=level,panels=len(triangles),complete=False,topwater_seconds=perf_counter()-begin,
                topwater_load_ratio=ratio,topwater_diagnostics=t.diagnostics)
            r['levels'].append(record);r['completed_topwater_levels']+=1
            body=body_from_triangles(triangles,data['world_cg'],mean.body.mass)
            radius_error=float(abs(body.mesh.faces_radii-np.linalg.norm(triangles-body.mesh.faces_centers[:,None],axis=2).max(axis=1)).max())
            if radius_error>1e-12:raise ValueError('Capytaine panel radii differ from farthest-vertex control')
            record['maximum_radius_error_m']=radius_error
            solver=cpt.BEMSolver(green_function=green)
            if solver.engine.exportable_settings['linear_solver']!='lu_decomposition':raise ValueError('non-overwriting default LU required')
            phi=np.zeros((len(triangles),6),complex);sigma=phi.copy();q=phi.copy();A=np.zeros((6,6))
            for name,value in (('capytaine_aux',phi),('capytaine_sources',sigma),('capytaine_bc',q),('capytaine_A',A)):
                a[f'{name}_{level}']=value
            for j,mode in enumerate(MODES):
                budget();pb=cpt.RadiationProblem(body=body,radiating_dof=mode,omega=np.inf,
                    rho=rho,g=c['gravity_m_s2'],water_depth=np.inf,forward_speed=0.)
                if j==0:
                    S,K=solver.engine.build_matrices(body.mesh,body.mesh,free_surface=pb.free_surface,
                        water_depth=pb.water_depth,wavenumber=pb.wavenumber,
                        adjoint_double_layer=True,diagonal_term_in_double_layer=True)
                    budget()
                result=solver.solve(pb,method='indirect',keep_details=True);r['completed_capytaine_calls']+=1;budget()
                factor=-1j*pb.encounter_omega
                phi[:,j]=np.asarray(result.potential/factor);sigma[:,j]=np.asarray(result.sources/factor)
                q[:,j]=np.asarray(pb.boundary_condition/factor)
                A[:,j]=[result.added_mass[m] for m in MODES]
                if any(result.radiation_damping[m]!=0 for m in MODES):raise ValueError('nonzero Ainf damping')
                save(args.output,r,a)
            ratios=dict(rigid_normal=scaled_check(q,-N,1.,c['trace_scaled_tolerance']),
                potential=scaled_check(S@sigma,phi,1.,c['trace_scaled_tolerance']),
                normal=scaled_check(K@sigma,q,1.,c['trace_scaled_tolerance']),
                load=scaled_check(rho*N.T@(panels.areas[:,None]*phi),A,scale,c['trace_scaled_tolerance']))
            del S,K
            r['trace_checks'].append(dict(level=level,**ratios));reference.append(A.copy())
            record['complete']=True;save(args.output,r,a);budget()
            print(f'paired Ainf: {level+1}/4 levels, {r["completed_capytaine_calls"]}/24 Capytaine calls',flush=True)
        if r['completed_capytaine_calls']!=24 or r['completed_topwater_levels']!=4 or not all(x['complete'] for x in r['levels']):
            raise ValueError('missing required paired physics cases')
        C=np.stack(top);A=np.stack(reference);a.update(topwater_Ainf=C,capytaine_Ainf=A)
        U,gate=mesh_quality(A/scale,dict(reference_mesh_relative=c['reference_relative'],reference_mesh_absolute_normalized=c['reference_absolute_normalized']))
        candidate=entrywise(C[-1]/scale,A[-1]/scale,U,c['candidate_absolute_normalized'],c['candidate_relative'])
        mesh=entrywise(C[-2]/scale,C[-1]/scale,np.zeros((6,6)),c['candidate_absolute_normalized'],c['candidate_relative'])
        a.update(capytaine_uncertainty_normalized=U,capytaine_uncertainty_physical=U*scale,
            interfamily_fine_difference_normalized=(A[-1]-data['previous_fine_Ainf'][-1])/scale,
            capytaine_previous_increment=abs(A[-2]-A[-3]),
            capytaine_latest_increment=abs(A[-1]-A[-2]),candidate_difference_normalized=(C[-1]-A[-1])/scale,
            candidate_latest_increment=abs(C[-1]-C[-2]),topwater_same_mesh_difference_normalized=(C-A)/scale)
        entries,failures=entry_diagnostics(C,A,scale,c);a.update(entries)
        parity=np.array([1,-1,1,-1,1,-1]);odd=parity[:,None]!=parity[None,:]
        r.update(complete=True,passed=False,elapsed_seconds=perf_counter()-started,
            added_mass_mesh_gate_passed=gate['passed'],instantaneous_body_gates_passed=bool(gate['passed'] and candidate['passed'] and mesh['passed']),
            added_mass_reference_qualified=False,
            gates=dict(reference=gate,candidate_mesh=mesh,candidate_comparison=candidate),
            dof_order=list(MODES),failed_entries=failures,
            maximum_odd_parity_normalized=dict(capytaine=float(np.max(abs(A/scale)[:,odd])),topwater=float(np.max(abs(C/scale)[:,odd]))),
            qualification_status='Ainf-only study; full finite-wave reference, domain/initialization and TD simulation remain unqualified')
        budget()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=('freeze','study'))
    p.add_argument('--controls',type=Path,default=Path('validation/paired-body-controls.json'))
    p.add_argument('--parent-bundle',type=Path);p.add_argument('--frozen',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();c=controls(args.controls);r=prepare(args.output,args.controls,'paired-body-'+args.stage);a={};save(args.output,r,a)
    try:
        if args.stage=='freeze':freeze(args,c,r,a)
        else:study(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a)
    print(json.dumps({k:r.get(k) for k in ('complete','passed','levels','gates','completed_capytaine_calls','error')}))
    return 0 if r['complete'] and r['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
