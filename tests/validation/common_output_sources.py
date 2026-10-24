"""Fixed output targets isolate archived source-resolution increments; no promotion."""
from time import perf_counter
STARTED=perf_counter()
import argparse,json,shutil
from io import BytesIO
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
from importlib import import_module
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from .common import prepare,save,read_owned,digest,current_sources,mesh_quality
from .paired_body import scaled_check,capytaine_radius_check
from .reference_mesh import verify_surface
from .frozen_source_quadrature import evaluate_trace,owned_points,pressure_matrix,mirrored_gauss_targets


def source_quality(values,quadrature_uncertainty,scale,relative=.001,absolute=1e-6):
    A=np.asarray(values)/scale;U=np.asarray(quadrature_uncertainty)/scale
    if A.shape!=(3,6,6) or U.shape!=A.shape or not np.isfinite(A).all() or not np.isfinite(U).all() or np.any(U<0):raise ValueError('three complete source levels and integration uncertainties required')
    previous=abs(A[1]-A[0]);latest=abs(A[2]-A[1])
    ep=U[0]+U[1];el=U[1]+U[2];allow=absolute+relative*abs(A[2])
    robust_previous=np.maximum(previous-ep,0)
    arrays=dict(source_previous_increment=(A[1]-A[0])*scale,source_latest_increment=(A[2]-A[1])*scale,
        source_increment_envelope=np.maximum(previous,latest)*scale,
        previous_quadrature_charge=ep*scale,latest_quadrature_charge=el*scale,
        source_previous_ratio=(previous+ep)/allow,source_latest_ratio=(latest+el)/allow,
        source_trend_ratio=np.maximum(latest+el-robust_previous,0)/absolute,
        source_trend_failed=latest+el>robust_previous+absolute)
    passed=bool(np.all(arrays['source_previous_ratio']<=1)&np.all(arrays['source_latest_ratio']<=1)&~np.any(arrays['source_trend_failed']))
    return arrays,dict(passed=passed,previous_max_ratio=float(arrays['source_previous_ratio'].max()),
        latest_max_ratio=float(arrays['source_latest_ratio'].max()),robust_trend_passed=not bool(np.any(arrays['source_trend_failed'])))


def pinned_bytes(bundle,stem,path):
    def stream_hash(stream):
        h=sha256()
        while chunk:=stream.read(1024*1024):h.update(chunk)
        return h.hexdigest()
    with ZipFile(bundle) as z:
        raw=path.with_suffix('.json').read_bytes();report=json.loads(raw)
        if z.read('validation/'+stem+'.json')!=raw:raise ValueError('input manifest differs from pinned bundle')
        with path.with_suffix('.npz').open('rb') as stream:
            if stream_hash(stream)!=report['arrays_sha256']:raise ValueError('local input arrays differ')
        direct='validation/'+stem+'.npz'
        if direct in z.namelist():
            with z.open(direct) as stream:
                if stream_hash(stream)!=report['arrays_sha256']:raise ValueError('pinned input arrays differ')
        else:
            entry='validation/'+stem+'.payload.zip'
            with z.open(entry) as stream:
                if stream_hash(stream)!=report['payload_package']['sha256']:raise ValueError('pinned payload differs')
            with ZipFile(BytesIO(z.read(entry))) as payload:
                with payload.open('arrays.npz') as stream:
                    if stream_hash(stream)!=report['arrays_sha256']:raise ValueError('nested input arrays differ')


def generate(args,c,r,a):
    fixed=dict(format_version=1,source_levels=[0,1,2],source_counts=[54,216,864],gauss_orders=[12,16,20],
        target_counts=[124416,221184,345600],evaluation_chunk_points=100,cpu_threads=1,maximum_elapsed_seconds=360.,
        density_kg_m3=1025.,gravity_m_s2=9.81,trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,
        reference_relative=.001,reference_absolute_normalized=1e-6,full_reference_scope=False,
        study_bundle_sha256='052df7c50865340e7ef4ebaddfc48adfef191b24057646bd395a889515dc98f2',
        quadrature_bundle_sha256='7e999761722f68b52c9675abfb32bc02aff101c55f03ad76babf749154ab7a63',
        historical_bundle_sha256='28a35794d8c20d35e59f955e0b4ced8fec3cdebd4c6f8e01beffd1cb205fdbb0',
        quadrature_controls_sha256='a691b8e5724689917bcb212599471c3530af2b02a3208bf85f012e764bc337d9')
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('undeclared common-output source controls')
    for name,key in (('study_bundle','study_bundle_sha256'),('quadrature_bundle','quadrature_bundle_sha256'),('historical_bundle','historical_bundle_sha256')):
        if digest(getattr(args,name))!=c[key]:raise ValueError('pinned input bundle differs: '+name)
    p,data=read_owned(args.study,digest('validation/paired-body-controls.json'),'paired-body-study')
    current_sources(p,tuple(p['producer_sha256']))
    q,qa=read_owned(args.quadrature,c['quadrature_controls_sha256'],'frozen-source-pressure-quadrature')
    pinned_bytes(args.study_bundle,args.study.stem,args.study);pinned_bytes(args.quadrature_bundle,args.quadrature.stem,args.quadrature)
    if (q['parent_manifest_sha256']!=digest(args.study) or q['parent_arrays_sha256']!=digest(args.study.with_suffix('.npz'))
            or q['parent_bundle_sha256']!=c['study_bundle_sha256'] or q['reused_bundle_sha256']!=c['historical_bundle_sha256']
            or q['input_sha256']!=p['input_sha256'] or q['trace_normalization']!=p['trace_normalization']
            or q['completed_levels']!=4 or not q['quadrature_stabilization']['passed']
            or [x['rule_parameter'] for x in q['levels']]!=[0,12,16,20] or p['completed_capytaine_calls']!=24
            or q['new_bem_solves']!=0 or q['source_discretization_qualified'] is not False
            or p['dof_order']!=['Surge','Sway','Heave','Roll','Pitch','Yaw']
            or q['full_reference_qualified'] is not False or p['full_reference_qualified'] is not False):raise ValueError('source/quadrature identity or scope differs')
    archive=args.output.with_suffix('.source')
    for source,name in ((args.study_bundle,'study.zip'),(args.quadrature_bundle,'quadrature.zip'),(args.historical_bundle,'historical.zip')):shutil.copyfile(source,archive/name)
    physical=args.study.with_suffix('.source')/'inputs.json';shutil.copyfile(physical,archive/'inputs.json');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'inputs.json')!=p['input_sha256']:raise ValueError('physical input differs')
    mean=SimulationInputs.from_json(physical.read_text());scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:])
    for key,dkey in (('source_triangles','triangles_2'),('normalized_sources','capytaine_sources_2'),('world_cg','world_cg'),('dry_mass','dry_mass')):np.testing.assert_array_equal(qa[key],data[dkey])
    np.testing.assert_array_equal(mean.initial_state.position,data['world_cg'])
    a.update(world_cg=data['world_cg'],dry_mass=data['dry_mass'],target_triangles=qa['source_triangles'])
    r.update(scope=__doc__,completed_sources=0,source_records=[],new_bem_solves=0,full_reference_qualified=False,
        source_discretization_qualified=False,required_four_level_coverage=False,input_sha256=p['input_sha256'],
        study_manifest_sha256=digest(args.study),study_arrays_sha256=digest(args.study.with_suffix('.npz')),
        quadrature_manifest_sha256=digest(args.quadrature),quadrature_arrays_sha256=digest(args.quadrature.with_suffix('.npz')),
        input_bundle_sha256={name:c[name+'_bundle_sha256'] for name in ('study','quadrature','historical')},
        protocol_sha256=digest(c['protocol']),trace_normalization=p['trace_normalization'])
    save(args.output,r,a)
    def budget():
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('common-output source budget exhausted')
    budget();import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        if cpt.__version__!='3.0.0' or any(x['num_threads']!=1 for x in threadpool_info()):raise ValueError('pinned dependency/one native thread required')
        r.update(capytaine_version=cpt.__version__,native_thread_pools=threadpool_info())
        target_body=body_from_triangles(qa['source_triangles'],data['world_cg'],mean.body.mass)
        r['target_radius_error_m']=capytaine_radius_check(target_body,qa['source_triangles'])
        for i,order in enumerate(c['gauss_orders'],1):
            prefix=f'target_{order}_';points,owner,relative,_=mirrored_gauss_targets(qa['source_triangles'],order)
            weights=target_body.mesh.faces_areas[owner]*relative
            for key,value in (('points',points),('owners',owner),('weights',weights)):np.testing.assert_array_equal(value,qa[f'level{i}_'+key])
            if len(points)!=c['target_counts'][i-1]:raise ValueError('missing common targets')
            owned_points(points,owner,qa['source_triangles'],c['ownership_absolute_m'])
            a.update({prefix+'points':points,prefix+'owners':owner,prefix+'weights':weights,prefix+'normal':-target_body.mesh.faces_normals[owner]})
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        if green.exportable_settings!=q['green_settings'] or green.exportable_settings!=p['green_settings']:raise ValueError('Green settings differ')
        r['green_settings']=green.exportable_settings;engine=cpt.BEMSolver(green_function=green).engine
        source_A=[];quad_U=[];quad_pass=[]
        for level,count in zip(c['source_levels'],c['source_counts'],strict=True):
            budget();triangles=data[f'triangles_{level}'];sigma=data[f'capytaine_sources_{level}'];prefix=f'source_{level}_'
            if triangles.shape!=(count,3,3) or sigma.shape!=(count,6):raise ValueError('missing source solution')
            proof=verify_surface(data['parent_body_triangles'],triangles)
            record=dict(level=level,panels=count,surface_proof=proof,reused=level==2,orders=[])
            a.update({prefix+'triangles':triangles,prefix+'normalized_sources':sigma})
            body=body_from_triangles(triangles,data['world_cg'],mean.body.mass);record['radius_error_m']=capytaine_radius_check(body,triangles)
            if level<2:
                native=evaluate_trace(engine,body.mesh,body.mesh.faces_centers,sigma,c['evaluation_chunk_points'],budget)
                native_A=pressure_matrix(body.mesh.faces_centers,-body.mesh.faces_normals,body.mesh.faces_areas,native,data['world_cg'],1025.)
                record['native_potential_ratio']=scaled_check(native,data[f'capytaine_aux_{level}'],1.,c['trace_scaled_tolerance'])
                record['native_load_ratio']=scaled_check(native_A,data['capytaine_Ainf'][level],scale,c['trace_scaled_tolerance'])
                a.update({prefix+'native_potential':native,prefix+'native_A_complex':native_A})
            matrices=[]
            for i,order in enumerate(c['gauss_orders'],1):
                target=f'target_{order}_';key=prefix+f'order_{order}_';begin=perf_counter()
                phi=qa[f'level{i}_potential'] if level==2 else evaluate_trace(engine,body.mesh,a[target+'points'],sigma,c['evaluation_chunk_points'],budget)
                A=pressure_matrix(a[target+'points'],a[target+'normal'],a[target+'weights'],phi,data['world_cg'],1025.)
                imaginary=scaled_check(A.imag,np.zeros((6,6)),scale,c['trace_scaled_tolerance'])
                if level==2:scaled_check(A,qa[f'level{i}_A_complex'],scale,c['trace_scaled_tolerance'])
                else:a[key+'potential']=phi
                a[key+'A_complex']=A;matrices.append(A.real.copy())
                record['orders'].append(dict(order=order,seconds=perf_counter()-begin,imaginary_ratio=imaginary))
                save(args.output,r,a);budget();print(f'common-output sources: source{count}, order{order}',flush=True)
            matrices=np.stack(matrices);U,gate=mesh_quality(matrices/scale,dict(reference_mesh_relative=.001,reference_mesh_absolute_normalized=1e-6))
            record['quadrature_stabilization']=gate
            a[prefix+'quadrature_A']=matrices;a[prefix+'quadrature_uncertainty']=U*scale
            qp=abs(matrices[1]-matrices[0])/scale;ql=abs(matrices[2]-matrices[1])/scale;allow=1e-6+.001*abs(matrices[2]/scale)
            a.update({prefix+'quadrature_previous_ratio':qp/allow,prefix+'quadrature_latest_ratio':ql/allow,
                prefix+'quadrature_trend_ratio':np.maximum(ql-qp,0)/1e-6,
                prefix+'quadrature_previous_failed':qp>allow,prefix+'quadrature_latest_failed':ql>allow,
                prefix+'quadrature_trend_failed':ql>qp+1e-6})
            source_A.append(matrices[-1]);quad_U.append(U*scale);quad_pass.append(gate['passed'])
            r['source_records'].append(record);r['completed_sources']+=1;save(args.output,r,a);budget()
    diagnostics,gate=source_quality(np.stack(source_A),np.stack(quad_U),scale)
    for key,value in diagnostics.items():a[key]=value
    for key in ('source_previous_ratio','source_latest_ratio'):a[key+'_failed']=a[key]>1
    a.update(source_A=np.stack(source_A),quadrature_uncertainty=np.stack(quad_U))
    for key in ('source_previous_increment','source_latest_increment'):a[key+'_absolute']=abs(a[key])
    gate['all_quadrature_passed']=bool(all(quad_pass));gate['passed']=bool(gate['passed'] and all(quad_pass))
    r.update(complete=True,passed=False,source_diagnostic=gate,elapsed_seconds=perf_counter()-STARTED,
        interpretation='Three source solutions, common output targets, empirical quadrature uncertainty charged; not full four-level or finite-wave qualification')
    budget()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/common-output-source-controls.json'))
    for name in ('study','study-bundle','quadrature','quadrature-bundle','historical-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'common-output-source-diagnostic');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a);print(json.dumps({key:r.get(key) for key in ('complete','passed','completed_sources','source_diagnostic','elapsed_seconds','error')}))
    return 0 if r['complete'] and r['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
