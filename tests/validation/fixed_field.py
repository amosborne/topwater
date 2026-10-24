"""Cold fixed2772 infinite-frequency field integration and bounded continuation."""
from time import perf_counter
STARTED=perf_counter()
import argparse,json,shutil
from pathlib import Path
from importlib import import_module
from importlib.metadata import distribution,version
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from topwater.rankine import TrianglePanels
from .common import prepare,read_owned,digest,current_sources,PHYSICAL,mesh_quality
from .common_output_sources import pinned_bytes
from .common_output_preflight import available_memory
from .finest_cost_probe import owned_zip,process_memory
from .finest_consumer import checkpoint_arrays,write_report,verify_arrays,package_payload,consume_orders
from .frozen_source_quadrature import evaluate_trace,pressure_matrix,owned_points,mirrored_gauss_targets
from .paired_body import scaled_check,capytaine_radius_check
from .reference_mesh import verify_surface


def source_interval(coarse,fine,coarse_u,fine_u,scale):
    values=[np.asarray(v) for v in (coarse,fine,coarse_u,fine_u,scale)]
    if any(v.shape!=(6,6) or not np.isfinite(v).all() for v in values) or np.any(coarse_u<0) or np.any(fine_u<0) or np.any(scale<=0):
        raise ValueError('complete finite source interval and nonnegative uncertainties required')
    allowance=1e-6*scale+.001*abs(fine)
    charged=abs(fine-coarse)+coarse_u+fine_u
    return dict(source_difference=fine-coarse,charged_difference=charged,source_interval_ratio=charged/allowance,source_interval_failed=charged>allowance)


def continuation_plan(previous):
    pairs=2772*308*(24**2+28**2)
    rate=max(row['seconds']/(2772*row['points']) for row in previous['orders'])
    field_bytes=sum(308*q*q for q in (20,24,28))*6*16
    target_bytes=sum(308*q*q for q in (20,24,28))*64
    peak=previous['latest_memory']['peak_committed_bytes']
    return dict(new_field_pairs=pairs,maximum_measured_seconds_per_pair=rate,
        projected_upper_seconds=(rate*pairs+previous['native_seconds'])*1.25+60.,
        previous_peak_committed_bytes=peak,retained_field_bytes=field_bytes,target_bytes=target_bytes,
        selected_input_allowance_bytes=32*1024**2,io_allowance_bytes=128*1024**2,
        estimated_peak_bytes=peak+field_bytes+target_bytes+160*1024**2)


def generate(args,c,r,a):
    continuation=c.get('continuation',False)
    fixed=dict(format_version=1,maximum_elapsed_seconds=900.,maximum_peak_bytes=2147483648,maximum_additional_disk_bytes=2147483648,
        source_level=2,source_panels=2772,gauss_orders=[12,16,20],target_counts=[44352,78848,123200],cpu_threads=1,chunk_points=100,
        trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,reference_relative=.001,reference_absolute_normalized=1e-6,
        source_bundle_sha256='fa306c4fd38fbda08debad105816aaf1e21a270617ae15d00d9f524773d0fc31',full_reference_scope=False)
    if continuation:
        fixed.update(continuation=True,maximum_elapsed_seconds=1500.,gauss_orders=[20,24,28],target_counts=[123200,177408,241472],
            previous_bundle_sha256='5c24ecbfba6a8ff839c77e5fb33085b9f729c72b9c108b329a950050e4523924')
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('undeclared cold2772 controls')
    if digest(args.source_bundle)!=c['source_bundle_sha256']:raise ValueError('source packet differs')
    basis=['source_triangles','native_points','native_potential','normalized_sources','native_Ainf']
    selected=['world_cg','dry_mass']+['level2_'+name for name in basis]
    source,sd=read_owned(args.source,digest('validation/fixed-source-export-controls.json'),'fixed-seed-infinite-source-export',selected=selected)
    pinned_bytes(args.source_bundle,args.source.stem,args.source);current_sources(source,PHYSICAL)
    if not source['source_export_complete'] or source['completed_capytaine_calls']!=18 or source['total_source_cases_represented']!=24 or source['covered_source_levels']!=[1,2,3,4]:raise ValueError('complete four-level source proof required')
    src=args.source.with_suffix('.source');archive=args.output.with_suffix('.source')
    for name,h in source['input_bundle_sha256'].items():
        if digest(src/(name+'.zip'))!=h:raise ValueError('physical input bundle differs')
        shutil.copyfile(src/(name+'.zip'),archive/(name+'.zip'))
    for name in ('inputs.json','preflight.zip'):shutil.copyfile(src/name,archive/name)
    shutil.copyfile(args.source_bundle,archive/'source-export.zip');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'inputs.json')!=source['input_sha256'] or digest(archive/'preflight.zip')!=source['preflight_bundle_sha256']:raise ValueError('physical/preflight identity differs')
    pilot,pd=read_owned(args.pilot,digest('validation/fixed-infinity-pilot-controls.json'),'fixed-seed-infinity-pilot',selected=['world_cg','dry_mass','quadrature_A','quadrature_uncertainty'])
    pinned_bytes(src/'pilot.zip',args.pilot.stem,args.pilot)
    if digest(src/'pilot.zip')!=source['pilot_bundle_sha256'] or digest(args.pilot)!=source['pilot_manifest_sha256'] or not pilot['pilot_passed'] or not pilot['source_export_verified']:raise ValueError('coarse conditional integration proof differs')
    for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(sd[key],pd[key])
    previous=None;prior=None
    if continuation:
        if args.previous is None or args.previous_bundle is None:raise ValueError('exact previous field packet required')
        if digest(args.previous_bundle)!=c['previous_bundle_sha256']:raise ValueError('previous packet differs')
        selected=['world_cg','dry_mass','source_triangles','normalized_sources','target_triangles','order_20_potential','order_20_A_complex']+[f'target_20_{name}' for name in ('points','owners','weights','normal')]
        previous,prior=read_owned(args.previous,digest('validation/fixed-field-2772-controls.json'),'fixed-seed-cold-field',selected=selected)
        pinned_bytes(args.previous_bundle,args.previous.stem,args.previous);current_sources(previous,PHYSICAL)
        if (previous['source_manifest_sha256']!=digest(args.source) or previous['source_bundle_sha256']!=digest(args.source_bundle)
                or previous['input_sha256']!=source['input_sha256'] or previous['completed_orders']!=3 or not previous['native_calibration_passed']
                or previous['completed_capytaine_calls']!=0 or previous['integration_passed'] or previous['green_settings']!=source['green_settings']):raise ValueError('previous failed integration scope differs')
        for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(prior[key],sd[key])
        for key in ('source_triangles','normalized_sources'):np.testing.assert_array_equal(prior[key],sd['level2_'+key])
        plan=continuation_plan(previous)
        if plan['projected_upper_seconds']>c['maximum_elapsed_seconds'] or plan['estimated_peak_bytes']>c['maximum_peak_bytes']:raise ValueError('measured continuation plan exceeds controls')
        shutil.copyfile(args.previous_bundle,archive/'previous.zip')
        r.update(continuation=True,previous_manifest_sha256=digest(args.previous),previous_bundle_sha256=digest(args.previous_bundle),
            previous_integration_passed=False,resource_model=plan)
    elif args.previous is not None or args.previous_bundle is not None:raise ValueError('prior field reuse requires declared continuation controls')
    admitted=[order for order in c['gauss_orders'] if order in (12,16,20)]
    keys=['wetted_triangles','world_cg','dry_mass','reference_seed_triangles','reference_triangles_2']+[f'target_{order}_{name}' for order in admitted for name in ('points','owners','weights','normal')]
    preflight,data=owned_zip(archive/'preflight.zip','fixed-asset-preflight-002','fixed-seed-asset-preflight',digest('validation/fixed-asset-preflight-controls.json'),keys)
    if not preflight['assets_usable'] or preflight['input_sha256']!=source['input_sha256']:raise ValueError('preflight physical identity differs')
    for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(sd[key],data[key])
    triangles=sd['level2_source_triangles'];np.testing.assert_array_equal(triangles,data['reference_triangles_2'])
    if triangles.shape!=(2772,3,3):raise ValueError('wrong native source level')
    mean=SimulationInputs.from_json((archive/'inputs.json').read_text());np.testing.assert_array_equal(mean.initial_state.position,sd['world_cg'])
    package=distribution('capytaine')
    if version('capytaine')!='3.0.0':raise ValueError('pinned native version required')
    identities=source['dependency_sha256']
    for name,h in identities.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('native dependency changed')
    host=dict(available_physical_bytes=available_memory(),available_disk_bytes=shutil.disk_usage(Path.cwd()).free)
    if host['available_physical_bytes']<1.25*c['maximum_peak_bytes'] or host['available_disk_bytes']<2*c['maximum_additional_disk_bytes']:raise ValueError('host margins fail')
    r.update(scope=__doc__,resource_only=False,synthetic_fields=False,input_sha256=source['input_sha256'],input_bundle_sha256=source['input_bundle_sha256'],
        source_manifest_sha256=digest(args.source),source_bundle_sha256=digest(args.source_bundle),preflight_bundle_sha256=source['preflight_bundle_sha256'],
        pilot_manifest_sha256=digest(args.pilot),pilot_bundle_sha256=source['pilot_bundle_sha256'],dependency_sha256=identities,protocol_sha256=digest(c['protocol']),
        host_resources=host,completed_orders=0,orders=[],memory_records=[],dof_order=source['dof_order'],trace_normalization=source['trace_normalization'],
        source_level=2,source_panels=2772,native_pair_count=2772**2,
        field_pair_count=2772*sum(count for order,count in zip(c['gauss_orders'],c['target_counts'],strict=True) if not(continuation and order==20)),
        reused_field_pair_count=2772*123200 if continuation else 0)
    if continuation:r['input_bundle_sha256']={**source['input_bundle_sha256'],'previous':digest(args.previous_bundle)}
    a.update(world_cg=sd['world_cg'],dry_mass=sd['dry_mass'],source_triangles=triangles,normalized_sources=sd['level2_normalized_sources'],
        exported_native_potential=sd['level2_native_potential'],exported_native_Ainf=sd['level2_native_Ainf'],
        coarse_quadrature_A=pd['quadrature_A'],coarse_quadrature_uncertainty=pd['quadrature_uncertainty'])
    def budget():
        memory=process_memory();r['latest_memory']=memory
        if max(memory['peak_working_set_bytes'],memory['peak_committed_bytes'])>c['maximum_peak_bytes']:raise MemoryError('cold field process cap exceeded')
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('cold field total budget exceeded')
    def stage(name):
        budget();r['memory_records'].append(dict(stage=name,**process_memory()));write_report(args.output,r)
    checkpoint_arrays(args.output,r,a);stage('inputs')
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        if any(row['num_threads']!=1 for row in threadpool_info()):raise ValueError('one native thread required')
        r['native_thread_pools']=threadpool_info();r['surface_proof']=verify_surface(data['wetted_triangles'],triangles)
        body=body_from_triangles(triangles,a['world_cg'],mean.body.mass);r['radius_error_m']=capytaine_radius_check(body,triangles)
        np.testing.assert_array_equal(body.mesh.faces_centers,sd['level2_native_points'])
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None)
        if green.exportable_settings!=source['green_settings'] or green.exportable_settings!=pilot['green_settings']:raise ValueError('infinite Green settings differ')
        r['green_settings']=green.exportable_settings;engine=cpt.BEMSolver(green_function=green).engine
        begin=perf_counter();native=evaluate_trace(engine,body.mesh,body.mesh.faces_centers,a['normalized_sources'],100,budget)
        native_A=pressure_matrix(body.mesh.faces_centers,-body.mesh.faces_normals,body.mesh.faces_areas,native,a['world_cg'],1025.)
        scale=np.sqrt(a['dry_mass'].diagonal()[:,None]*a['dry_mass'].diagonal()[None,:])
        r.update(native_potential_ratio=scaled_check(native,a['exported_native_potential'],1.,1e-10),native_load_ratio=scaled_check(native_A,a['exported_native_Ainf'],scale,1e-10),
            native_seconds=perf_counter()-begin,native_evaluated_points=len(native),native_calibration_passed=True)
        a.update(native_points=body.mesh.faces_centers,native_normal=-body.mesh.faces_normals,native_weights=body.mesh.faces_areas,native_potential=native,native_A_complex=native_A)
        checkpoint_arrays(args.output,r,a);stage('native_calibration')
        seed=data['reference_seed_triangles'];panels=TrianglePanels(seed,-1)
        if continuation:np.testing.assert_array_equal(seed,prior['target_triangles'])
        for order,count in zip(c['gauss_orders'],c['target_counts'],strict=True):
            key=f'target_{order}_';pts,owners,relative,error=mirrored_gauss_targets(seed,order)
            if len(pts)!=count:raise ValueError('missing common targets')
            for name,value in [('points',pts),('owners',owners),('weights',panels.areas[owners]*relative),('normal',panels.normals[owners])]:
                if order in admitted:np.testing.assert_array_equal(value,data[key+name])
                else:data[key+name]=value
                if continuation and order==20:np.testing.assert_array_equal(value,prior[key+name])
            owned_points(pts,owners,seed,c['ownership_absolute_m']);budget()
        a['target_triangles']=seed
        def fill(order,query,out):
            begin=perf_counter()
            if continuation and order==20:
                np.testing.assert_array_equal(query,prior['target_20_points']);np.copyto(out,prior['order_20_potential'])
            else:evaluate_trace(engine,body.mesh,query,a['normalized_sources'],100,budget,out=out)
            if np.any(abs(out.imag)>1e-10*(1+abs(out.real))):raise ValueError('unexpected complex infinite trace')
            r['orders'].append(dict(order=order,points=len(query),seconds=perf_counter()-begin,reused=bool(continuation and order==20)));r['completed_orders']+=1
        def checkpoint(order):
            if continuation and order==20:np.testing.assert_array_equal(a['order_20_A_complex'],prior['order_20_A_complex'])
            checkpoint_arrays(args.output,r,a);stage(f'order{order}');print(f'cold2772: order{order} complete',flush=True)
        matrices=consume_orders(data,c['gauss_orders'],a,fill,budget,checkpoint)
        scaled_check(matrices.imag,np.zeros_like(matrices.real),scale,1e-10)
        real=matrices.real;U,gate=mesh_quality(real/scale,dict(reference_mesh_relative=.001,reference_mesh_absolute_normalized=1e-6))
        previous=abs(real[1]-real[0])/scale;latest=abs(real[2]-real[1])/scale;allow=1e-6+.001*abs(real[2]/scale)
        a.update(quadrature_A=real,quadrature_uncertainty=U*scale,previous_ratio=previous/allow,latest_ratio=latest/allow,
            previous_failed=previous>allow,latest_failed=latest>allow,trend_failed=latest>previous+1e-6)
        a.update(source_interval(pd['quadrature_A'][-1],real[-1],pd['quadrature_uncertainty'],U*scale,scale))
        r.update(quadrature_stabilization=gate,source_interval_diagnostic=dict(intervals=1,qualified=False,all_entries_within_allowance=bool(not a['source_interval_failed'].any()),maximum_ratio=float(a['source_interval_ratio'].max())))
        checkpoint_arrays(args.output,r,a);verify_arrays(args.output.with_suffix('.npz'),a,r['arrays_sha256']);stage('array_verification')
        r['payload_package']=package_payload(args.output,archive,r,budget);stage('package_verification')
    for name,h in identities.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('dependency changed during fields')
    disk=sum(p.stat().st_size for p in archive.rglob('*') if p.is_file())+args.output.with_suffix('.npz').stat().st_size+args.output.with_suffix('.payload.zip').stat().st_size+65536
    if disk>c['maximum_additional_disk_bytes']:raise OSError('field artifact disk cap exceeded')
    if r['completed_orders']!=3 or not r['native_calibration_passed'] or r['completed_capytaine_calls']!=0:raise ValueError('incomplete field scope')
    budget();r.update(complete=True,passed=False,integration_passed=bool(gate['passed']),accounted_disk_upper_bytes=disk,elapsed_seconds=perf_counter()-STARTED)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/fixed-field-2772-controls.json'))
    for name in ('source','source-bundle','pilot','output'):p.add_argument('--'+name,type=Path,required=True)
    for name in ('previous','previous-bundle'):p.add_argument('--'+name,type=Path)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'fixed-seed-cold-field');a={}
    r.update(completed_capytaine_calls=0,integration_passed=False,native_calibration_passed=False,full_reference_qualified=False,source_discretization_qualified=False,full_run_authorized=False)
    checkpoint_arrays(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,integration_passed=False,error=dict(type=type(error).__name__,message=str(error)))
    write_report(args.output,r);print(json.dumps({key:r.get(key) for key in ('complete','passed','integration_passed','native_calibration_passed','completed_orders','quadrature_stabilization','source_interval_diagnostic','latest_memory','elapsed_seconds','error')}))
    return 0 if r['complete'] and r['integration_passed'] else 1


if __name__=='__main__':raise SystemExit(main())
