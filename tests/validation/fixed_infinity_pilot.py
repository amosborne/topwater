"""Fresh coarsest fixed-seed Ainf export and common-output quadrature pilot."""
from time import perf_counter
STARTED=perf_counter()
import argparse,gc,json,shutil
from pathlib import Path
from zipfile import ZipFile
from importlib import import_module
from importlib.metadata import distribution,version
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from topwater.rankine import TrianglePanels
from .common import prepare,read_owned,digest,mesh_quality,current_sources,PHYSICAL
from .common_output_sources import pinned_bytes
from .common_output_preflight import available_memory
from .finest_cost_probe import owned_zip,process_memory
from .finest_consumer import checkpoint_arrays,write_report,verify_arrays,package_payload,consume_orders
from .frozen_source_quadrature import evaluate_trace,pressure_matrix,owned_points,mirrored_gauss_targets
from .paired_body import scaled_check,capytaine_radius_check
from .reference_mesh import verify_surface
from .ainf_sources import validate_infinite_basis,export_basis


def generate(args,c,r,a):
    fixed=dict(format_version=1,maximum_elapsed_seconds=600.,maximum_peak_bytes=2147483648,maximum_additional_disk_bytes=2147483648,
        source_panels=1232,source_level=1,gauss_orders=[12,16,20],target_counts=[44352,78848,123200],cpu_threads=1,
        chunk_points=100,trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,reference_relative=.001,reference_absolute_normalized=1e-6,
        preflight_bundle_sha256='cd66598c22357706b8517065a886e79742e887500a63c388c559fc70e7fc924f',full_reference_scope=False)
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('undeclared fixed-seed pilot controls')
    if digest(args.preflight_bundle)!=c['preflight_bundle_sha256']:raise ValueError('asset preflight differs')
    p,data=read_owned(args.preflight,digest('validation/fixed-asset-preflight-controls.json'),'fixed-seed-asset-preflight')
    pinned_bytes(args.preflight_bundle,args.preflight.stem,args.preflight);current_sources(p,PHYSICAL)
    if not p['assets_usable'] or not p['fresh_ainf_export_required'] or p['full_run_authorized'] or p['native_evaluations']!=0:raise ValueError('asset scope differs')
    src=args.preflight.with_suffix('.source');archive=args.output.with_suffix('.source')
    for name,h in p['bundle_sha256'].items():
        if digest(src/(name+'.zip'))!=h:raise ValueError('input bundle changed')
        shutil.copyfile(src/(name+'.zip'),archive/(name+'.zip'))
    shutil.copyfile(args.preflight_bundle,archive/'preflight.zip');shutil.copyfile(src/'inputs.json',archive/'inputs.json');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'inputs.json')!=p['input_sha256']:raise ValueError('actual physical input changed')
    raw,rd=owned_zip(archive/'raw.zip','reference-study-003','reference',digest('validation/canonical-boat.json'),['Ainf','triangles_1'])
    triangles=data['reference_triangles_1'];np.testing.assert_array_equal(triangles,rd['triangles_1'])
    if triangles.shape!=(1232,3,3) or raw['input_sha256']!=p['input_sha256']:raise ValueError('coarsest physical source differs')
    mean=SimulationInputs.from_json((archive/'inputs.json').read_text());np.testing.assert_array_equal(mean.initial_state.position,data['world_cg'])
    package=distribution('capytaine')
    if version('capytaine')!='3.0.0':raise ValueError('pinned Capytaine required')
    dependencies={str(name):digest(Path(package.locate_file(name))) for name in package.files
        if ('Delhommeau_float64' in str(name) and str(name).endswith(('.pyd','.so'))) or str(name) in ('capytaine/bem/solver.py','capytaine/bem/engines.py','capytaine/green_functions/delhommeau.py')}
    if not any('Delhommeau_float64' in name for name in dependencies):raise ValueError('native library identity missing')
    host=dict(available_physical_bytes=available_memory(),available_disk_bytes=shutil.disk_usage(Path.cwd()).free)
    if host['available_physical_bytes']<1.25*c['maximum_peak_bytes'] or host['available_disk_bytes']<2*c['maximum_additional_disk_bytes']:raise ValueError('host margins fail')
    r.update(scope=__doc__,resource_only=False,synthetic_fields=False,input_sha256=p['input_sha256'],preflight_manifest_sha256=digest(args.preflight),
        preflight_bundle_sha256=digest(args.preflight_bundle),input_bundle_sha256=p['bundle_sha256'],dependency_sha256=dependencies,
        protocol_sha256=digest(c['protocol']),host_resources=host,completed_capytaine_calls=0,completed_orders=0,orders=[],
        source_export_verified=False,full_reference_qualified=False,source_discretization_qualified=False,full_run_authorized=False,
        original_raw_passed=False,original_native_gates=raw['gates'],memory_records=[],new_source_solves_required=6,
        dof_order=['Surge','Sway','Heave','Roll','Pitch','Yaw'],
        trace_normalization='potential/sources/bc divided by -i*encounter_omega before symbolic infinity cancellation; unit CG acceleration potentials',
        auxiliary_potential_units=['m','m','m','m^2','m^2','m^2'])
    a.update(source_triangles=triangles,world_cg=data['world_cg'],dry_mass=data['dry_mass'],original_native_Ainf=rd['Ainf'][0])
    def budget():
        memory=process_memory();r['latest_memory']=memory
        if max(memory['peak_working_set_bytes'],memory['peak_committed_bytes'])>c['maximum_peak_bytes']:raise MemoryError('pilot process memory ceiling exceeded')
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('pilot total wall budget exhausted')
    def stage(name):
        budget();r['memory_records'].append(dict(stage=name,**process_memory()));write_report(args.output,r)
    stage('inputs');checkpoint_arrays(args.output,r,a)
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        if any(row['num_threads']!=1 for row in threadpool_info()):raise ValueError('one native thread required')
        r['native_thread_pools']=threadpool_info();r['surface_proof']=verify_surface(data['wetted_triangles'],triangles)
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        if green.exportable_settings!=p['green_identity']['recorded']['infinite']:raise ValueError('infinite Green settings changed')
        r['green_settings']=green.exportable_settings
        scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:])
        def native_progress(arrays,completed):
            a.update(arrays);r['completed_capytaine_calls']=completed;checkpoint_arrays(args.output,r,a)
        def on_solved(completed):r.update(completed_capytaine_calls=completed);write_report(args.output,r)
        native_arrays,record=export_basis(cpt,green,triangles,data['world_cg'],mean.body.mass,rd['Ainf'][0],scale,budget,native_progress,on_solved)
        a.update(native_arrays);r.update(record);r['source_export_verified']=True
        sigma=native_arrays['normalized_sources'];A=native_arrays['native_Ainf']
        gc.collect();stage('source_export_verified');checkpoint_arrays(args.output,r,a)
        body=body_from_triangles(triangles,data['world_cg'],mean.body.mass)
        engine=cpt.BEMSolver(green_function=green).engine;seed=data['reference_seed_triangles'];panels=TrianglePanels(seed,-1)
        for order,count in zip(c['gauss_orders'],c['target_counts'],strict=True):
            key=f'target_{order}_';pts,owners,relative,error=mirrored_gauss_targets(seed,order)
            if len(pts)!=count:raise ValueError('required common targets missing')
            for name,value in [('points',pts),('owners',owners),('weights',panels.areas[owners]*relative),('normal',panels.normals[owners])]:np.testing.assert_array_equal(value,data[key+name])
            owned_points(pts,owners,seed,1e-10);budget()
        def fill(order,query,out):
            begin=perf_counter();evaluate_trace(engine,body.mesh,query,sigma,100,budget,out=out)
            if np.any(abs(out.imag)>1e-10*(1+abs(out.real))):raise ValueError('unexpected complex infinite trace')
            r['orders'].append(dict(order=order,points=len(query),seconds=perf_counter()-begin));r['completed_orders']+=1
        def checkpoint(order):
            checkpoint_arrays(args.output,r,a);stage(f'order{order}');print(f'fixed Ainf pilot: order{order} complete',flush=True)
        matrices=consume_orders(data,c['gauss_orders'],a,fill,budget,checkpoint)
        scaled_check(matrices.imag,np.zeros_like(matrices.real),scale,1e-10)
        real=matrices.real;U,gate=mesh_quality(real/scale,dict(reference_mesh_relative=.001,reference_mesh_absolute_normalized=1e-6))
        previous=abs(real[1]-real[0])/scale;latest=abs(real[2]-real[1])/scale;allow=1e-6+.001*abs(real[2]/scale)
        a.update(quadrature_A=real,quadrature_uncertainty=U*scale,previous_ratio=previous/allow,latest_ratio=latest/allow,
            trend_ratio=np.maximum(latest-previous,0)/1e-6,
            previous_failed=previous>allow,latest_failed=latest>allow,trend_failed=latest>previous+1e-6,
            output_change=real[-1]-A,output_change_absolute=abs(real[-1]-A))
        r['quadrature_stabilization']=gate;checkpoint_arrays(args.output,r,a)
        verify_arrays(args.output.with_suffix('.npz'),a,r['arrays_sha256']);stage('array_verification')
        r['payload_package']=package_payload(args.output,archive,r,budget);stage('package_verification')
    for name,h in dependencies.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('dependency changed during pilot')
    disk=sum(file.stat().st_size for file in archive.rglob('*') if file.is_file())+args.output.with_suffix('.npz').stat().st_size+args.output.with_suffix('.payload.zip').stat().st_size+65536
    if disk>c['maximum_additional_disk_bytes']:raise OSError('pilot disk ceiling exceeded')
    if r['completed_orders']!=3 or r['completed_capytaine_calls']!=6 or not r['source_export_verified']:raise ValueError('incomplete pilot cannot pass')
    budget();r.update(complete=True,passed=False,pilot_passed=bool(gate['passed']),accounted_disk_upper_bytes=disk,elapsed_seconds=perf_counter()-STARTED,
        interpretation='One source mesh only; conditional output integration, no source-mesh or finite-wave qualification')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/fixed-infinity-pilot-controls.json'))
    for name in ('preflight','preflight-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'fixed-seed-infinity-pilot');a={};checkpoint_arrays(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,pilot_passed=False,error=dict(type=type(error).__name__,message=str(error)))
    write_report(args.output,r);print(json.dumps({key:r.get(key) for key in ('complete','passed','pilot_passed','source_export_verified','completed_capytaine_calls','quadrature_stabilization','latest_memory','elapsed_seconds','error')}))
    return 0 if r['complete'] and r.get('pilot_passed') else 1


if __name__=='__main__':raise SystemExit(main())
