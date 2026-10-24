"""Shared finest-source buffer/load/checkpoint consumer; resource rehearsal only for now."""
from time import perf_counter
STARTED=perf_counter()
import argparse,json,shutil
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
from importlib import import_module
from importlib.metadata import distribution
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from .common import prepare,read_owned,digest
from .common_output_sources import pinned_bytes
from .common_output_preflight import available_memory
from .finest_cost_probe import owned_zip,process_memory
from .frozen_source_quadrature import evaluate_trace,owned_points,pressure_matrix
from .paired_body import scaled_check,capytaine_radius_check
from .reference_mesh import verify_surface


def stream_hash(file):
    h=sha256()
    while data:=file.read(1024*1024):h.update(data)
    return h.hexdigest()


def write_report(output,r):output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')


def checkpoint_arrays(output,r,arrays):
    path=output.with_suffix('.npz');np.savez(path,**arrays)
    with path.open('rb') as file:r['arrays_sha256']=stream_hash(file)
    write_report(output,r)


def verify_arrays(path,arrays,expected_hash):
    with path.open('rb') as file:
        if stream_hash(file)!=expected_hash:raise ValueError('checkpoint array bytes changed')
    with np.load(path) as restored:
        if set(restored.files)!=set(arrays):raise ValueError('checkpoint fields missing')
        for name,value in arrays.items():np.testing.assert_array_equal(restored[name],value)


def consume_orders(data,orders,arrays,fill,budget,checkpoint):
    """Same retained buffers and full load work for actual or rehearsal providers."""
    matrices=[]
    for order in orders:
        budget();key=f'target_{order}_';prefix=f'order_{order}_'
        points=data[key+'points'];normal=data[key+'normal'];weights=data[key+'weights']
        for name in ('points','owners','weights','normal'):arrays[key+name]=data[key+name]
        field=np.empty((len(points),6),complex);arrays[prefix+'potential']=field
        fill(order,points,field)
        if not np.isfinite(field).all():raise ValueError('incomplete/nonfinite retained output')
        matrix=pressure_matrix(points,normal,weights,field,arrays['world_cg'],1025.)
        arrays[prefix+'A_complex']=matrix;matrices.append(matrix);budget();checkpoint(order)
    return np.stack(matrices)


def synthetic_fill(order,points,out,budget):
    # High-entropy synthetic complex data touches every byte and stresses packaging.
    rng=np.random.default_rng(24000+order)
    for start in range(0,len(points),100):
        end=min(start+100,len(points));out[start:end].real=rng.standard_normal((end-start,6));out[start:end].imag=rng.standard_normal((end-start,6))
    budget()


def package_payload(output,archive,r,budget):
    """Exercise the real archive writer and streaming byte verification."""
    path=output.with_suffix('.payload.zip');snapshot=output.read_bytes();npz=output.with_suffix('.npz')
    with ZipFile(path,'x',ZIP_DEFLATED) as z:
        z.write(npz,'arrays.npz');budget();z.writestr('payload-manifest.json',snapshot)
        for file in sorted(archive.rglob('*')):
            if file.is_file() and file.name not in {name+'.zip' for name in r['input_bundle_sha256']}:z.write(file,'source/'+file.relative_to(archive).as_posix())
        z.writestr('external-inputs.json',json.dumps(r['input_bundle_sha256'],indent=2))
    budget()
    with ZipFile(path) as z:
        if z.testzip() is not None:raise ValueError('packaged CRC failed')
        with z.open('arrays.npz') as file:
            if stream_hash(file)!=r['arrays_sha256']:raise ValueError('packaged arrays changed')
        if z.read('payload-manifest.json')!=snapshot:raise ValueError('packaged metadata changed')
    with path.open('rb') as file:h=stream_hash(file)
    budget();return dict(path=path.name,bytes=path.stat().st_size,sha256=h,snapshot_sha256=sha256(snapshot).hexdigest())


def generate(args,c,r,a):
    fixed=dict(format_version=1,mode='resource-rehearsal',maximum_elapsed_seconds=60.,maximum_peak_bytes=2147483648,
        maximum_additional_disk_bytes=2147483648,source_panels=3456,gauss_orders=[12,16,20],target_counts=[124416,221184,345600],
        chunk_points=100,cpu_threads=1,trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,
        probe_bundle_sha256='cd5756817dcd38a99ae4931fa9f3befe5c8ba733a4cb1ace74cbe9852bb084cc',full_reference_scope=False)
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('no full native run authorized by these controls')
    if digest(args.probe_bundle)!=c['probe_bundle_sha256']:raise ValueError('probe bundle differs')
    probe,probe_arrays=read_owned(args.probe,digest('validation/finest-cost-probe-controls.json'),'finest-source-cost-probe')
    pinned_bytes(args.probe_bundle,args.probe.stem,args.probe)
    if not probe['probe_feasible'] or probe['full_run_authorized']:raise ValueError('limited probe scope differs')
    src=args.probe.with_suffix('.source');archive=args.output.with_suffix('.source')
    for name,h in probe['input_bundle_sha256'].items():
        if digest(src/(name+'.zip'))!=h:raise ValueError('actual input bundle differs')
        shutil.copyfile(src/(name+'.zip'),archive/(name+'.zip'))
    shutil.copyfile(args.probe_bundle,archive/'probe.zip');shutil.copyfile(src/'preflight.zip',archive/'preflight.zip')
    if digest(archive/'preflight.zip')!=probe['preflight_bundle_sha256']:raise ValueError('embedded preflight identity differs')
    shutil.copyfile(src/'inputs.json',archive/'inputs.json');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'inputs.json')!=probe['input_sha256']:raise ValueError('physical input differs')
    target_keys=['world_cg','dry_mass','target_triangles','source_A','quadrature_uncertainty']+[f'target_{order}_{name}' for order in c['gauss_orders'] for name in ('points','owners','weights','normal')]
    common,data=owned_zip(archive/'common.zip','common-output-sources-001','common-output-source-diagnostic',digest('validation/common-output-source-controls.json'),target_keys)
    source_keys=['world_cg','dry_mass','parent_body_triangles','triangles_3','capytaine_sources_3','capytaine_aux_3','capytaine_Ainf']
    study,source=owned_zip(archive/'study.zip','paired-body-study-002','paired-body-study',digest('validation/paired-body-controls.json'),source_keys)
    if common['input_sha256']!=probe['input_sha256'] or study['input_sha256']!=probe['input_sha256']:raise ValueError('physical identity differs')
    for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(data[key],source[key]);np.testing.assert_array_equal(probe_arrays[key],source[key])
    for key,old in (('triangles_3','source_triangles'),('capytaine_sources_3','normalized_sources')):np.testing.assert_array_equal(source[key],probe_arrays[old])
    mean=SimulationInputs.from_json((archive/'inputs.json').read_text());np.testing.assert_array_equal(mean.initial_state.position,source['world_cg'])
    package=distribution('capytaine')
    identities={**probe['inspected_dependency_sha256'],**probe['native_binary_sha256']}
    for name,h in identities.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('dependency identity changed')
    with ZipFile(src/'preflight.zip') as z:p=json.loads(z.read('validation/common-output-finest-preflight-001.json'))
    host=dict(available_physical_bytes=available_memory(),available_disk_bytes=shutil.disk_usage(Path.cwd()).free)
    if host['available_physical_bytes']<1.25*p['resource_model']['estimated_peak_bytes'] or host['available_disk_bytes']<2*c['maximum_additional_disk_bytes']:raise ValueError('host margins fail')
    r.update(mode=c['mode'],resource_only=True,synthetic_fields=True,rehearsal_passed=False,full_run_authorized=False,
        full_reference_qualified=False,source_discretization_qualified=False,new_bem_solves=0,native_evaluated_points=0,
        input_sha256=probe['input_sha256'],probe_manifest_sha256=digest(args.probe),probe_bundle_sha256=digest(args.probe_bundle),
        preflight_bundle_sha256=probe['preflight_bundle_sha256'],
        input_bundle_sha256=probe['input_bundle_sha256'],dependency_sha256=identities,host_resources=host,
        protocol_sha256=digest(c['protocol']),completed_stages=[],memory_records=[],scope=__doc__)
    def budget():
        memory=process_memory();r['latest_memory']=memory
        if max(memory['peak_working_set_bytes'],memory['peak_committed_bytes'])>c['maximum_peak_bytes']:raise MemoryError('consumer process memory cap exceeded')
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('consumer total budget exhausted')
    def stage(name):
        budget();r['completed_stages'].append(name);r['memory_records'].append(dict(stage=name,**process_memory()));write_report(args.output,r)
    stage('inputs')
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        if any(row['num_threads']!=1 for row in threadpool_info()):raise ValueError('one native thread required')
        r['native_thread_pools']=threadpool_info();r['surface_proof']=verify_surface(source['parent_body_triangles'],source['triangles_3'])
        body=body_from_triangles(source['triangles_3'],source['world_cg'],mean.body.mass);r['radius_error_m']=capytaine_radius_check(body,source['triangles_3'])
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        if green.exportable_settings!=study['green_settings'] or green.exportable_settings!=common['green_settings']:raise ValueError('Green settings differ')
        r['green_settings']=green.exportable_settings;engine=cpt.BEMSolver(green_function=green).engine
        a.update(world_cg=source['world_cg'],dry_mass=source['dry_mass'],source_triangles=source['triangles_3'],normalized_sources=source['capytaine_sources_3'])
        native=evaluate_trace(engine,body.mesh,body.mesh.faces_centers,source['capytaine_sources_3'],100,budget)
        native_A=pressure_matrix(body.mesh.faces_centers,-body.mesh.faces_normals,body.mesh.faces_areas,native,source['world_cg'],1025.)
        scale=np.sqrt(source['dry_mass'].diagonal()[:,None]*source['dry_mass'].diagonal()[None,:])
        r['native_potential_ratio']=scaled_check(native,source['capytaine_aux_3'],1.,1e-10);r['native_load_ratio']=scaled_check(native_A,source['capytaine_Ainf'][3],scale,1e-10)
        r['native_evaluated_points']=len(native);a.update(native_points=body.mesh.faces_centers,native_potential=native,native_A_complex=native_A);stage('native')
        for order,count in zip(c['gauss_orders'],c['target_counts'],strict=True):
            key=f'target_{order}_'
            if data[key+'points'].shape!=(count,3):raise ValueError('full common target grid missing')
            owned_points(data[key+'points'],data[key+'owners'],data['target_triangles'],c['ownership_absolute_m']);budget()
        def checkpoint(order):
            checkpoint_arrays(args.output,r,a);stage(f'checkpoint{order}')
        matrices=consume_orders(data,c['gauss_orders'],a,lambda order,points,out:synthetic_fill(order,points,out,budget),budget,checkpoint)
        a['resource_load_matrices']=matrices
        checkpoint_arrays(args.output,r,a);stage('final_checkpoint')
        verify_arrays(args.output.with_suffix('.npz'),a,r['arrays_sha256']);stage('array_verification')
        r['payload_package']=package_payload(args.output,archive,r,budget);stage('package_verification')
    for name,h in identities.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('dependency changed during consumer')
    required=['inputs','native','checkpoint12','checkpoint16','checkpoint20','final_checkpoint','array_verification','package_verification']
    if r['completed_stages']!=required:raise ValueError('required consumer stage missing')
    disk=sum(file.stat().st_size for file in archive.rglob('*') if file.is_file())+args.output.with_suffix('.npz').stat().st_size+args.output.with_suffix('.payload.zip').stat().st_size+65536
    if disk>c['maximum_additional_disk_bytes']:raise OSError('consumer actual additional disk cap exceeded')
    budget();r.update(complete=True,passed=False,rehearsal_passed=True,accounted_additional_disk_upper_bytes=disk,
        retained_field_bytes=sum(a[f'order_{order}_potential'].nbytes for order in c['gauss_orders']),elapsed_seconds=perf_counter()-STARTED)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/finest-consumer-rehearsal-controls.json'))
    for name in ('probe','probe-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'finest-consumer-resource-rehearsal');a={};checkpoint_arrays(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,rehearsal_passed=False,error=dict(type=type(error).__name__,message=str(error)))
    write_report(args.output,r);print(json.dumps({key:r.get(key) for key in ('complete','passed','rehearsal_passed','retained_field_bytes','latest_memory','elapsed_seconds','error')}))
    return 0 if r['complete'] and r.get('rehearsal_passed') else 1


if __name__=='__main__':raise SystemExit(main())
