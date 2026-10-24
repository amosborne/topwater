"""Read-only resource model for the complete 3456-source common-grid continuation."""
from time import perf_counter
STARTED=perf_counter()
import argparse,ctypes,json,shutil
from importlib.metadata import distribution,version
from pathlib import Path
from zipfile import ZipFile
import numpy as np
from .common import prepare,save,read_owned,digest
from .common_output_sources import pinned_bytes


def work_model(sources,targets,chunk=100,modes=6):
    if sources<1 or chunk<1 or modes!=6 or len(targets)!=3 or any(n<1 for n in targets):raise ValueError('complete positive three-order six-mode work required')
    total=sum(targets);native=sources
    return dict(source_panels=sources,target_counts=list(targets),target_points=total,
        field_pairs=sources*total,native_calibration_pairs=sources*native,
        all_pairs=sources*(total+native),new_field_bytes=(total+native)*modes*16,
        new_native_point_bytes=native*3*8,
        # Installed Capytaine default early-dot path allocates complex S AND K.
        native_matrix_bytes=2*min(chunk,max(targets))*sources*16,
        load_work_bytes=max(targets)*modes*(8+16),
        block_calls=sum((n+chunk-1)//chunk for n in targets)+(native+chunk-1)//chunk)


def resource_model(work,rates,resident,largest_payload,input_bundles,producer_bytes,c):
    if not rates or min(rates)<=0 or not np.isfinite(rates).all():raise ValueError('finite measured positive pair rates required')
    field_seconds=[work['all_pairs']*rate for rate in rates]
    upper=max(field_seconds)*c['time_safety_factor']+c['preparation_allowance_seconds']
    output=work['new_field_bytes']+work['new_native_point_bytes']+resident
    memory=(resident+2*largest_payload+work['new_field_bytes']+work['new_native_point_bytes']+
        c['native_scratch_multiplier']*work['native_matrix_bytes']+work['load_work_bytes']+c['runtime_overhead_bytes'])
    disk=2*input_bundles+4*output+2*producer_bytes
    return dict(measured_rate_seconds_per_pair=list(rates),calibrated_native_seconds=field_seconds,
        estimated_upper_wall_seconds=upper,estimated_peak_bytes=memory,estimated_raw_output_upper_bytes=output,
        estimated_additional_disk_bytes=disk,
        ceiling_checks=dict(time=upper<=c['future_maximum_elapsed_seconds'],memory=memory<=c['maximum_peak_bytes'],disk=disk<=c['maximum_additional_disk_bytes']))


def available_memory():
    class Status(ctypes.Structure):
        _fields_=[('length',ctypes.c_uint32),('load',ctypes.c_uint32)]+[(name,ctypes.c_uint64) for name in ('total','available','page_total','page_available','virtual_total','virtual_available','extended')]
    value=Status();value.length=ctypes.sizeof(value)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(value)):raise OSError('available physical memory unavailable')
    return int(value.available)


def generate(args,c,r,a):
    fixed=dict(format_version=1,source_level=3,source_panels=3456,target_counts=[124416,221184,345600],gauss_orders=[12,16,20],
        evaluation_chunk_points=100,modes=6,maximum_elapsed_seconds=60.,future_maximum_elapsed_seconds=3600.,
        maximum_peak_bytes=2147483648,maximum_additional_disk_bytes=2147483648,time_safety_factor=1.25,
        preparation_allowance_seconds=60.,native_scratch_multiplier=4,runtime_overhead_bytes=536870912,
        available_memory_margin=1.25,available_disk_margin=2.,full_reference_scope=False,
        common_bundle_sha256='c208dffe223dffb1db87af8b8bf0eabd96d8ad6a40b358e924c26208dd7fbef1',
        study_bundle_sha256='052df7c50865340e7ef4ebaddfc48adfef191b24057646bd395a889515dc98f2',
        quadrature_bundle_sha256='7e999761722f68b52c9675abfb32bc02aff101c55f03ad76babf749154ab7a63',
        historical_bundle_sha256='28a35794d8c20d35e59f955e0b4ced8fec3cdebd4c6f8e01beffd1cb205fdbb0')
    if any(c.get(key)!=value for key,value in fixed.items()):raise ValueError('undeclared preflight controls')
    for name in ('common','study','quadrature','historical'):
        if digest(getattr(args,name+'_bundle'))!=c[name+'_bundle_sha256']:raise ValueError('input bundle identity differs')
    common,data=read_owned(args.common,digest('validation/common-output-source-controls.json'),'common-output-source-diagnostic')
    study,source=read_owned(args.study,digest('validation/paired-body-controls.json'),'paired-body-study')
    quad,qdata=read_owned(args.quadrature,digest('validation/frozen-source-gauss-continuation-controls.json'),'frozen-source-pressure-quadrature')
    for name in ('common','study','quadrature'):pinned_bytes(getattr(args,name+'_bundle'),getattr(args,name).stem,getattr(args,name))
    if (common['completed_sources']!=3 or not common['source_diagnostic']['all_quadrature_passed']
            or common['study_manifest_sha256']!=digest(args.study) or common['study_arrays_sha256']!=digest(args.study.with_suffix('.npz'))
            or common['quadrature_manifest_sha256']!=digest(args.quadrature) or common['quadrature_arrays_sha256']!=digest(args.quadrature.with_suffix('.npz'))
            or common['input_sha256']!=study['input_sha256'] or quad['input_sha256']!=study['input_sha256']
            or any(x['full_reference_qualified'] is not False for x in (common,study,quad))):raise ValueError('prior scope/physical/provenance mismatch')
    for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(data[key],source[key]);np.testing.assert_array_equal(qdata[key],source[key])
    triangles=source['triangles_3'];sigma=source['capytaine_sources_3']
    if triangles.shape!=(3456,3,3) or sigma.shape!=(3456,6):raise ValueError('missing complete finest source solution')
    for i,(order,n) in enumerate(zip(c['gauss_orders'],c['target_counts'],strict=True),1):
        prefix=f'target_{order}_'
        if data[prefix+'points'].shape!=(n,3):raise ValueError('missing common target grid')
        for key in ('points','owners','weights'):np.testing.assert_array_equal(data[prefix+key],qdata[f'level{i}_'+key])
    archive=args.output.with_suffix('.source')
    for name in ('common','study','quadrature','historical'):shutil.copyfile(getattr(args,name+'_bundle'),archive/(name+'.zip'))
    shutil.copyfile(c['protocol'],archive/'protocol.md');physical=args.study.with_suffix('.source')/'inputs.json';shutil.copyfile(physical,archive/'inputs.json')
    if digest(archive/'inputs.json')!=study['input_sha256']:raise ValueError('actual physical input differs')
    # Record the inspected external allocation path, without importing/executing Capytaine.
    package=distribution('capytaine');dependencies={}
    if version('capytaine')!='3.0.0':raise ValueError('pinned external allocation model required')
    for name in ('bem/engines.py','green_functions/delhommeau.py','green_functions/abstract_green_function.py'):
        file=Path(package.locate_file('capytaine/'+name));dependencies['capytaine/'+name]=digest(file)
    r.update(input_sha256=study['input_sha256'],input_bundle_sha256={name:c[name+'_bundle_sha256'] for name in ('common','study','quadrature','historical')},
        common_manifest_sha256=digest(args.common),study_manifest_sha256=digest(args.study),quadrature_manifest_sha256=digest(args.quadrature),
        protocol_sha256=digest(c['protocol']),inspected_dependency_sha256=dependencies,capytaine_version=version('capytaine'),
        native_evaluations=0,new_bem_solves=0,full_reference_qualified=False,source_discretization_qualified=False,scope=__doc__)
    a.update(source_triangles=triangles,normalized_sources=sigma,world_cg=source['world_cg'],dry_mass=source['dry_mass'])
    work=work_model(3456,c['target_counts'],c['evaluation_chunk_points'])
    observations=[]
    for row in common['source_records'][:2]:
        for measured,n in zip(row['orders'],c['target_counts'],strict=True):
            observations.append(dict(source_panels=row['panels'],targets=n,seconds=measured['seconds'],origin='common-output'))
    observations.append(dict(source_panels=864,targets=345600,seconds=quad['levels'][-1]['seconds'],origin='order20-continuation'))
    rates=[x['seconds']/(x['source_panels']*x['targets']) for x in observations]
    resident=sum(x.nbytes for arrays in (data,source,qdata) for x in arrays.values())
    payload=max(getattr(args,name).with_suffix('.npz').stat().st_size for name in ('common','study','quadrature'))
    bundle_bytes=sum(getattr(args,name+'_bundle').stat().st_size for name in ('common','study','quadrature','historical'))
    producer_bytes=sum((archive/name).stat().st_size for name in r['producer_sha256'])
    model=resource_model(work,rates,resident,payload,bundle_bytes,producer_bytes,c)
    free_memory=available_memory();free_disk=shutil.disk_usage(Path.cwd()).free
    host=dict(available_physical_bytes=free_memory,available_disk_bytes=free_disk,
        memory_available=free_memory>=model['estimated_peak_bytes']*c['available_memory_margin'],
        disk_available=free_disk>=model['estimated_additional_disk_bytes']*c['available_disk_margin'])
    accepted=all(model['ceiling_checks'].values()) and host['memory_available'] and host['disk_available']
    r.update(complete=True,passed=False,work=work,timing_observations=observations,resource_model=model,host_resources=host,
        resource_preflight_accepted=accepted,full_run_authorized=False,
        decision='Proceed only to a separately bounded representative finest-source cost screen' if accepted else 'Reject full continuation under this resource model; resolve failed resource gates',
        elapsed_seconds=perf_counter()-STARTED)
    if r['elapsed_seconds']>c['maximum_elapsed_seconds']:raise TimeoutError('read-only preflight budget exhausted')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/common-output-preflight-controls.json'))
    for name in ('common','common-bundle','study','study-bundle','quadrature','quadrature-bundle','historical-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'common-output-finest-preflight');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except Exception as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a);print(json.dumps({key:r.get(key) for key in ('complete','passed','resource_preflight_accepted','resource_model','decision','elapsed_seconds','error')}))
    return 0 if r['complete'] and r.get('resource_preflight_accepted') else 1


if __name__=='__main__':raise SystemExit(main())
