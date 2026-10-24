"""Bounded finest-source timing/memory probe, never a full field run or qualification."""
from time import perf_counter
STARTED=perf_counter()
import argparse,ctypes,json,shutil
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
from hashlib import sha256
from importlib.metadata import distribution,version
from importlib import import_module
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from .common import prepare,save,read_owned,digest
from .common_output_sources import pinned_bytes
from .common_output_preflight import available_memory
from .frozen_source_quadrature import evaluate_trace,owned_points,pressure_matrix
from .paired_body import scaled_check,capytaine_radius_check
from .reference_mesh import facets,verify_surface


def sample_starts(facet_ids,chunk=100,maximum=24):
    ids=np.asarray(facet_ids);n=len(ids);half=n//2
    if ids.ndim!=1 or n%2 or half<chunk or np.any(ids<0):raise ValueError('complete paired target grid required')
    starts=set()
    for base in (0,half):
        starts.update((base,base+(half-chunk)//2,base+half-chunk))
    for base in (0,half):
        for group in np.unique(ids[base:base+half]):
            covered=np.concatenate([np.arange(s,s+chunk) for s in sorted(starts) if base<=s<base+half])
            if np.any(ids[covered]==group):continue
            members=np.flatnonzero(ids[base:base+half]==group)+base
            starts.add(int(np.clip(members[len(members)//2]-chunk//2,base,base+half-chunk)))
    if len(starts)>maximum:raise ValueError('facet coverage exceeds declared probe blocks')
    return sorted(starts)


def projection(block_seconds,source_panels,field_pairs,native_seconds,safety=1.25,reserve=60.,chunk=100):
    values=np.asarray(block_seconds,float)
    if not len(values) or not np.isfinite(values).all() or np.any(values<=0):raise ValueError('complete finite positive sampled times required')
    rate=float(values.max()/(source_panels*chunk))
    return dict(maximum_sampled_seconds_per_pair=rate,projected_field_seconds=rate*field_pairs,
        projected_upper_seconds=(rate*field_pairs+native_seconds)*safety+reserve)


def process_memory():
    class Counters(ctypes.Structure):
        _fields_=[('cb',ctypes.c_uint32),('faults',ctypes.c_uint32)]+[(name,ctypes.c_size_t) for name in
            ('peak_working','working','peak_paged','paged','peak_nonpaged','nonpaged','pagefile','peak_pagefile','private')]
    kernel=ctypes.windll.kernel32;kernel.GetCurrentProcess.restype=ctypes.c_void_p
    api=ctypes.windll.psapi.GetProcessMemoryInfo;api.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_uint32];api.restype=ctypes.c_int
    counters=Counters();counters.cb=ctypes.sizeof(counters)
    if not api(kernel.GetCurrentProcess(),ctypes.byref(counters),counters.cb):raise OSError('process memory counters unavailable')
    return dict(peak_working_set_bytes=int(counters.peak_working),working_set_bytes=int(counters.working),
        peak_committed_bytes=int(counters.peak_pagefile),private_bytes=int(counters.private))


def owned_zip(bundle,stem,kind,control_hash):
    base='validation/'+stem
    with ZipFile(bundle) as z:
        r=json.loads(z.read(base+'.json'));payload=z.read(base+'.npz')
        if (r['kind']!=kind or r['complete'] is not True or r['contract_sha256']!=control_hash
                or sha256(payload).hexdigest()!=r['arrays_sha256']
                or sha256(z.read(base+'.source/contract.json')).hexdigest()!=control_hash):raise ValueError('corrupt or wrong-kind input graph')
        for name,h in r['producer_sha256'].items():
            if sha256(z.read(base+'.source/'+name)).hexdigest()!=h:raise ValueError('corrupt historical producer')
        with np.load(BytesIO(payload)) as a:arrays={name:a[name].copy() for name in a.files}
    return r,arrays


def generate(args,c,r,a):
    fixed=dict(format_version=1,maximum_elapsed_seconds=60.,source_panels=3456,gauss_orders=[12,16,20],cpu_threads=1,
        chunk_points=100,maximum_blocks_per_order=24,trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,
        maximum_peak_bytes=2147483648,future_maximum_elapsed_seconds=3600.,time_safety_factor=1.25,preparation_allowance_seconds=60.,
        preflight_bundle_sha256='11ea95f6e3c3a7a0249521b82f734bde21cb85cb679898f9e97c6f5d80dbcb3f',full_reference_scope=False)
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('undeclared native probe controls')
    if digest(args.preflight_bundle)!=c['preflight_bundle_sha256']:raise ValueError('preflight bundle differs')
    p,pa=read_owned(args.preflight,digest('validation/common-output-preflight-controls.json'),'common-output-finest-preflight')
    pinned_bytes(args.preflight_bundle,args.preflight.stem,args.preflight)
    if not p['resource_preflight_accepted'] or p['full_run_authorized'] or p['native_evaluations']!=0:raise ValueError('preflight scope differs')
    src=args.preflight.with_suffix('.source');archive=args.output.with_suffix('.source')
    for name,h in p['input_bundle_sha256'].items():
        if digest(src/(name+'.zip'))!=h:raise ValueError('actual input bundle differs')
        shutil.copyfile(src/(name+'.zip'),archive/(name+'.zip'))
    shutil.copyfile(args.preflight_bundle,archive/'preflight.zip');shutil.copyfile(src/'inputs.json',archive/'inputs.json');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'inputs.json')!=p['input_sha256']:raise ValueError('physical input differs')
    common,data=owned_zip(archive/'common.zip','common-output-sources-001','common-output-source-diagnostic',digest('validation/common-output-source-controls.json'))
    study,source=owned_zip(archive/'study.zip','paired-body-study-002','paired-body-study',digest('validation/paired-body-controls.json'))
    if common['input_sha256']!=p['input_sha256'] or study['input_sha256']!=p['input_sha256']:raise ValueError('source/target physical identity differs')
    for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(data[key],source[key]);np.testing.assert_array_equal(pa[key],source[key])
    triangles=source['triangles_3'];sigma=source['capytaine_sources_3']
    np.testing.assert_array_equal(triangles,pa['source_triangles']);np.testing.assert_array_equal(sigma,pa['normalized_sources'])
    mean=SimulationInputs.from_json((archive/'inputs.json').read_text());np.testing.assert_array_equal(mean.initial_state.position,source['world_cg'])
    package=distribution('capytaine');dependencies={name:digest(Path(package.locate_file(name))) for name in p['inspected_dependency_sha256']}
    if dependencies!=p['inspected_dependency_sha256'] or version('capytaine')!='3.0.0':raise ValueError('inspected allocation dependency changed')
    binaries={str(name):digest(Path(package.locate_file(name))) for name in package.files if 'Delhommeau_float64' in str(name) and str(name).endswith(('.pyd','.so'))}
    if not binaries:raise ValueError('actual native library identity missing')
    host=dict(available_physical_bytes=available_memory(),available_disk_bytes=shutil.disk_usage(Path.cwd()).free)
    model=p['resource_model']
    if host['available_physical_bytes']<1.25*model['estimated_peak_bytes'] or host['available_disk_bytes']<2*model['estimated_additional_disk_bytes']:raise ValueError('refreshed resource margins fail')
    r.update(scope=__doc__,input_sha256=p['input_sha256'],preflight_manifest_sha256=digest(args.preflight),
        preflight_bundle_sha256=digest(args.preflight_bundle),input_bundle_sha256=p['input_bundle_sha256'],
        inspected_dependency_sha256=dependencies,native_binary_sha256=binaries,host_resources=host,
        new_bem_solves=0,native_evaluated_points=0,full_run_authorized=False,full_reference_qualified=False,
        source_discretization_qualified=False,completed_orders=0,orders=[],memory_records=[],protocol_sha256=digest(c['protocol']))
    a.update(source_triangles=triangles,normalized_sources=sigma,world_cg=source['world_cg'],dry_mass=source['dry_mass'])
    def budget():
        memory=process_memory();r['latest_memory']=memory
        if max(memory['peak_working_set_bytes'],memory['peak_committed_bytes'])>c['maximum_peak_bytes']:raise MemoryError('native probe process memory ceiling exceeded')
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('native probe total budget exhausted')
    budget();r['memory_records'].append(dict(stage='inputs',**process_memory()));save(args.output,r,a)
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        if any(x['num_threads']!=1 for x in threadpool_info()):raise ValueError('one native thread required')
        r['native_thread_pools']=threadpool_info();r['surface_proof']=verify_surface(source['parent_body_triangles'],triangles)
        body=body_from_triangles(triangles,source['world_cg'],mean.body.mass);r['radius_error_m']=capytaine_radius_check(body,triangles)
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        if green.exportable_settings!=study['green_settings'] or green.exportable_settings!=common['green_settings']:raise ValueError('Green settings changed')
        r['green_settings']=green.exportable_settings;engine=cpt.BEMSolver(green_function=green).engine
        scale=np.sqrt(source['dry_mass'].diagonal()[:,None]*source['dry_mass'].diagonal()[None,:])
        r['preparation_seconds']=perf_counter()-STARTED;begin=perf_counter()
        phi=evaluate_trace(engine,body.mesh,body.mesh.faces_centers,sigma,c['chunk_points'],budget)
        native_A=pressure_matrix(body.mesh.faces_centers,-body.mesh.faces_normals,body.mesh.faces_areas,phi,source['world_cg'],1025.)
        r['native_seconds']=perf_counter()-begin;r['native_evaluated_points']+=len(phi)
        r['native_potential_ratio']=scaled_check(phi,source['capytaine_aux_3'],1.,c['trace_scaled_tolerance'])
        r['native_load_ratio']=scaled_check(native_A,source['capytaine_Ainf'][3],scale,c['trace_scaled_tolerance'])
        a.update(native_points=body.mesh.faces_centers,native_potential=phi,native_A_complex=native_A)
        r['memory_records'].append(dict(stage='native',**process_memory()));save(args.output,r,a);budget()
        physical_normals=np.array([f['normal'] for f in facets(source['parent_body_triangles'])]);times=[]
        for order in c['gauss_orders']:
            key=f'target_{order}_';points=data[key+'points'];normal=data[key+'normal']
            distance=np.linalg.norm(normal[:,None,:]+physical_normals[None,:,:],axis=2)
            labels=distance.argmin(axis=1)
            if np.max(distance.min(axis=1))>1e-8:raise ValueError('common output physical normals differ')
            starts=sample_starts(labels,c['chunk_points'],c['maximum_blocks_per_order']);record=dict(order=order,starts=starts,block_seconds=[])
            indices=np.array([np.arange(s,s+c['chunk_points']) for s in starts]);selected=indices.ravel()
            prefix=f'order_{order}_'
            proof=owned_points(points[selected],data[key+'owners'][selected],data['target_triangles'],c['ownership_absolute_m'])
            a.update({prefix+'indices':indices,prefix+'points':points[selected],prefix+'owners':data[key+'owners'][selected],
                prefix+'weights':data[key+'weights'][selected],prefix+'facet_ids':labels[selected]})
            fields=[]
            for indices_block in indices:
                begin=perf_counter();value=evaluate_trace(engine,body.mesh,points[indices_block],sigma,c['chunk_points'],budget)
                seconds=perf_counter()-begin;times.append(seconds);record['block_seconds'].append(seconds);fields.append(value)
                r['native_evaluated_points']+=len(value)
                if np.any(abs(value.imag)>1e-10*(1+abs(value.real))):raise ValueError('unexpected complex infinity trace')
            a[prefix+'potential']=np.concatenate(fields);record.update(ownership=proof,facet_ids=np.unique(labels[selected]).tolist())
            r['orders'].append(record);r['completed_orders']+=1;r['memory_records'].append(dict(stage=f'order{order}',**process_memory()))
            save(args.output,r,a);budget();print(f'finest cost probe: order{order}, {len(starts)} blocks',flush=True)
    estimate=projection(times,3456,p['work']['field_pairs'],r['native_seconds'],c['time_safety_factor'],max(c['preparation_allowance_seconds'],r['preparation_seconds']))
    for name,h in {**dependencies,**binaries}.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('dependency changed during probe')
    budget();r.update(complete=True,passed=False,projection=estimate,
        probe_feasible=estimate['projected_upper_seconds']<=c['future_maximum_elapsed_seconds'],elapsed_seconds=perf_counter()-STARTED)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/finest-cost-probe-controls.json'))
    for name in ('preflight','preflight-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'finest-source-cost-probe');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a);print(json.dumps({key:r.get(key) for key in ('complete','passed','probe_feasible','native_evaluated_points','projection','latest_memory','elapsed_seconds','error')}))
    return 0 if r['complete'] and r.get('probe_feasible') else 1


if __name__=='__main__':raise SystemExit(main())
