"""Remaining fixed-seed Ainf source exports only; no common-grid field evaluations."""
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
from .common import prepare,read_owned,digest,current_sources,PHYSICAL
from .common_output_sources import pinned_bytes
from .common_output_preflight import available_memory
from .finest_cost_probe import owned_zip,process_memory
from .finest_consumer import checkpoint_arrays,write_report,verify_arrays,package_payload
from .reference_mesh import verify_surface
from .ainf_sources import export_basis


def dense_plan(n):
    matrix=n*n*16
    # Two cache matrices, two pristine copies, LU, one full-size transient.
    return dict(panels=n,one_complex_matrix_bytes=matrix,dense_matrices=6,dense_bytes=6*matrix,
        baseline_allowance_bytes=2*1024**3,workspace_allowance_bytes=256*1024**2,basis_allowance_bytes=8*n*6*16,
        estimated_peak_bytes=6*matrix+2*1024**3+256*1024**2+8*n*6*16)


def generate(args,c,r,a):
    fixed=dict(format_version=1,levels=[2,3,4],counts=[2772,4928,7700],new_source_calls=18,maximum_elapsed_seconds=600.,
        maximum_peak_bytes=8589934592,maximum_additional_disk_bytes=2147483648,baseline_allowance_bytes=2147483648,cpu_threads=1,
        trace_scaled_tolerance=1e-10,pilot_bundle_sha256='e4d43719636db7ac67731b5f98d418b247371b585f2bd89069b3405eaaee2e79',full_reference_scope=False)
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('undeclared source-only export controls')
    if digest(args.pilot_bundle)!=c['pilot_bundle_sha256']:raise ValueError('pilot bundle differs')
    pilot,pd=read_owned(args.pilot,digest('validation/fixed-infinity-pilot-controls.json'),'fixed-seed-infinity-pilot')
    pinned_bytes(args.pilot_bundle,args.pilot.stem,args.pilot);current_sources(pilot,PHYSICAL)
    if not pilot['source_export_verified'] or not pilot['pilot_passed'] or pilot['completed_capytaine_calls']!=6 or pilot['synthetic_fields']:raise ValueError('actual coarsest source proof missing')
    src=args.pilot.with_suffix('.source');archive=args.output.with_suffix('.source')
    for name,h in pilot['input_bundle_sha256'].items():
        if digest(src/(name+'.zip'))!=h:raise ValueError('input changed')
        shutil.copyfile(src/(name+'.zip'),archive/(name+'.zip'))
    shutil.copyfile(args.pilot_bundle,archive/'pilot.zip');shutil.copyfile(src/'preflight.zip',archive/'preflight.zip');shutil.copyfile(src/'inputs.json',archive/'inputs.json');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'inputs.json')!=pilot['input_sha256'] or digest(archive/'preflight.zip')!=pilot['preflight_bundle_sha256']:raise ValueError('physical/asset graph changed')
    keys=['wetted_triangles','world_cg','dry_mass','body_mass']+[f'reference_triangles_{level}' for level in [1,2,3,4]]
    frozen,fd=owned_zip(archive/'frozen.zip','frozen-fixed-001','freeze',digest('validation/canonical-boat.json'),keys)
    raw,rd=owned_zip(archive/'raw.zip','reference-study-003','reference',digest('validation/canonical-boat.json'),['Ainf'])
    mean=SimulationInputs.from_json((archive/'inputs.json').read_text())
    for key in ('world_cg','dry_mass'):np.testing.assert_array_equal(pd[key],fd[key])
    np.testing.assert_array_equal(pd['source_triangles'],fd['reference_triangles_1']);np.testing.assert_array_equal(pd['native_Ainf'],rd['Ainf'][0])
    np.testing.assert_array_equal(mean.initial_state.position,fd['world_cg'])
    package=distribution('capytaine')
    if version('capytaine')!='3.0.0':raise ValueError('pinned allocation path required')
    inspected=['capytaine/bem/engines.py','capytaine/bem/solver.py','capytaine/bodies/bodies.py','capytaine/tools/block_circulant_matrices.py','capytaine/green_functions/abstract_green_function.py','capytaine/green_functions/delhommeau.py']
    identities={name:digest(Path(package.locate_file(name))) for name in inspected}
    identities.update({str(name):digest(Path(package.locate_file(name))) for name in package.files if 'Delhommeau_float64' in str(name) and str(name).endswith(('.pyd','.so'))})
    plans=[dense_plan(n) for n in c['counts']]
    if any(plan['estimated_peak_bytes']>c['maximum_peak_bytes'] for plan in plans):raise MemoryError('explicit dense model exceeds native cap')
    host=dict(available_physical_bytes=available_memory(),available_disk_bytes=shutil.disk_usage(Path.cwd()).free)
    if host['available_physical_bytes']<1.25*max(x['estimated_peak_bytes'] for x in plans) or host['available_disk_bytes']<2*c['maximum_additional_disk_bytes']:raise ValueError('prelaunch host margins fail')
    r.update(scope=__doc__,resource_only=False,synthetic_fields=False,source_export_complete=False,source_discretization_qualified=False,
        full_reference_qualified=False,full_run_authorized=False,common_grid_evaluations=0,completed_capytaine_calls=0,levels=[],memory_records=[],
        allocation_model=plans,host_resources=host,dependency_sha256=identities,input_sha256=pilot['input_sha256'],input_bundle_sha256=pilot['input_bundle_sha256'],
        pilot_manifest_sha256=digest(args.pilot),pilot_bundle_sha256=digest(args.pilot_bundle),preflight_bundle_sha256=pilot['preflight_bundle_sha256'],
        protocol_sha256=digest(c['protocol']),dof_order=pilot['dof_order'],trace_normalization=pilot['trace_normalization'],auxiliary_potential_units=pilot['auxiliary_potential_units'],
        reused_coarsest=dict(current_commit=pilot['current_commit'],cases=6,source_checks=pilot['source_checks']),covered_source_levels=[1])
    a.update(world_cg=fd['world_cg'],dry_mass=fd['dry_mass'])
    basis_keys=['source_triangles','native_points','native_potential','normalized_sources','native_boundary','native_Ainf','operator_potential','operator_normal','potential_residual','normal_residual','native_integrated_A_complex']
    for name in basis_keys:a['level1_'+name]=pd[name]
    del pd;gc.collect()
    def budget():
        memory=process_memory();r['latest_memory']=memory
        if max(memory['peak_working_set_bytes'],memory['peak_committed_bytes'])>c['maximum_peak_bytes']:raise MemoryError('native source process cap exceeded')
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('source-only total wall budget exceeded')
    def stage(name):
        budget();r['memory_records'].append(dict(stage=name,**process_memory()));write_report(args.output,r)
    stage('preflight');checkpoint_arrays(args.output,r,a)
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        if any(row['num_threads']!=1 for row in threadpool_info()):raise ValueError('one native thread required')
        r['native_thread_pools']=threadpool_info();green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        if green.exportable_settings!=pilot['green_settings']:raise ValueError('infinite Green changed')
        r['green_settings']=green.exportable_settings;scale=np.sqrt(fd['dry_mass'].diagonal()[:,None]*fd['dry_mass'].diagonal()[None,:])
        for level,n,plan in zip(c['levels'],c['counts'],plans,strict=True):
            gc.collect();memory=process_memory()
            if max(memory['private_bytes'],memory['working_set_bytes'])>c['baseline_allowance_bytes']:raise MemoryError('dense-stage baseline exceeds modeled allowance')
            if available_memory()<1.25*plan['estimated_peak_bytes']:raise MemoryError('fresh per-mesh host margin fails')
            triangles=fd[f'reference_triangles_{level}']
            if triangles.shape!=(n,3,3):raise ValueError('full source mesh missing')
            record=dict(level=level,panels=n,complete=False,solved_cases=0,modes_completed=0,surface_proof=verify_surface(fd['wetted_triangles'],triangles));r['levels'].append(record)
            before=r['completed_capytaine_calls']
            def progress(arrays,completed):
                for key,value in arrays.items():a[f'level{level}_'+key]=value
                record['modes_completed']=completed
                checkpoint_arrays(args.output,r,a)
            def on_solved(completed):
                record['solved_cases']=completed;r['completed_capytaine_calls']=before+completed;write_report(args.output,r)
            arrays,checks=export_basis(cpt,green,triangles,fd['world_cg'],mean.body.mass,rd['Ainf'][level-1],scale,budget,progress,on_solved)
            for key,value in arrays.items():a[f'level{level}_'+key]=value
            record.update(checks,complete=True);r['covered_source_levels'].append(level)
            del arrays;gc.collect();checkpoint_arrays(args.output,r,a);stage(f'level{level}_verified')
            print(f'fixed source export: level{level}, {r["completed_capytaine_calls"]}/18 new cases',flush=True)
        if r['completed_capytaine_calls']!=18 or r['covered_source_levels']!=[1,2,3,4] or not all(row['complete'] for row in r['levels']):raise ValueError('missing required source cases')
        verify_arrays(args.output.with_suffix('.npz'),a,r['arrays_sha256']);stage('array_verification')
        r['payload_package']=package_payload(args.output,archive,r,budget);stage('package_verification')
    for name,h in identities.items():
        if digest(Path(package.locate_file(name)))!=h:raise ValueError('dependency changed during export')
    disk=sum(file.stat().st_size for file in archive.rglob('*') if file.is_file())+args.output.with_suffix('.npz').stat().st_size+args.output.with_suffix('.payload.zip').stat().st_size+65536
    if disk>c['maximum_additional_disk_bytes']:raise OSError('source artifact disk cap exceeded')
    budget();r.update(complete=True,passed=False,source_export_complete=True,total_source_cases_represented=24,accounted_disk_upper_bytes=disk,elapsed_seconds=perf_counter()-STARTED)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/fixed-source-export-controls.json'))
    for name in ('pilot','pilot-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'fixed-seed-infinite-source-export');a={}
    r.update(source_export_complete=False,completed_capytaine_calls=0,common_grid_evaluations=0);checkpoint_arrays(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,source_export_complete=False,error=dict(type=type(error).__name__,message=str(error)))
    write_report(args.output,r);print(json.dumps({key:r.get(key) for key in ('complete','passed','source_export_complete','completed_capytaine_calls','total_source_cases_represented','latest_memory','elapsed_seconds','error')}))
    return 0 if r['complete'] and r.get('source_export_complete') else 1


if __name__=='__main__':raise SystemExit(main())
