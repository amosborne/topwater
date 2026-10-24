"""Read-only study003 physical/geometry/source-asset audit; no hydrodynamic promotion."""
from time import perf_counter
STARTED=perf_counter()
import argparse,json,shutil
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile
from importlib.metadata import version
import numpy as np
from topwater import SimulationInputs,WaveComponents
from topwater.frames import rotation_matrix
from topwater.rankine import TrianglePanels
from tests import reference_boat
from .common import prepare,save,digest,study_preflight
from .finest_cost_probe import owned_zip,process_memory
from .common_output_preflight import available_memory
from .reference_mesh import verify_surface
from .reference_transfinite import shape_statistics,facet_statistics,verify_patch,validate_spec
from .frozen_source_quadrature import mirrored_gauss_targets,owned_points


def source_inventory(arrays,counts):
    records=[]
    for level,n in enumerate(counts,1):
        fields={}
        for name in ('sources','potential','pressure','normal_velocity','source_potential','source_normal_velocity'):
            key=f'level{level}_{name}'
            if key not in arrays or arrays[key].shape!=(3,n,7) or not np.isfinite(arrays[key]).all():raise ValueError('missing complete finite-frequency basis: '+key)
            fields[name]=dict(shape=list(arrays[key].shape),dtype=str(arrays[key].dtype),bytes=arrays[key].nbytes)
        records.append(dict(level=level,panels=n,finite_basis=fields))
    infinite=[key for key in arrays if 'sources' in key.lower() and ('ainf' in key.lower() or 'infinite' in key.lower())]
    return dict(levels=records,ainf_density_keys=infinite,has_ainf_densities=bool(infinite),
        radiation_frequency_cases=3*6,diffraction_frequency_heading_cases=3,original_required_diffraction_cases=9,
        missing_diffraction_frequency_heading_cases=6,
        scope='Six radiation modes at each frequency; diffraction ONLY at its paired case heading. raw_Ainf is coefficients, not densities.')


def require_calm_state(current,archived):
    for key in ('position','orientation','linear_velocity','angular_velocity'):np.testing.assert_array_equal(getattr(current,key),getattr(archived,key))
    if np.any(current.linear_velocity) or np.any(current.angular_velocity):raise ValueError('canonical calm zero-velocity state required')


def generate(args,c,r,a):
    fixed=dict(format_version=1,maximum_elapsed_seconds=60.,maximum_peak_bytes=2147483648,maximum_additional_disk_bytes=2147483648,
        counts=[1232,2772,4928,7700],subdivisions=[2,3,4,5],seed_panels=308,gauss_orders=[12,16,20],
        canonical_contract_sha256='1d0be3bc7fd20c08bffdd6b7241a58ada74fb74467ce0dd8921eb089d2621ad1',
        frozen_bundle_sha256='fcaf28963be08e4443410b8c19d680f1fb81cff3c44a291c34aa69c6e1a9847d',
        raw_bundle_sha256='e18bf7366eb651766edc5487addd84c798c1ae53183f58a0e88b3b623ed6729d',
        sources_bundle_sha256='4a317cbb91cc51348b0059b40b4a6437315e3656624ad2903743dec4721214e5',full_reference_scope=False)
    if any(c.get(key)!=value for key,value in fixed.items()):raise ValueError('undeclared fixed-asset preflight controls')
    if digest('validation/canonical-boat.json')!=c['canonical_contract_sha256']:raise ValueError('canonical contract changed')
    archive=args.output.with_suffix('.source');identities={}
    for name in ('frozen','raw','sources'):
        file=getattr(args,name+'_bundle');h=digest(file)
        if h!=c[name+'_bundle_sha256']:raise ValueError('fixed asset bundle differs')
        shutil.copyfile(file,archive/(name+'.zip'));identities[name]=h
    frozen,fd=owned_zip(args.frozen_bundle,'frozen-fixed-001','freeze',c['canonical_contract_sha256'])
    raw,rd=owned_zip(args.raw_bundle,'reference-study-003','reference',c['canonical_contract_sha256'])
    sources,sd=owned_zip(args.sources_bundle,'periodic-sources-003','periodic-source-recovery',digest('validation/periodic-source-controls.json'))
    with ZipFile(args.frozen_bundle) as z:
        frozen_bytes=z.read('validation/frozen-fixed-001.json');physical=z.read('validation/frozen-fixed-001.source/inputs.json')
        contract=json.loads(z.read('validation/frozen-fixed-001.source/contract.json'));proof=json.loads(z.read('validation/frozen-fixed-001-proof.json'))
    with ZipFile(args.raw_bundle) as z:
        raw_bytes=z.read('validation/reference-study-003.json')
        if sha256(z.read('bound-freeze.zip')).hexdigest()!=identities['frozen']:raise ValueError('raw bound freeze differs')
    with ZipFile(args.sources_bundle) as z:
        if sha256(z.read('validation/periodic-sources-003.source/raw-reference.zip')).hexdigest()!=identities['raw']:raise ValueError('source bound reference differs')
        case=z.read('validation/periodic-sources-003.source/mild-case.json')
    if (not frozen['passed'] or raw['passed'] or sources['passed'] is not True or sources['full_reference_qualified'] is not False
            or raw['input_sha256']!=frozen['input_sha256'] or sha256(physical).hexdigest()!=frozen['input_sha256']
            or sources['raw_manifest_sha256']!=sha256(raw_bytes).hexdigest()
            or sources['frozen_manifest_sha256']!=sha256(frozen_bytes).hexdigest()
            or raw['frozen_physical_identity_sha256']!=frozen['physical_identity_sha256']
            or sources['frozen_physical_identity_sha256']!=frozen['physical_identity_sha256']
            or sha256(case).hexdigest()!=sources['case_sha256'] or sources['case_sha256']!=digest('validation/mild-periodic-case.json')
            or raw['completed_physics_calls']!=132 or sources['completed_physics_calls']!=84):raise ValueError('asset identity/scope/coverage differs')
    (archive/'inputs.json').write_bytes(physical);(archive/'mild-case.json').write_bytes(case);shutil.copyfile(c['protocol'],archive/'protocol.md')
    if contract['physical_problem']!=frozen['physical_problem']:raise ValueError('physical conventions differ')
    if (sources['time_convention']!=contract['physical_problem']['time_convention'] or sources['density_kg_m3']!=1025.
            or sources['gravity_m_s2']!=9.81 or sources['water_depth']!='infinite' or sources['forward_speed_m_s']!=0.
            or sources['basis_columns']!=contract['physical_problem']['mode_order']+['Diffraction']):raise ValueError('finite source conventions differ')
    if any(sources['green_settings'].get(key)!=value for key,value in contract['reference_numerics']['finite_green_settings'].items()):raise ValueError('finite Green controls differ')
    if version('capytaine')!='3.0.0':raise ValueError('archived radius convention requires pinned dependency')
    r.update(scope=__doc__,bundle_sha256=identities,input_sha256=frozen['input_sha256'],
        physical_identity_sha256=frozen['physical_identity_sha256'],protocol_sha256=digest(c['protocol']),
        historical_commits={name:report['current_commit'] for name,report in [('frozen',frozen),('raw',raw),('sources',sources)]},
        physical_problem=contract['physical_problem'],original_gates=raw['gates'],original_raw_passed=False,
        source_reproduction_passed=sources['passed'],full_reference_qualified=False,source_discretization_qualified=False,
        native_evaluations=0,new_bem_solves=0,levels=[],memory_records=[])
    def budget(stage):
        memory=process_memory();r['memory_records'].append(dict(stage=stage,**memory))
        if max(memory['peak_working_set_bytes'],memory['peak_committed_bytes'])>c['maximum_peak_bytes']:raise MemoryError('read-only preflight memory cap exceeded')
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('fixed asset preflight budget exhausted')
    budget('inputs')
    body=reference_boat.reference_body();pose,hydro=reference_boat.equilibrium();R=rotation_matrix(pose.orientation)
    old_inputs=SimulationInputs.from_json(physical.decode());old_body=old_inputs.body
    require_calm_state(pose,old_inputs.initial_state)
    for key in ('vertices','faces','opening_rim','artificial_faces'):np.testing.assert_array_equal(getattr(body.mesh,key),getattr(old_body.mesh,key))
    for key in ('mass','cg','inertia_cg'):np.testing.assert_array_equal(getattr(body,key),getattr(old_body,key))
    mass=np.zeros((6,6));mass[:3,:3]=body.mass*np.eye(3);mass[3:,3:]=R@body.inertia_cg@R.T
    for key,value in dict(world_cg=pose.position,orientation=pose.orientation,dry_mass=mass,
            hydrostatic_stiffness=hydro.stiffness,wetted_triangles=hydro.wetted_triangles,body_mass=np.array(body.mass)).items():np.testing.assert_array_equal(value,fd[key])
    for key in ('world_cg','dry_mass','hydrostatic_stiffness'):np.testing.assert_array_equal(sd[key],fd[key])
    for key in ('Ainf','A','B'):np.testing.assert_array_equal(sd['raw_'+key],rd[key])
    np.testing.assert_array_equal(sd['angular_frequencies_rad_s'],2*np.pi*np.asarray(contract['cases']['frequencies_hz']))
    case_data=json.loads(case);np.testing.assert_array_equal(sd['headings_rad'],np.deg2rad(case_data['sea']['propagation_headings_degrees']))
    r['current_physical_arrays_bitwise_equal']=True;r['current_input_json_bytes_equal']=SimulationInputs(body,WaveComponents([],[],[],[]),pose).to_json().encode()==physical
    r['finite_source_convention']={key:sources[key] for key in ('time_convention','basis_columns','basis_normalization','source_convention','green_settings')}
    r['source_inventory']=source_inventory(sd,c['counts'])
    if r['source_inventory']['has_ainf_densities']:raise ValueError('unexpected unscoped Ainf densities in finite archive')
    spec=contract['reference_numerics'];validate_spec(spec);seed=fd['reference_seed_triangles']
    if seed.shape!=(308,3,3):raise ValueError('seed incomplete')
    verify_surface(fd['wetted_triangles'],seed);a.update({key:fd[key] for key in ('wetted_triangles','world_cg','orientation','dry_mass','hydrostatic_stiffness','body_mass','hull_length','hull_beam')});a['reference_seed_triangles']=seed
    r['seed_quality']=shape_statistics(seed)
    for level,(count,n) in enumerate(zip(c['counts'],c['subdivisions'],strict=True),1):
        mesh=fd[f'reference_triangles_{level}']
        if mesh.shape!=(count,3,3) or count!=len(seed)*n*n:raise ValueError('complete fixed family missing')
        np.testing.assert_array_equal(mesh,rd[f'triangles_{level}']);np.testing.assert_array_equal(mesh,sd[f'level{level}_triangles'])
        geometry=verify_surface(fd['wetted_triangles'],mesh);half=len(seed)//2
        errors=[verify_patch(seed[i],mesh[i*n*n:(i+1)*n*n],n,spec) for i in range(half)]
        _,_,_,reflection=mirrored_gauss_targets(mesh,1)
        distances=np.linalg.norm(mesh-mesh.mean(axis=1)[:,None],axis=2);radius_error=float(np.max(abs(distances[:,0]-distances.max(axis=1))))
        if radius_error>1e-12 or proof['levels'][level-1]['radius_error_m']!=0:raise ValueError('first-vertex radius convention differs')
        r['levels'].append(dict(level=level,subdivisions=n,panels=count,geometry=geometry,quality=shape_statistics(mesh),facets=facet_statistics(fd['wetted_triangles'],mesh),
            maximum_patch_lattice_error=max(errors),maximum_reflection_error_m=reflection,maximum_first_vertex_radius_error_m=radius_error,archived_actual_radius_error_m=0.))
        a[f'reference_triangles_{level}']=mesh;budget(f'mesh{level}')
    panels=TrianglePanels(seed,-1);r['common_output_grids']=[]
    for order in c['gauss_orders']:
        points,owners,relative,error=mirrored_gauss_targets(seed,order);weights=panels.areas[owners]*relative
        proof_points=owned_points(points,owners,seed,1e-10)
        if np.any(weights<=0) or np.max(abs(np.bincount(owners,weights,minlength=308)-panels.areas))>1e-10:raise ValueError('common-grid area weights differ')
        key=f'target_{order}_';a.update({key+'points':points,key+'owners':owners,key+'weights':weights,key+'normal':panels.normals[owners]})
        r['common_output_grids'].append(dict(order=order,targets=len(points),ownership=proof_points,maximum_reflection_error_m=error));budget(f'grid{order}')
    r['source_export_workspace_plan']=study_preflight(308,[1,2,3,4],spec['maximum_estimated_workspace_bytes'],panel_counts=c['counts'])
    r['future_host_observation']=dict(available_physical_bytes=available_memory(),available_disk_bytes=shutil.disk_usage(Path.cwd()).free)
    targets=[row['targets'] for row in r['common_output_grids']]
    r['work_inventory']=dict(target_counts=targets,all_four_source_field_pairs=sum(targets)*sum(c['counts']),
        finest_field_pairs=sum(targets)*c['counts'][-1],fresh_ainf_export_calls=24,missing_finite_diffraction_export_calls=24)
    r.update(complete=True,passed=False,assets_usable=True,fresh_ainf_export_required=True,full_run_authorized=False,
        decision='Declare a fresh Ainf source-export/reproduction pilot before any controlled-output hydrodynamic claim',elapsed_seconds=perf_counter()-STARTED)
    save(args.output,r,a);budget('saved')
    disk=sum(file.stat().st_size for file in archive.rglob('*') if file.is_file())+args.output.with_suffix('.npz').stat().st_size+65536
    if disk>c['maximum_additional_disk_bytes']:raise OSError('preflight disk cap exceeded')
    r['accounted_disk_upper_bytes']=disk;r['elapsed_seconds']=perf_counter()-STARTED


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--controls',type=Path,default=Path('validation/fixed-asset-preflight-controls.json'))
    for name in ('frozen-bundle','raw-bundle','sources-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'fixed-seed-asset-preflight');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,assets_usable=False,error=dict(type=type(error).__name__,message=str(error)))
    args.output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps({key:r.get(key) for key in ('complete','passed','assets_usable','fresh_ainf_export_required','work_inventory','elapsed_seconds','error')}))
    return 0 if r['complete'] and r.get('assets_usable') else 1


if __name__=='__main__':raise SystemExit(main())
