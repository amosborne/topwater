"""Offline exterior periodic fields and checked surface-map candidate; no TD run."""
import argparse
from importlib import import_module
import json
from pathlib import Path
import shutil
from time import perf_counter
from zipfile import ZipFile
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from topwater.frames import rotation_matrix
from .common import prepare,save,read_owned,digest,current_sources,REFERENCE
from .mild import load_case,raw_inputs,initial_body
from .periodic_grid import hull_planes,radius,mean_angles,annulus,exterior_points,periodic_interpolate,mapped_surface


def incident(points,omega,heading,gravity,zeta):
    k=omega**2/gravity;d=np.array([np.cos(heading),np.sin(heading)])
    return -1j*gravity/omega*np.exp(k*points[:,2]+1j*k*(points[:,:2]@d))*zeta


def boundary_limits(values):
    """Adjacent exterior samples at descending offsets with a fixed ratio of two."""
    return 2*values[...,1:,:]-values[...,:-1,:]


def envelope(values,floor):
    previous=abs(values[-2]-values[-3]);latest=abs(values[-1]-values[-2])
    return np.maximum(previous,latest).sum(axis=0),bool(np.all(latest<=previous+floor))


def world_vertices(body,state):
    return (body.mesh.vertices-body.cg)@rotation_matrix(state.orientation).T+state.position


def generate(args,c,case,report,out):
    if (c['format_version']!=1 or c['angular_divisions']!=[32,64,128]
            or c['outer_radius_m']!=24. or c['candidate_depth_m']!=24.
            or c['radial_offsets_m']!=[0.,.15,.4,1.,2.5,6.,12.]
            or c['boundary_offsets_m']!=[.02,.01,.005,.0025]
            or c['evaluation_chunk_points']!=100 or c['cpu_threads']!=1
            or c['maximum_elapsed_seconds']!=1800):
        raise ValueError('unsupported predeclared field/grid controls')
    contract,raw,r,frozen,data=raw_inputs(args,report)
    sr,s=read_owned(args.sources,digest('validation/periodic-source-controls.json'),'periodic-source-recovery',passed=True)
    current_sources(sr,REFERENCE+('tests/validation/mild.py','tests/validation/periodic_sources.py'))
    if (sr['completed_physics_calls']!=84 or sr['raw_manifest_sha256']!=digest(args.raw)
            or sr['frozen_manifest_sha256']!=digest(args.frozen) or sr['case_sha256']!=digest(c['case'])
            or digest(args.sources_bundle)!=c['source_bundle_sha256']):
        raise ValueError('source/case/raw identities differ from archived successful recovery')
    archive=args.output.with_suffix('.source')
    with ZipFile(args.sources_bundle) as bundle:
        for suffix in ('.json','.npz'):
            if bundle.read('validation/'+args.sources.stem+suffix)!=args.sources.with_suffix(suffix).read_bytes():
                raise ValueError('source bytes differ from pinned bundle')
    shutil.copyfile(args.sources_bundle,archive/'periodic-sources.zip')
    for origin,target in ((c['case'],'mild-case.json'),(c['protocol'],'protocol.md'),
            (args.frozen.with_suffix('.source')/'inputs.json','frozen-inputs.json')):
        shutil.copyfile(origin,archive/target)
    mean=SimulationInputs.from_json((archive/'frozen-inputs.json').read_text())
    origin=data['world_cg'][:2];planes=hull_planes(world_vertices(mean.body,mean.initial_state),mean.body.mesh.faces)
    angles=mean_angles(data['wetted_triangles'],origin,128);n=len(angles)
    inner=radius(planes,angles,0.,origin)
    vertices,faces,radii=annulus(angles,inner,c['outer_radius_m'],c['radial_offsets_m'],origin)
    offsets=np.asarray(c['boundary_offsets_m'])
    boundary_xy=origin+(inner[None,:]+offsets[:,None])[:,:,None]*np.c_[np.cos(angles),np.sin(angles)][None,:,:]
    boundary=np.c_[boundary_xy.reshape(-1,2),np.zeros(4*n)]
    points=np.r_[boundary,vertices[n:]];exterior_points(points,planes)
    out.update(angles_rad=angles,mean_vertices=vertices,mean_faces=faces,reference_points=points,
        mean_hull_planes=planes,world_cg=data['world_cg'],radial_coordinates_m=radii,
        full_hull_faces=mean.body.mesh.faces,mean_full_hull_vertices=world_vertices(mean.body,mean.initial_state))
    nr=len(c['radial_offsets_m'])+1
    omega=s['angular_frequencies_rad_s'];heading=s['headings_rad'];zeta=s['incident_phasors_m'];g=sr['gravity_m_s2']
    out.update(angular_frequencies_rad_s=omega,incident_phasors_m=zeta)
    scattered=np.zeros((4,3,len(points)),complex);out['scattered_potential_samples']=scattered
    inc=np.stack([incident(points,w,h,g,z) for w,h,z in zip(omega,heading,zeta,strict=True)])
    out['incident_potential_samples']=inc
    report.update(source_manifest_sha256=digest(args.sources),source_bundle_sha256=digest(args.sources_bundle),
        case_sha256=digest(c['case']),protocol_sha256=digest(c['protocol']),
        point_count=len(points),angular_points=n,completed_component_evaluations=0,
        full_reference_qualified=False,nonlinear_restart_ready=False,original_raw_passed=raw['passed'],
        scope=__doc__,density_kg_m3=sr['density_kg_m3'],gravity_m_s2=g,
        maximum_green_block_shape=[c['evaluation_chunk_points'],max(len(data[f'reference_triangles_{j}']) for j in range(1,5))],
        field_uncertainty_scope='empirical field diagnostics only; no field-to-motion uncertainty propagation')
    # Complete actual physical graph and invalid-point checks precede native work.
    save(args.output,report,out);started=perf_counter()
    def budget():
        if perf_counter()-started>c['maximum_elapsed_seconds']:raise TimeoutError('field evaluation budget exhausted')
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        report['native_thread_pools']=threadpool_info()
        if cpt.__version__!='3.0.0' or any(p['num_threads']!=1 for p in threadpool_info()):
            raise ValueError('pinned Capytaine and single native thread required')
        save(args.output,report,out)
        green=cpt.Delhommeau(tabulation_cache_dir=None);budget()
        if green.exportable_settings!=sr['green_settings']:raise ValueError('reference Green settings differ')
        report['green_settings']=green.exportable_settings
        engine=cpt.BEMSolver(green_function=green).engine
        for level in range(4):
            triangles=data[f'reference_triangles_{level+1}']
            np.testing.assert_array_equal(triangles,s[f'level{level+1}_triangles'])
            body=body_from_triangles(triangles,data['world_cg'],float(data['body_mass']))
            for f,w in enumerate(omega):
                sigma=s[f'level{level+1}_scattered_sources'][f]
                for start in range(0,len(points),c['evaluation_chunk_points']):
                    budget();end=min(start+c['evaluation_chunk_points'],len(points))
                    matrix=engine.build_S_matrix(points[start:end],body.mesh,free_surface=0.,water_depth=np.inf,wavenumber=w*w/g)
                    scattered[level,f,start:end]=matrix@sigma;del matrix;budget()
                if not np.isfinite(scattered[level,f]).all():raise ValueError('nonfinite exterior potential')
                report['completed_component_evaluations']+=1;save(args.output,report,out)
                print(f'periodic fields: {report["completed_component_evaluations"]}/12',flush=True)
    total=scattered+inc[None,:,:];out['total_potential_samples']=total
    limits=boundary_limits(total[:,:,:4*n].reshape(4,3,4,n));out['boundary_limit_potential']=limits
    phi=np.concatenate((limits[:,:,-1:, :],total[:,:,4*n:].reshape(4,3,nr-1,n)),axis=2)
    eta=1j*omega[None,:,None,None]*phi/g
    out.update(potential_phasors=phi,elevation_phasors=eta,
        initial_potential=phi[-1].sum(axis=0).real,initial_elevation=eta[-1].sum(axis=0).real)
    floor=c['field_trend_floor_SI'];diagnostics={}
    for name,value,unitlimit in (('potential',phi,c['field_phi_diagnostic_limit_m2_s']),('elevation',eta,c['field_eta_diagnostic_limit_m'])):
        mesh,trend=envelope(value,floor)
        boundary_value=limits if name=='potential' else limits*1j*omega[None,:,None,None]/g
        prev=abs(boundary_value[-1,:,1]-boundary_value[-1,:,0]);last=abs(boundary_value[-1,:,2]-boundary_value[-1,:,1])
        bnd=np.maximum(prev,last).sum(axis=0)
        out[name+'_mesh_envelope']=mesh;out[name+'_boundary_limit_envelope']=bnd
        combined=mesh.copy();combined[0]+=bnd
        diagnostics[name]=dict(mesh_trend=trend,boundary_trend=bool(np.all(last<=prev+floor)),
            maximum_mesh_envelope=float(mesh.max()),maximum_boundary_limit_envelope=float(bnd.max()),
            maximum_combined_field_ratio=float(combined.max()/unitlimit))
    report['field_diagnostics']=diagnostics
    # Fixed nested angular controls measure representation error; no field-to-
    # motion implication follows from this interpolation comparison.
    for name,value in (('potential',phi[-1]),('elevation',eta[-1])):
        coarse=[]
        for divisions in c['angular_divisions'][:-1]:
            a=mean_angles(data['wetted_triangles'],origin,divisions)
            indices=np.array([np.argmin(abs(angles-x)) for x in a])
            if np.max(abs(angles[indices]-a))>1e-12:raise ValueError('missing nested angular samples')
            coarse.append(periodic_interpolate(a,value[:,:,indices],angles))
        previous=abs(coarse[1]-coarse[0]);last=abs(value-coarse[1])
        out[name+'_angular_envelope']=np.maximum(previous,last).sum(axis=0)
        diagnostics[name]['maximum_angular_envelope']=float(out[name+'_angular_envelope'].max())
    motion=np.array([r['RAO'][-1,f,f]*zeta[f] for f in range(3)])
    report['surface_maps']=[]
    for scale in case['amplitude_controls']['scales']:
        displacement=(motion*scale).sum(axis=0).real
        velocity=(-1j*omega[:,None]*motion*scale).sum(axis=0).real
        state=initial_body(mean.initial_state,displacement,velocity)
        actual=hull_planes(world_vertices(mean.body,state),mean.body.mesh.faces)
        entry=dict(amplitude_scale=scale,geometry_passed=False)
        try:
            if len(mean.body.mesh.opening_rim):
                rim=world_vertices(mean.body,state)[mean.body.mesh.opening_rim]
                clearance=float(rim[:,2].min()-np.max(out['initial_elevation']*scale))
                entry['conservative_opening_clearance_m']=clearance
                if clearance<=0:raise ValueError('opening rim may intersect represented surface; no artificial-cap contact allowed')
            a,v,f,checks=mapped_surface(angles,out['initial_elevation']*scale,planes,actual,origin,c['outer_radius_m'],c['radial_offsets_m'])
            prefix='scale_'+str(scale)+'_';out.update({prefix+'angles_rad':a,prefix+'vertices':v,prefix+'faces':f,
                prefix+'hull_planes':actual,prefix+'full_hull_vertices':world_vertices(mean.body,state),
                prefix+'initial_potential':periodic_interpolate(angles,out['initial_potential']*scale,a)})
            entry.update(geometry_passed=True,**checks)
        except ValueError as error:entry['error']=str(error)
        report['surface_maps'].append(entry)
    budget()
    report.update(complete=True,passed=False,elapsed_seconds=perf_counter()-started,
        required_readiness=dict(exterior_samples=True,surface_map_geometry=all(e['geometry_passed'] for e in report['surface_maps']),
            field_to_motion_uncertainty=False,body_boundary_compatibility=False,complete_fluid_domain=False,linear_limit=False),
        qualification_status='unqualified field/map candidate: field-to-motion, body compatibility, full domain and linear-limit evidence missing')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,default=Path('validation/periodic-grid-controls.json'))
    for name in ('raw','frozen','sources','sources-bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());case=load_case(c['case'])
    report=prepare(args.output,args.controls,'periodic-field-map-candidate');out={};save(args.output,report,out)
    try:generate(args,c,case,report,out)
    except (Exception,KeyboardInterrupt) as error:report.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,report,out)
    print(json.dumps({k:report.get(k) for k in ('complete','passed','completed_component_evaluations','surface_maps','error')}))
    return 0 if report['complete'] and report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
