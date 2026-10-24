"""Frozen Capytaine-source pressure quadrature; no BEM solve or promotion."""
from time import perf_counter
STARTED=perf_counter()
import argparse
import json
from pathlib import Path
import shutil
from zipfile import ZipFile
from importlib import import_module
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from .common import prepare,save,read_owned,digest,current_sources,midpoint_children,mesh_quality
from .paired_body import scaled_check,capytaine_radius_check


def integration_cells(triangles,level):
    cells=np.asarray(triangles).copy();owners=np.arange(len(cells))
    for _ in range(level):cells=midpoint_children(cells);owners=np.tile(owners,4)
    return cells,owners


def duffy_rule(order):
    """Positive triangle weights normalized to unit area, from [0,1]^2."""
    if not isinstance(order,(int,np.integer)) or isinstance(order,bool) or order<1:raise ValueError('positive Gauss order required')
    x,w=np.polynomial.legendre.leggauss(order);x=(x+1)/2;w=w/2
    u,v=np.meshgrid(x,x,indexing='ij');wu,wv=np.meshgrid(w,w,indexing='ij')
    bary=np.stack(((1-u)*(1-v),u,(1-u)*v),axis=-1).reshape(-1,3)
    weights=(2*(1-u)*wu*wv).ravel()
    return bary,weights


def mirrored_gauss_targets(triangles,order):
    """Mirror target geometry only; caller must evaluate all full source fields."""
    t=np.asarray(triangles);half=len(t)//2
    if t.ndim!=3 or t.shape[1:]!=(3,3) or not half or len(t)!=2*half:raise ValueError('paired reflected source triangles required')
    def canonical(values):return np.array([p[np.lexsort((p[:,2],p[:,1],p[:,0]))] for p in values])
    error=float(np.max(abs(canonical(t[:half]*np.array([1,-1,1]))-canonical(t[half:]))))
    if error>1e-12 or np.min(t[:half,:,1])<-1e-12:raise ValueError('verified y0 reflected parent layout required')
    bary,relative=duffy_rule(order)
    port=np.einsum('qk,nkj->qnj',bary,t[:half]).reshape(-1,3)
    points=np.r_[port,port*np.array([1,-1,1])]
    owner=np.tile(np.arange(half),len(bary));owners=np.r_[owner,owner+half]
    weights=np.tile(np.repeat(relative,half),2)
    return points,owners,weights,error


def owned_points(points,owners,triangles,tolerance):
    p=np.asarray(points);t=np.asarray(triangles)[owners]
    normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);twice=np.linalg.norm(normal,axis=1);normal/=twice[:,None]
    distance=abs(np.einsum('ni,ni->n',p-t[:,0],normal))
    bary=np.stack([np.einsum('ni,ni->n',np.cross(t[:,i]-p,t[:,j]-p),normal)/twice
                   for i,j in ((1,2),(2,0),(0,1))],axis=1)
    if (not np.isfinite(p).all() or np.any(p[:,2]>=0) or np.any(distance>tolerance)
            or np.any(bary<=0) or np.max(abs(bary.sum(axis=1)-1))>1e-6):
        raise ValueError('strict interior points of an owned hull triangle below z0 required')
    return dict(maximum_plane_distance_m=float(distance.max()),minimum_barycentric_weight=float(bary.min()))


def pressure_matrix(points,normal,weights,potential,center,density):
    N=np.c_[normal,np.cross(points-center,normal)]
    return density*N.T@(weights[:,None]*potential)


def telescoping(coarse,fine,quadrature,scale,tolerance):
    total=fine-coarse;output=quadrature-coarse;remainder=fine-quadrature
    scaled_check(output+remainder,total,scale,tolerance)
    return total,output,remainder


def reused_trace(points,owners,weights,arrays,prefix):
    """Reuse fields only on exactly identical physical integration targets."""
    for key,value in (('points',points),('owners',owners),('weights',weights)):
        np.testing.assert_array_equal(value,arrays[prefix+key])
    phi=arrays[prefix+'potential']
    if phi.shape!=(len(points),6) or not np.isfinite(phi).all():raise ValueError('invalid reused potential')
    return phi.copy()


def generate(args,c,r,a):
    fixed=dict(format_version=1,source_level=2,finer_level=3,source_panels=864,evaluation_chunk_points=100,
        density_kg_m3=1025.,gravity_m_s2=9.81,trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,
        reference_relative=.001,reference_absolute_normalized=1e-6,cpu_threads=1,full_reference_scope=False)
    rule=c.get('integration_rule','composite-centroid')
    if rule=='composite-centroid':
        fixed.update(integration_refinements=[0,1,2,3],required_target_counts=[864,3456,13824,55296],maximum_elapsed_seconds=180.)
        parameters=c.get('integration_refinements')
    elif rule=='duffy-gauss-mirrored':
        if c.get('gauss_orders')==[0,12,16,20]:
            fixed.update(gauss_orders=[0,12,16,20],required_target_counts=[864,124416,221184,345600],maximum_elapsed_seconds=400.,
                reuse_orders=[12,16],reuse_bundle_sha256='28a35794d8c20d35e59f955e0b4ced8fec3cdebd4c6f8e01beffd1cb205fdbb0',
                reuse_controls_sha256='bbe5c441da37a98188a213f4eef996afd99eab7be1ad7d42224ec6d95283825b')
        else:fixed.update(gauss_orders=[0,8,12,16],required_target_counts=[864,55296,124416,221184],maximum_elapsed_seconds=600.)
        parameters=c.get('gauss_orders')
    else:raise ValueError('undeclared integration rule')
    if any(c.get(k)!=v for k,v in fixed.items()):raise ValueError('unsupported frozen-source quadrature controls')
    parent,data=read_owned(args.study,digest('validation/paired-body-controls.json'),'paired-body-study')
    current_sources(parent,tuple(parent['producer_sha256']))
    if (digest(args.bundle)!=c['parent_bundle_sha256'] or parent['completed_capytaine_calls']!=24
            or parent['full_reference_qualified'] is not False or parent['dof_order']!=['Surge','Sway','Heave','Roll','Pitch','Yaw']):
        raise ValueError('parent identity/coverage/scope mismatch')
    archive=args.output.with_suffix('.source')
    with ZipFile(args.bundle) as z:
        for suffix in ('.json','.npz'):
            if z.read('validation/'+args.study.stem+suffix)!=args.study.with_suffix(suffix).read_bytes():raise ValueError('parent bytes differ from pinned bundle')
        physical=z.read('validation/'+args.study.stem+'.source/inputs.json')
    (archive/'inputs.json').write_bytes(physical)
    if digest(archive/'inputs.json')!=parent['input_sha256']:raise ValueError('actual physical input differs')
    shutil.copyfile(args.bundle,archive/'parent.zip');shutil.copyfile(c['protocol'],archive/'protocol.md')
    mean=SimulationInputs.from_json(physical.decode());triangles=data['triangles_2'];sigma=data['capytaine_sources_2']
    if triangles.shape!=(864,3,3) or sigma.shape!=(864,6):raise ValueError('required frozen source layout missing')
    np.testing.assert_array_equal(mean.initial_state.position,data['world_cg'])
    a.update(source_triangles=triangles,normalized_sources=sigma,world_cg=data['world_cg'],dry_mass=data['dry_mass'],
        native_coarse_A=data['capytaine_Ainf'][2],native_finer_A=data['capytaine_Ainf'][3],native_coarse_potential=data['capytaine_aux_2'])
    r.update(parent_manifest_sha256=digest(args.study),parent_arrays_sha256=digest(args.study.with_suffix('.npz')),
        parent_bundle_sha256=digest(args.bundle),input_sha256=digest(archive/'inputs.json'),protocol_sha256=digest(c['protocol']),
        completed_levels=0,levels=[],new_bem_solves=0,full_reference_qualified=False,source_discretization_qualified=False,
        scope=__doc__,trace_normalization=parent['trace_normalization'],integration_rule=rule)
    reuse={};reused={}
    if c.get('reuse_orders'):
        if args.reuse_study is None or args.reuse_bundle is None:raise ValueError('required reused study/bundle missing')
        prior,reused=read_owned(args.reuse_study,c['reuse_controls_sha256'],'frozen-source-pressure-quadrature')
        if digest(args.reuse_bundle)!=c['reuse_bundle_sha256']:raise ValueError('reused bundle identity differs')
        with ZipFile(args.reuse_bundle) as z:
            for suffix in ('.json','.npz'):
                if z.read('validation/'+args.reuse_study.stem+suffix)!=args.reuse_study.with_suffix(suffix).read_bytes():raise ValueError('reused study differs from pinned bundle')
        for key in ('parent_manifest_sha256','parent_arrays_sha256','parent_bundle_sha256','input_sha256','trace_normalization'):
            if prior[key]!=r[key]:raise ValueError('reused physical/source identity differs: '+key)
        for key in ('source_triangles','normalized_sources','world_cg','dry_mass','native_coarse_A','native_finer_A','native_coarse_potential'):
            np.testing.assert_array_equal(reused[key],a[key])
        if (prior['completed_levels']!=4 or prior['new_bem_solves']!=0 or prior['full_reference_qualified'] is not False
                or prior['source_discretization_qualified'] is not False or prior['integration_rule']!=rule
                or [item['rule_parameter'] for item in prior['levels']]!=[0,8,12,16]):raise ValueError('reused study scope/coverage differs')
        reuse={order:'level'+str(next(item['level'] for item in prior['levels'] if item['rule_parameter']==order))+'_' for order in c['reuse_orders']}
        shutil.copyfile(args.reuse_bundle,archive/'reuse.zip')
        r.update(reused_manifest_sha256=digest(args.reuse_study),reused_arrays_sha256=digest(args.reuse_study.with_suffix('.npz')),
            reused_bundle_sha256=digest(args.reuse_bundle),reused_orders=c['reuse_orders'])
    elif args.reuse_study is not None or args.reuse_bundle is not None:raise ValueError('undeclared reused input')
    save(args.output,r,a)
    def budget():
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('frozen-source quadrature wall budget exhausted')
    budget();import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles
        if cpt.__version__!='3.0.0' or any(p['num_threads']!=1 for p in threadpool_info()):raise ValueError('pinned dependency and one native thread required')
        r['native_thread_pools']=threadpool_info();r['capytaine_version']=cpt.__version__
        body=body_from_triangles(triangles,data['world_cg'],mean.body.mass)
        r['maximum_radius_error_m']=capytaine_radius_check(body,triangles)
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None);budget()
        if green.exportable_settings!=parent['green_settings']:raise ValueError('frozen-source Green settings differ')
        if reuse and green.exportable_settings!=prior['green_settings']:raise ValueError('reused Green settings differ')
        r['green_settings']=green.exportable_settings;engine=cpt.BEMSolver(green_function=green).engine
        scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:]);matrices=[]
        for level,(parameter,count) in enumerate(zip(parameters,c['required_target_counts'],strict=True)):
            budget();geometry_error=0.
            if level==0:
                points=body.mesh.faces_centers.copy();owner=np.arange(len(triangles));relative=np.ones(len(triangles))
            elif rule=='composite-centroid':
                cells,owner=integration_cells(triangles,parameter);points=cells.mean(axis=1);relative=np.full(len(points),1/4**parameter)
            else:points,owner,relative,geometry_error=mirrored_gauss_targets(triangles,parameter)
            if len(points)!=count:raise ValueError('missing required integration targets')
            proof=owned_points(points,owner,triangles,c['ownership_absolute_m'])
            weights=body.mesh.faces_areas[owner]*relative
            if abs(weights.sum()-body.mesh.faces_areas.sum())>1e-10:raise ValueError('quadrature area conservation failed')
            face_error=float(np.max(abs(np.bincount(owner,weights,minlength=len(triangles))-body.mesh.faces_areas)))
            if face_error>1e-10:raise ValueError('individual face quadrature weights differ')
            phi=reused_trace(points,owner,weights,reused,reuse[parameter]) if parameter in reuse else np.zeros((len(points),6),complex)
            prefix=f'level{level}_'
            a.update({prefix+'points':points,prefix+'owners':owner,prefix+'weights':weights,prefix+'potential':phi})
            save(args.output,r,a);begin=perf_counter()
            for start in range(0,0 if parameter in reuse else len(points),c['evaluation_chunk_points']):
                budget();end=min(start+c['evaluation_chunk_points'],len(points))
                S=engine.build_S_matrix(points[start:end],body.mesh,free_surface=0.,water_depth=np.inf,wavenumber=np.inf)
                phi[start:end]=S@sigma;del S;budget()
            if not np.isfinite(phi).all():raise ValueError('nonfinite frozen-source trace')
            matrix=pressure_matrix(points,-body.mesh.faces_normals[owner],weights,phi,data['world_cg'],c['density_kg_m3'])
            imaginary=scaled_check(matrix.imag,np.zeros((6,6)),scale,c['trace_scaled_tolerance'])
            record=dict(level=level,rule_parameter=parameter,points=len(points),seconds=perf_counter()-begin,
                reused=parameter in reuse,
                ownership=proof,imaginary_ratio=imaginary,maximum_face_area_error_m2=face_error,
                maximum_reflection_geometry_error_m=geometry_error)
            if parameter in reuse:record['reused_load_ratio']=scaled_check(matrix,reused[reuse[parameter]+'A_complex'],scale,c['trace_scaled_tolerance'])
            if level==0:
                record['native_potential_ratio']=scaled_check(phi,data['capytaine_aux_2'],1.,c['trace_scaled_tolerance'])
                record['native_load_ratio']=scaled_check(matrix,data['capytaine_Ainf'][2],scale,c['trace_scaled_tolerance'])
            total,output,remainder=telescoping(a['native_coarse_A'],a['native_finer_A'],matrix,scale,c['trace_scaled_tolerance'])
            a.update({prefix+'A_complex':matrix,prefix+'measured_total_increment':total,
                prefix+'output_quadrature_change':output,prefix+'unresolved_remainder':remainder})
            matrices.append(matrix.real.copy());r['levels'].append(record);r['completed_levels']+=1
            save(args.output,r,a);budget();print(f'frozen-source quadrature: {level+1}/4 levels ({count} points)',flush=True)
    if r['completed_levels']!=4:raise ValueError('missing pressure integration levels')
    A=np.stack(matrices);U,gate=mesh_quality(A/scale,dict(reference_mesh_relative=c['reference_relative'],reference_mesh_absolute_normalized=c['reference_absolute_normalized']))
    previous=abs(A[-2]-A[-3])/scale;latest=abs(A[-1]-A[-2])/scale
    allow=c['reference_absolute_normalized']+c['reference_relative']*abs(A[-1]/scale)
    a.update(quadrature_A=A,quadrature_stabilization_envelope=U*scale,
        previous_ratio=previous/allow,latest_ratio=latest/allow,trend_failed=latest>previous+c['reference_absolute_normalized'])
    r.update(complete=True,passed=False,elapsed_seconds=perf_counter()-STARTED,quadrature_stabilization=gate,
        failed_entries=dict(previous=np.argwhere(previous>allow).tolist(),latest=np.argwhere(latest>allow).tolist(),
            trend=np.argwhere(a['trend_failed']).tolist()),
        interpretation='Fixed-source output-quadrature diagnosis only; remainder also contains changed collocation solution and finer-output quadrature, not pure source error')
    budget()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,default=Path('validation/frozen-source-quadrature-controls.json'))
    p.add_argument('--reuse-study',type=Path);p.add_argument('--reuse-bundle',type=Path)
    for name in ('study','bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'frozen-source-pressure-quadrature');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a);print(json.dumps({k:r.get(k) for k in ('complete','passed','completed_levels','quadrature_stabilization','elapsed_seconds','error')}))
    return 0 if r['complete'] and r['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
