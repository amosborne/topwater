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


def generate(args,c,r,a):
    fixed=dict(format_version=1,source_level=2,finer_level=3,source_panels=864,
        integration_refinements=[0,1,2,3],required_target_counts=[864,3456,13824,55296],evaluation_chunk_points=100,
        density_kg_m3=1025.,gravity_m_s2=9.81,trace_scaled_tolerance=1e-10,ownership_absolute_m=1e-10,
        reference_relative=.001,reference_absolute_normalized=1e-6,maximum_elapsed_seconds=180.,cpu_threads=1,full_reference_scope=False)
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
        scope=__doc__,trace_normalization=parent['trace_normalization'])
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
        r['green_settings']=green.exportable_settings;engine=cpt.BEMSolver(green_function=green).engine
        scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:]);matrices=[]
        for level,count in zip(c['integration_refinements'],c['required_target_counts'],strict=True):
            budget();cells,owner=integration_cells(triangles,level)
            points=body.mesh.faces_centers.copy() if level==0 else cells.mean(axis=1)
            if len(points)!=count:raise ValueError('missing required integration targets')
            proof=owned_points(points,owner,triangles,c['ownership_absolute_m'])
            weights=body.mesh.faces_areas[owner]/4**level
            if abs(weights.sum()-body.mesh.faces_areas.sum())>1e-10:raise ValueError('quadrature area conservation failed')
            phi=np.zeros((len(points),6),complex);prefix=f'level{level}_'
            a.update({prefix+'points':points,prefix+'owners':owner,prefix+'weights':weights,prefix+'potential':phi})
            save(args.output,r,a);begin=perf_counter()
            for start in range(0,len(points),c['evaluation_chunk_points']):
                budget();end=min(start+c['evaluation_chunk_points'],len(points))
                S=engine.build_S_matrix(points[start:end],body.mesh,free_surface=0.,water_depth=np.inf,wavenumber=np.inf)
                phi[start:end]=S@sigma;del S;budget()
            if not np.isfinite(phi).all():raise ValueError('nonfinite frozen-source trace')
            matrix=pressure_matrix(points,-body.mesh.faces_normals[owner],weights,phi,data['world_cg'],c['density_kg_m3'])
            imaginary=scaled_check(matrix.imag,np.zeros((6,6)),scale,c['trace_scaled_tolerance'])
            record=dict(level=level,points=len(points),seconds=perf_counter()-begin,ownership=proof,imaginary_ratio=imaginary)
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
    for name in ('study','bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'frozen-source-pressure-quadrature');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a);print(json.dumps({k:r.get(k) for k in ('complete','passed','completed_levels','quadrature_stabilization','elapsed_seconds','error')}))
    return 0 if r['complete'] and r['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
