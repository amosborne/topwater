"""Actual waterline cut rows: preserve boundary work and the original P1 field.

Float integration-reference geometry only. Strong on-panel queries own their
original source parent; all normal-flux contractions are common-plane identities.
No hidden physical mass, cut-state evolution or body-pressure power claim.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor,lu_solve
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.constants import GRAVITY,WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.capytaine_reference import body_from_triangles
from tests.benchmark_rankine import peak_rss_bytes
from tests.laplace_reference import triangle_samples
from .boat_weak_exterior_moments import owned
from .boat_local_radiation_trace import phases,amplitude,far_trace_spectrum
from .compact_p1_weak import WeakPatch
from .compact_p1_trace import CompactP1Trace
from .boat_radiation_time import harmonic_maps
from .boat_spatial_work import source_matrix
from .cut_hat_geometry import cut_hat
from .plane_p1_field import PlaneP1Field
from .panel_gradient_reference import checked_influence
from .cuda_operator import cp


def hat_value(point,centre,hx,hy):
    a,b=(point-centre)/[hx,hy]
    return max(0.,1-max(abs(a),abs(b),abs(a-b)))


def strongest_node(point,x,y,hx,hy):
    col=int(np.searchsorted(x,point[0])-1);row=int(np.searchsorted(y,point[1])-1)
    candidates=[(j,i) for j in (row,row+1) for i in (col,col+1)
                if 0<j<len(y)-1 and 0<i<len(x)-1]
    if not candidates: raise ValueError('waterline must fit inside the local patch')
    return max(candidates,key=lambda rc:hat_value(point,np.array([x[rc[1]],y[rc[0]]]),hx,hy))


def triangle_rule(parts,order):
    if not len(parts['triangles']):
        return np.empty((0,2)),np.empty(0),np.empty(0),np.empty((0,2)),np.empty(0,dtype=int)
    triangles=np.concatenate((parts['triangles'],np.zeros((*parts['triangles'].shape[:2],1))),axis=2)
    points,weights=triangle_samples(triangles,order)
    bary,_=triangle_samples(np.eye(3)[None],order)
    values=np.einsum('qj,tj->tq',bary[0],parts['values'])
    return (points[:,:,:2].reshape(-1,2),(parts['areas'][:,None]*weights).ravel(),values.ravel(),
        np.repeat(parts['gradients'],len(weights),axis=0),np.repeat(parts['parent_id'],len(weights)))


def fourier_values(spectrum,points,origin,length):
    size=spectrum.shape[0];columns=spectrum.shape[2]
    k=2*np.pi*np.fft.fftfreq(size,d=length/size)
    result=np.empty((len(points),columns))
    for first in range(0,len(points),128):
        p=points[first:first+128]
        ey=np.exp(1j*(p[:,1]-origin[1])[:,None]*k)
        ex=np.exp(1j*(p[:,0]-origin[0])[:,None]*k)
        rows=(ey@spectrum.reshape(size,-1)).reshape(len(p),size,columns)
        result[first:first+len(p)]=np.sum(rows*ex[:,:,None],axis=1).real/(size*size)
    return result


def run(args,report):
    import capytaine as cpt
    parent=owned(args.snapshot);field_meta=owned(args.field);audit=owned(args.audit)
    if (parent['parent_field_sha256']!=sha256(args.field.read_bytes()).hexdigest()
            or parent['parent_audit_sha256']!=sha256(args.audit.read_bytes()).hexdigest()):
        raise ValueError('snapshot/reference ownership differs')
    inputs=SimulationInputs.from_json((args.audit.with_suffix('.source')/'inputs.json').read_text())
    identity=sha256(inputs.to_json().encode()).hexdigest()
    if any(record['input_sha256']!=identity for record in (parent,field_meta,audit)):
        raise ValueError('physical inputs differ')
    with np.load(args.snapshot.with_suffix('.npz')) as z:
        x=z['patch_x'].copy();y=z['patch_y'].copy();delta=z['residual'].copy();density=z['enriched_density'].copy()
    with np.load(args.audit.with_suffix('.npz')) as z:
        triangles=z['triangles'].copy();external_sources=z['sources_0'].copy()
    panel_hash=sha256(triangles.tobytes()).hexdigest()
    if any(record['panel_sha256']!=panel_hash for record in (parent,field_meta,audit)):
        raise ValueError('hull panels differ')
    grid=parent['grid'];length=parent['length_m'];cg=inputs.initial_state.position
    with np.load(args.field.with_suffix('.npz')) as z: trace=z['surface_0'].reshape(grid,grid,2)
    hx,hy=float(x[1]-x[0]),float(y[1]-y[0]);shape=(len(y)-2,len(x)-2)
    if (not np.allclose(np.diff(x),hx,atol=128*np.finfo(float).eps,rtol=0.)
            or not np.allclose(np.diff(y),hy,atol=128*np.finfo(float).eps,rtol=0.)):
        raise ValueError('uniform local P1 mesh required')
    trace_field=PlaneP1Field(x,y,delta,args.maximum_pairs)
    weak=WeakPatch(shape,hx,hy)
    local_full=np.stack([cp.asnumpy(weak.apply(cp.asarray(delta[1:-1,1:-1,c]))) for c in range(4)],axis=-1)
    origin=cg[:2]-length/2;beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    inner=np.max(abs(triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam
    far=far_trace_spectrum(trace,origin,length,cg[:2],inner,beam)
    panels=TrianglePanels(triangles,-1)
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    q=np.column_stack((modes[:,[2,3]],np.zeros((len(panels),2))))
    _,body_normal,_=harmonic_maps(panels.centroids,panels.normals,grid,length/grid,origin)
    local_normal=np.column_stack([np.sum(CompactP1Trace(x,y,delta[:,:,c]).bulk(panels.centroids)[1]
        *panels.normals,axis=1) for c in range(4)])
    reconstructed=lu_solve(lu_factor(source_matrix(panels)),q-(body_normal@far.reshape(-1,4)).real-local_normal)
    np.testing.assert_allclose(reconstructed,density,atol=1e-11,rtol=1e-10)
    kx,ky=np.meshgrid(2*np.pi*np.fft.fftfreq(grid,d=length/grid),2*np.pi*np.fft.fftfreq(grid,d=length/grid))
    vertical_spectrum=far*np.hypot(kx,ky)[:,:,None]
    wet=hydrostatics(inputs.body,inputs.initial_state,density=WATER_DENSITY,gravity=GRAVITY)
    footprint=ConvexHull(wet.waterline_segments[:,:,:2].reshape(-1,2))
    polygon=footprint.points[footprint.vertices]
    edge=np.roll(polygon,-1,axis=0)-polygon
    lengths=np.linalg.norm(edge,axis=1);longest=int(np.argmax(lengths))
    midpoint=(polygon[longest]+polygon[(longest+1)%len(polygon)])/2
    normal=np.array([edge[longest,1],-edge[longest,0]])/lengths[longest]
    incoming=polygon-np.roll(polygon,1,axis=0)
    turns=np.arctan2(incoming[:,0]*edge[:,1]-incoming[:,1]*edge[:,0],np.sum(incoming*edge,axis=1))
    corner=polygon[int(np.argmax(turns))]
    straight_rc=strongest_node(midpoint,x,y,hx,hy)
    corner_rc=strongest_node(corner,x,y,hx,hy)
    exterior_rc=strongest_node(midpoint+3*max(hx,hy)*normal,x,y,hx,hy)
    xx,yy=np.meshgrid(x[1:-1],y[1:-1]);centres=np.column_stack((xx.ravel(),yy.ravel()))
    corners=np.array([[hx,0.],[hx,hy],[0.,hy],[-hx,0.],[-hx,-hy],[0.,-hy]])
    signed=(centres[:,None,:]+corners[None])@footprint.equations[:,:2].T+footprint.equations[:,2]
    hidden_candidates=np.flatnonzero(np.max(signed,axis=(1,2))<0)
    if not len(hidden_candidates): raise ValueError('no fully covered verification hat')
    hid=hidden_candidates[np.argmin(np.linalg.norm(centres[hidden_candidates]-cg[:2],axis=1))]
    hidden_rc=(int(hid//shape[1])+1,int(hid%shape[1])+1)
    cases=[('straight waterline',straight_rc),('waterline corner',corner_rc),
           ('fully exterior',exterior_rc),('fully covered',hidden_rc)]
    geometries=[]
    for name,(row,col) in cases:
        geometry=cut_hat([x[col],y[row]],hx,hy,polygon)
        parents=trace_field.source_parents(row,col,geometry)
        if name.startswith('fully'):
            expected=hx*hy if name=='fully exterior' else 0.
            np.testing.assert_allclose(geometry['physical_hat_mass'],expected,atol=1e-14,rtol=1e-12)
        elif not 0<geometry['physical_hat_mass']<hx*hy:
            raise ValueError('selected contact hat does not straddle the actual waterline')
        geometries.append((geometry,parents))
    omega=2*np.pi*field_meta['frequency_hz'];body=body_from_triangles(triangles,cg,inputs.body.mass)
    settings=dict(audit['rows'][0]['green_settings']);settings.pop('green_function')
    cache=Path('artifacts/capytaine-surface-audit-cache');cache.mkdir(exist_ok=True)
    green=cpt.Delhommeau(**settings,tabulation_cache_dir=str(cache))
    solver=cpt.BEMSolver(green_function=green)
    # Check the owned native free-space Green convention before transplanting
    # its source coefficients into our positive-Rankine halfspace map.
    sign_probe=np.array([[cg[0]+2*beam,cg[1],0.]])
    native_rankine,_=green.evaluate(sign_probe,body.mesh,free_surface=np.inf,
        early_dot_product=False,diagonal_term_in_double_layer=False)
    positive_rankine=checked_influence(sign_probe,panels).single_layer
    if not (np.all(native_rankine.real<0) and np.all(positive_rankine>0)):
        raise ValueError('external native Rankine sign does not match the source conversion')
    external_halfspace_density=phases(-external_sources/(-1j*omega))
    external=[SimpleNamespace(body=body,sources=external_sources[:,c],free_surface=0.,water_depth=np.inf,
        encounter_wavenumber=omega**2/GRAVITY,forward_speed=0.) for c in range(2)]
    report.update(input_sha256=identity,panel_sha256=panel_hash,patch_shape=shape,
        body_density_reconstruction_max_error=float(np.max(abs(reconstructed-density))),
        external_rankine_sign='native negative S; positive halfspace unit-velocity density=-sources/(-i*omega)',
        native_rankine_relative_difference=float(np.linalg.norm(native_rankine+positive_rankine)/np.linalg.norm(positive_rankine)),
        maximum_source_pairs=args.maximum_pairs,orders=args.orders,cases=[])
    saved=dict(polygon=polygon,x=x,y=y)
    for case_id,((name,(row,col)),(geometry,parents)) in enumerate(zip(cases,geometries,strict=True)):
        record=dict(name=name,node=[row,col],centre=[float(x[col]),float(y[row])],
            physical_hat_mass=geometry['physical_hat_mass'],full_hat_mass=hx*hy,
            interior_triangles=len(geometry['interior']['triangles']),exterior_triangles=len(geometry['exterior']['triangles']),
            waterline_segments=len(geometry['waterline']['segments']),zero_area_remnants=len(geometry['tangencies']),quadrature=[])
        report['cases'].append(record)
        full=local_full[row-1,col-1]
        report['active_case']=name
        for order in args.orders:
            started=perf_counter();print(json.dumps(dict(case=name,order=order,start=True)),flush=True)
            ip,iw,iv,ig,iids=triangle_rule(geometry['interior'],order)
            ep,ew,ev,_,eids=triangle_rule(geometry['exterior'],order)
            if len(ip):
                a=trace_field.vector_potential(ip,parents[iids])
                hidden_gradient=2*np.einsum('m,mi,mic->c',iw,ig,a)
                hidden_direct=(iw*iv)@trace_field.normal_flux(ip,parents[iids])
            else: hidden_gradient=np.zeros(4);hidden_direct=np.zeros(4)
            line=geometry['waterline'];waterline=np.zeros(4)
            if len(line['segments']):
                nodes,weights=leggauss(order);t=(nodes+1)/2;weights=weights/2
                line_points=(line['segments'][:,:1]+t[None,:,None]*(line['segments'][:,1:]-line['segments'][:,:1])).reshape(-1,2)
                line_values=((1-t)*line['values'][:,:1]+t*line['values'][:,1:]).ravel()
                line_length=np.linalg.norm(line['segments'][:,1]-line['segments'][:,0],axis=1)
                w=(line_length[:,None]*weights).ravel()*line_values
                a=trace_field.vector_potential(line_points,np.repeat(parents[line['parent_id']],order))
                n=np.repeat(line['hole_outward_normals'],order,axis=0)
                waterline=2*np.einsum('m,mi,mic->c',w,n,a)
            corrected=full-hidden_gradient+waterline
            source=np.zeros(4);transplanted_source=np.zeros(4)
            global_far=np.zeros(4);reference=np.zeros(2,dtype=complex)
            if len(ep):
                direct=(ew*ev)@trace_field.normal_flux(ep,parents[eids])
                global_far=(ew*ev)@fourier_values(vertical_spectrum,ep,origin,length)
                for first in range(0,len(ep),128):
                    last=min(first+128,len(ep));p=np.column_stack((ep[first:last],np.zeros(last-first)))
                    vertical=2*checked_influence(p,panels).gradient_single_layer[:,:,2]
                    source+=(ew[first:last]*ev[first:last])@(vertical@density)
                    transplanted_source+=(ew[first:last]*ev[first:last])@(vertical@external_halfspace_density)
                with threadpool_limits(1):
                    for c,result in enumerate(external):
                        velocity=np.concatenate([solver.compute_velocity(np.column_stack((ep[first:first+256],np.zeros(len(ep[first:first+256])))),result)[:,2]
                            for first in range(0,len(ep),256)])/(-1j*omega)
                        reference[c]=(ew*ev)@velocity
            else: direct=np.zeros(4)
            difference=amplitude((corrected-direct)[None])[0]
            total=amplitude((global_far+corrected+source)[None])[0]
            reference_source=amplitude(transplanted_source[None])[0]
            reference_lift=reference-reference_source
            lift_defect=amplitude((global_far+corrected)[None])[0]-reference_lift
            density_defect=amplitude((source-transplanted_source)[None])[0]
            np.testing.assert_allclose(lift_defect+density_defect,total-reference,atol=1e-14,rtol=1e-12)
            entry=dict(order=order,seconds=perf_counter()-started,
                global_far=global_far.tolist(),current_source=source.tolist(),
                transplanted_reference_source=transplanted_source.tolist(),
                inferred_reference_lift_real=reference_lift.real.tolist(),inferred_reference_lift_imag=reference_lift.imag.tolist(),
                lift_moment_defect_real=lift_defect.real.tolist(),lift_moment_defect_imag=lift_defect.imag.tolist(),
                density_response_defect_real=density_defect.real.tolist(),density_response_defect_imag=density_defect.imag.tolist(),
                lift_moment_defect_absolute=abs(lift_defect).tolist(),density_response_defect_absolute=abs(density_defect).tolist(),
                full_local=full.tolist(),hidden_gradient=hidden_gradient.tolist(),waterline=waterline.tolist(),
                corrected_local=corrected.tolist(),direct_local=direct.tolist(),hidden_direct=hidden_direct.tolist(),
                local_restriction_absolute_difference=abs(difference).tolist(),
                hidden_ibp_absolute_difference=abs(amplitude((hidden_direct-hidden_gradient+waterline)[None])[0]).tolist(),
                full_direct_closure_absolute_difference=abs(amplitude((hidden_direct+direct-full)[None])[0]).tolist(),
                omission_of_waterline_absolute_effect=abs(amplitude(waterline[None])[0]).tolist(),
                total_reference_moment_real=reference.real.tolist(),total_reference_moment_imag=reference.imag.tolist(),
                total_absolute_moment_defect=abs(total-reference).tolist(),
                area_normalized_moment_defect=(abs(total-reference)/geometry['physical_hat_mass']).tolist()
                    if geometry['physical_hat_mass']>0 else None,
                source_pairs_so_far=trace_field.pairs)
            record['quadrature'].append(entry)
            saved[f'geometry_{case_id}_parents']=geometry['parents']
            saved[f'geometry_{case_id}_interior']=geometry['interior']['triangles']
            saved[f'geometry_{case_id}_exterior']=geometry['exterior']['triangles']
            saved[f'geometry_{case_id}_waterline']=line['segments']
            saved[f'case_{case_id}_q{order}_local']=np.array([full,hidden_gradient,waterline,corrected,direct,hidden_direct])
            saved[f'case_{case_id}_q{order}_global_source']=np.array([global_far,source,transplanted_source])
            saved[f'case_{case_id}_q{order}_reference']=np.array([reference,reference_lift])
            np.savez(args.output.with_suffix('.npz'),**saved)
            report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
            print(json.dumps(dict(case=name,**entry)),flush=True)
    report.pop('active_case',None)
    report.update(complete=True,qualified=False,scope=__doc__,source_pairs=trace_field.pairs,
        reference_field_seconds=trace_field.seconds,peak_process_rss_bytes=peak_rss_bytes(),
        limits=['Float cut integration geometry, not a moving-contact topology certificate',
            'Reference lift is inferred by subtracting transplanted source work from the actual reference velocity; this decomposition is not an independent lift oracle',
            'Waterline aligned with original P1 source edges remains explicitly unsupported by strict-interior query API',
            'No physical hidden mass, no nodal/mass-inverse rate or moving-waterline state transfer',
            'Source/weak restriction checks do not establish full pressure/body power conservation',
            'Independent fine Poisson-core validation and external reference convergence remain open',
            'No nonlinear trajectory,60-second body accuracy or throughput qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['snapshot','field','audit','output']: parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--orders',type=int,nargs='+',default=[4,8])
    parser.add_argument('--maximum-pairs',type=int,default=500_000_000)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or not args.orders or min(args.orders)<2 or max(args.orders)>16 or args.maximum_pairs<1:
        parser.error('fresh output/source archive, orders2..16 and positive source-pair budget required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={})
    for name in ['snapshot','field','audit']:
        path=getattr(args,name);shutil.copyfile(path,archive/(name+'.json'))
        report['parents'][name]=dict(path=str(path),sha256=sha256(path.read_bytes()).hexdigest())
    shutil.copyfile(args.audit.with_suffix('.source')/'inputs.json',archive/'inputs.json')
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
