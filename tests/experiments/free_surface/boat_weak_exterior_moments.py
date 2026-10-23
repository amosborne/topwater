"""Physical-exterior weak moments of frozen actual-hull heave/roll snapshots.

Only complete P1 hats separated from the convex footprint are tested. No cut
contact supports, mass inversion, hidden water or time advancement are added.
External reference velocities retain their known finite-panel uncertainty.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.constants import WATER_DENSITY,GRAVITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from tests.capytaine_reference import body_from_triangles
from tests.laplace_reference import triangle_samples
from .boat_local_radiation_trace import phases,amplitude,far_trace_spectrum
from .boat_local_trace_enrichment import tensor_values
from .boat_radiation_time import harmonic_maps
from .boat_spatial_work import source_matrix
from .compact_p1_trace import CompactP1Trace
from .compact_p1_weak import WeakPatch,hat,hat_fourier
from .panel_gradient_reference import checked_influence
from .cuda_operator import cp


def owned(path):
    result=json.loads(path.read_text())
    if not result['complete'] or sha256(path.with_suffix('.npz').read_bytes()).hexdigest()!=result['arrays_sha256']:
        raise ValueError('complete owned numerical arrays required: '+str(path))
    return result


def relative_columns(value,reference):
    return (np.linalg.norm(value-reference,axis=0)/np.linalg.norm(reference,axis=0)).tolist()


def body_poisson_moments(centres,hx,hy,panels,order):
    """Integrated body-panel adjoint; no centroid replacement."""
    samples,weights=triangle_samples(panels.triangles,order)
    result=np.empty((len(centres),len(panels)))
    nodal=np.zeros((3,3));nodal[1,1]=1.
    for row,centre in enumerate(centres):
        test_hat=CompactP1Trace(centre[0]+np.array([-hx,0.,hx]),centre[1]+np.array([-hy,0.,hy]),nodal)
        image=test_hat.bulk(samples.reshape(-1,3))[0].reshape(len(panels),-1)
        result[row]=-panels.areas*(image@weights)
    return result


def run(args,report):
    import capytaine as cpt
    parent=owned(args.snapshot);field=owned(args.field);audit=owned(args.audit)
    if (parent['parent_field_sha256']!=sha256(args.field.read_bytes()).hexdigest()
            or parent['parent_audit_sha256']!=sha256(args.audit.read_bytes()).hexdigest()):
        raise ValueError('snapshot does not own these external references')
    inputs=SimulationInputs.from_json((args.audit.with_suffix('.source')/'inputs.json').read_text())
    identity=sha256(inputs.to_json().encode()).hexdigest()
    if any(r['input_sha256']!=identity for r in (parent,field,audit)):
        raise ValueError('physical inputs differ')
    with np.load(args.snapshot.with_suffix('.npz')) as z:
        x=z['patch_x'].copy();y=z['patch_y'].copy();residual=z['residual'].copy()
        density=z['enriched_density'].copy();base_density=z['base_density'].copy()
    with np.load(args.audit.with_suffix('.npz')) as z:
        triangles=z['triangles'].copy();external_sources=z['sources_0'].copy()
    panel_hash=sha256(triangles.tobytes()).hexdigest()
    if any(r['panel_sha256']!=panel_hash for r in (parent,field,audit)):
        raise ValueError('actual hull panel ownership differs')
    grid=parent['grid'];length=parent['length_m'];cg=inputs.initial_state.position
    with np.load(args.field.with_suffix('.npz')) as z:
        trace=z['surface_0'].reshape(grid,grid,2)
    panels=TrianglePanels(triangles,-1)
    hx,hy=float(x[1]-x[0]),float(y[1]-y[0]);shape=(len(y)-2,len(x)-2)
    if (not np.allclose(np.diff(x),hx,atol=128*np.finfo(float).eps,rtol=0.)
            or not np.allclose(np.diff(y),hy,atol=128*np.finfo(float).eps,rtol=0.)
            or residual.shape!=(len(y),len(x),4) or density.shape!=(len(panels),4)
            or not np.isfinite(residual).all() or not np.isfinite(density).all()):
        raise ValueError('uniform owned patch and four real phase columns required')
    origin=cg[:2]-length/2;beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    inner=np.max(abs(triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam
    spectrum=np.fft.fft2(phases(trace.reshape(-1,2)).reshape(grid,grid,4),axes=(0,1))
    far_spectrum=far_trace_spectrum(trace,origin,length,cg[:2],inner,beam)
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    q=np.column_stack((modes[:,[2,3]],np.zeros((len(panels),2))))
    _,body_normal,_=harmonic_maps(panels.centroids,panels.normals,grid,length/grid,origin)
    matrix=source_matrix(panels)
    factor=lu_factor(matrix)
    # Reconstruct the original current BVP to detect a changed trace/phase
    # policy. This test is independent of the weak-row assembly below.
    reconstructed=lu_solve(factor,q-(body_normal@spectrum.reshape(-1,4)).real)
    np.testing.assert_allclose(reconstructed,base_density,atol=1e-11,rtol=1e-10)
    local_normal=[]
    for column in range(4):
        patch=CompactP1Trace(x,y,residual[:,:,column])
        velocity=patch.bulk(panels.centroids)[1]
        local_normal.append(np.sum(velocity*panels.normals,axis=1))
    reconstructed=lu_solve(factor,q-(body_normal@far_spectrum.reshape(-1,4)).real-np.array(local_normal).T)
    np.testing.assert_allclose(reconstructed,density,atol=1e-11,rtol=1e-10)

    star,_,_,_,_=hat(hx,hy)
    xx,yy=np.meshgrid(x[1:-1],y[1:-1]);centres=np.column_stack((xx.ravel(),yy.ravel()))
    wet=hydrostatics(inputs.body,inputs.initial_state,density=WATER_DENSITY,gravity=GRAVITY)
    footprint=ConvexHull(wet.waterline_segments[:,:,:2].reshape(-1,2))
    vertices=centres[:,None,:]+star.triangles[None,:,:,:2].reshape(1,-1,2)
    # One separating halfplane must contain the ENTIRE complete star outside
    # the footprint. Merely having exterior vertices is insufficient.
    signed=vertices@footprint.equations[:,:2].T+footprint.equations[:,2]
    clearance=np.max(np.min(signed,axis=1),axis=1)
    candidates=np.flatnonzero((clearance>=.03*beam)&(clearance<=.5*beam)
        &np.all(abs(centres-cg[:2])<inner,axis=1))
    if not len(candidates): raise ValueError('no wholly physical exterior hats')
    selected=candidates[np.linspace(0,len(candidates)-1,min(48,len(candidates)),dtype=int)]
    centres=centres[selected]
    kx,ky=np.meshgrid(2*np.pi*np.fft.fftfreq(grid,d=length/grid),2*np.pi*np.fft.fftfreq(grid,d=length/grid))
    multiplier=np.hypot(kx,ky)*hat_fourier(hx,hy,kx,ky)
    coarse_all=tensor_values(spectrum*multiplier[:,:,None],x[1:-1],y[1:-1],origin,length)
    far_all=tensor_values(far_spectrum*multiplier[:,:,None],x[1:-1],y[1:-1],origin,length)
    coarse=coarse_all.reshape(-1,4)[selected];far=far_all.reshape(-1,4)[selected]
    operator=WeakPatch(shape,hx,hy)
    local=np.column_stack([cp.asnumpy(operator.apply(cp.asarray(residual[1:-1,1:-1,c]))).ravel()[selected] for c in range(4)])

    # Independent direct test of the closed-form hat Fourier transform,
    # including oblique modes that distinguish the diagonal's sign.
    samples,weights=triangle_samples(star.triangles,12)
    a,b=samples[:,:,0]/hx,samples[:,:,1]/hy
    test_values=np.maximum(0.,1-np.maximum.reduce([abs(a),abs(b),abs(a-b)]))
    test_weight=(star.areas[:,None]*weights*test_values).ravel()
    modes_test=np.array([[0.,0.],[1.3/hx,.7/hy],[-.9/hx,1.4/hy]])
    direct=np.exp(1j*samples.reshape(-1,3)[:,:2]@modes_test.T).T@test_weight
    formula=hat_fourier(hx,hy,modes_test[:,0],modes_test[:,1])
    np.testing.assert_allclose(direct,formula,atol=2e-14,rtol=2e-12)

    omega=2*np.pi*field['frequency_hz']
    external_body=body_from_triangles(triangles,cg,inputs.body.mass)
    settings=dict(audit['rows'][0]['green_settings']);settings.pop('green_function')
    cache=Path('artifacts/capytaine-surface-audit-cache');cache.mkdir(exist_ok=True)
    solver=cpt.BEMSolver(green_function=cpt.Delhommeau(**settings,tabulation_cache_dir=str(cache)))
    external=[SimpleNamespace(body=external_body,sources=external_sources[:,c],free_surface=0.,
        water_depth=np.inf,encounter_wavenumber=omega**2/GRAVITY,forward_speed=0.) for c in range(2)]
    report.update(input_sha256=identity,panel_sha256=panel_hash,patch_shape=shape,spacing_m=[hx,hy],
        tested_hats=len(selected),support_clearance_min_m=float(clearance[selected].min()),
        whole_support_certificate='One convex-footprint halfplane separates every vertex of the complete six-triangle star',
        body_density_reconstruction_max_error=float(np.max(abs(reconstructed-density))),
        hat_fourier_max_error=float(np.max(abs(direct-formula))),
        operator_diagnostics=operator.diagnostics,quadrature=[])
    saved=dict(centres=centres,selected=selected,clearance=clearance[selected],coarse=coarse,far=far,local=local)
    source_results=[];reference_results=[];body_results=[]
    for order in [4,8]:
        started=perf_counter()
        samples,weights=triangle_samples(star.triangles,order)
        a,b=samples[:,:,0]/hx,samples[:,:,1]/hy
        test_values=np.maximum(0.,1-np.maximum.reduce([abs(a),abs(b),abs(a-b)]))
        quadrature_weights=(star.areas[:,None]*weights*test_values).ravel()
        source_moment=np.empty((len(centres),len(panels)))
        reference_moment=np.empty((len(centres),2),dtype=complex)
        potential_moment=np.empty_like(reference_moment)
        body_moment=body_poisson_moments(centres,hx,hy,panels,order)
        for row,centre in enumerate(centres):
            points=samples.reshape(-1,3)+np.r_[centre,0.]
            source=np.empty((len(points),len(panels)))
            for first in range(0,len(points),128):
                last=min(first+128,len(points))
                source[first:last]=2*checked_influence(points[first:last],panels).gradient_single_layer[:,:,2]
            source_moment[row]=quadrature_weights@source
            with threadpool_limits(1):
                for column,result in enumerate(external):
                    velocity=np.concatenate([solver.compute_velocity(points[first:first+256],result)[:,2]
                        for first in range(0,len(points),256)])/(-1j*omega)
                    potential=np.concatenate([solver.compute_potential(points[first:first+256],result)
                        for first in range(0,len(points),256)])/(-1j*omega)
                    reference_moment[row,column]=quadrature_weights@velocity
                    potential_moment[row,column]=quadrature_weights@potential
            if row%12==0: print(json.dumps(dict(order=order,completed_hats=row+1,total=len(centres))),flush=True)
        source_results.append(source_moment);reference_results.append(reference_moment);body_results.append(body_moment)
        source_phase=source_moment@density
        reference=reference_moment
        row=dict(order=order,seconds=perf_counter()-started,
            source_poisson_reciprocity_max_difference=float(np.max(abs(source_moment-body_moment))),
            source_poisson_reciprocity_density_contracted_max_difference=np.max(abs(amplitude((source_moment-body_moment)@density)),axis=0).tolist(),
            reference_rms_volume_flux=np.sqrt(np.mean(abs(reference)**2,axis=0)).tolist(),cases=[])
        row['reference_robin_relation_relative_l2']=relative_columns(reference,omega**2/GRAVITY*potential_moment)
        for name,result in [('coarse global',amplitude(coarse+source_moment@base_density)),
                            ('weak local plus coarse far plus body',amplitude(far+local+source_phase))]:
            row['cases'].append(dict(method=name,relative_l2=relative_columns(result,reference),
                versus_reference_robin_relation_relative_l2=relative_columns(result,omega**2/GRAVITY*potential_moment),
                rms_absolute_volume_flux_error=np.sqrt(np.mean(abs(result-reference)**2,axis=0)).tolist(),
                rms_area_normalized_error_m_per_s=(np.sqrt(np.mean(abs(result-reference)**2,axis=0))/(hx*hy)).tolist(),
                max_absolute_volume_flux_error=np.max(abs(result-reference),axis=0).tolist()))
        report['quadrature'].append(row)
        saved.update({f'source_{order}':source_moment,f'reference_{order}':reference_moment,
            f'reference_potential_{order}':potential_moment,f'body_poisson_{order}':body_moment})
        np.savez(args.output.with_suffix('.npz'),**saved)
        report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
        print(json.dumps(row),flush=True)
    started=perf_counter()
    refined_body=body_poisson_moments(centres,hx,hy,panels,16)
    report['body_only_refinement']=dict(order=16,seconds=perf_counter()-started,
        surface_order=8,
        reciprocity_max_difference=float(np.max(abs(source_results[-1]-refined_body))),
        reciprocity_density_contracted_max_difference=np.max(abs(amplitude((source_results[-1]-refined_body)@density)),axis=0).tolist(),
        change_from_body_order8=np.max(abs(amplitude((refined_body-body_results[-1])@density)),axis=0).tolist())
    saved['body_poisson_16']=refined_body
    np.savez(args.output.with_suffix('.npz'),**saved)
    report.update(quadrature_sensitivity=dict(
        reference_relative_l2=relative_columns(reference_results[0],reference_results[1]),
        source_density_contracted_max_difference=np.max(abs(amplitude((source_results[0]-source_results[1])@density)),axis=0).tolist(),
        body_poisson_density_contracted_max_difference=np.max(abs(amplitude((body_results[0]-body_results[1])@density)),axis=0).tolist()),
        contribution_rms_volume_flux=dict(coarse=np.sqrt(np.mean(abs(amplitude(far))**2,axis=0)).tolist(),
            local=np.sqrt(np.mean(abs(amplitude(local))**2,axis=0)).tolist(),
            body=np.sqrt(np.mean(abs(amplitude(source_results[-1]@density))**2,axis=0)).tolist()))
    report.update(complete=True,qualified=False,scope=__doc__,
        arrays_sha256=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest(),peak_process_rss_bytes=peak_rss_bytes(),
        limits=['Wholly exterior test hats only; no cut waterline/contact rows or global mass closure',
            'No mass inverse: area normalization is a diagnostic, not a nodal evolution rate',
            'External evaluated field retains known finite-panel/derivative uncertainty',
            'Gauss4/8 differences are sensitivity, not certified integration-error bounds',
            'Spatial RMS of complex amplitudes is not time RMS or a body-motion error',
            'Fine independent Poisson-kernel validation remains open; requested lag tolerances are not certified cancelled-moment bounds',
            'Source reciprocity alone is not full pressure/body power conservation',
            'No nonlinear trajectory, full-minute motion accuracy or throughput qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['snapshot','field','audit','output']: parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists(): parser.error('fresh output/source archive required')
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
