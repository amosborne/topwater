"""Local trace enrichment of frozen external heave/roll FIELD snapshots.

All body solves use real time-phase columns. The external field's actual
evaluated velocity is the spatial reference; its known finite-panel Robin
defect is reported separately and is not treated as exact open-water physics.
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
from scipy.spatial import ConvexHull,cKDTree
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import GRAVITY,WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.capytaine_reference import body_from_triangles
from tests.benchmark_rankine import peak_rss_bytes
from .boat_hos_continuation import window
from .boat_local_trace_enrichment import tensor_values,norms
from .boat_radiation_time import harmonic_maps
from .boat_spatial_work import source_matrix
from .compact_p1_trace import CompactP1Trace
from .panel_gradient_reference import checked_influence


def phases(values):
    return np.column_stack((values.real,values.imag))


def amplitude(values):
    return values[:,:2]+1j*values[:,2:]


def run(args,report):
    import capytaine as cpt
    field=json.loads(args.field.read_text());audit=json.loads(args.audit.read_text())
    for path,metadata in [(args.field,field),(args.audit,audit)]:
        if not metadata['complete'] or sha256(path.with_suffix('.npz').read_bytes()).hexdigest()!=metadata['arrays_sha256']:
            raise ValueError('complete owned external arrays required')
    inputs=SimulationInputs.from_json((args.audit.with_suffix('.source')/'inputs.json').read_text())
    identity=sha256(inputs.to_json().encode()).hexdigest()
    if field['input_sha256']!=identity or audit['input_sha256']!=identity or field['panel_sha256']!=audit['panel_sha256']:
        raise ValueError('external snapshots differ in physical input or panels')
    grid_rows=[r for r in field['rows'] if r['trace_control']=='full reference auxiliary trace']
    grid=grid_rows[0]['grid_size'];length=grid_rows[0]['length_m']
    with np.load(args.field.with_suffix('.npz')) as z:
        trace=z['surface_0'].reshape(grid,grid,2)
        reference_force=z['reference_force'].copy()
        reference_body=z['reference_body'].copy()
    with np.load(args.audit.with_suffix('.npz')) as z:
        triangles=z['triangles'].copy();sources=z['sources_0'].copy()
        np.testing.assert_allclose(z['body_0'],reference_body,atol=1e-12,rtol=1e-12)
    if sha256(triangles.tobytes()).hexdigest()!=audit['panel_sha256']:
        raise ValueError('actual external triangles hash mismatch')
    omega=2*np.pi*.1
    body=body_from_triangles(triangles,inputs.initial_state.position,inputs.body.mass)
    settings=dict(audit['rows'][0]['green_settings']);settings.pop('green_function')
    cache=Path('artifacts/capytaine-surface-audit-cache');cache.mkdir(exist_ok=True)
    green=cpt.Delhommeau(**settings,tabulation_cache_dir=str(cache))
    solver=cpt.BEMSolver(green_function=green)
    results=[SimpleNamespace(body=body,sources=sources[:,c],free_surface=0.,water_depth=np.inf,
        encounter_wavenumber=omega**2/GRAVITY,forward_speed=0.) for c in range(2)]

    def field_at(points,velocity=False):
        values=[]
        with threadpool_limits(1):
            for result in results:
                blocks=[]
                for first in range(0,len(points),256):
                    function=solver.compute_velocity if velocity else solver.compute_potential
                    blocks.append(function(points[first:first+256],result)/(-1j*omega))
                values.append(np.concatenate(blocks))
        return np.stack(values,axis=1)

    started=perf_counter();cg=inputs.initial_state.position
    panels=TrianglePanels(triangles,-1)
    beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    weighted=panels.areas[:,None]*modes
    q=np.column_stack((modes[:,[2,3]],np.zeros((len(panels),2))))
    origin=cg[:2]-length/2
    spectrum=np.fft.fft2(phases(trace.reshape(-1,2)).reshape(grid,grid,4),axes=(0,1))
    pbody,nbody,_=harmonic_maps(panels.centroids,panels.normals,grid,length/grid,origin)
    lift=(pbody@spectrum.reshape(-1,4)).real
    normal=(nbody@spectrum.reshape(-1,4)).real
    factor=lu_factor(source_matrix(panels));single=assemble_dirichlet_halfspace(panels).single_layer
    base_density=lu_solve(factor,q-normal)
    base_body=amplitude(lift+single@base_density)
    wet=hydrostatics(inputs.body,inputs.initial_state,density=WATER_DENSITY,gravity=GRAVITY)
    inner=np.max(abs(triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam
    outer=inner+beam
    # Remove unresolved hull-scale trace detail BEFORE sampling the global
    # Fourier part. The local Poisson lift restores it. Both windows act on
    # surface DATA only; no three-dimensional harmonic field is windowed.
    coarse_axis=np.arange(grid)*length/grid
    gx,gy=np.meshgrid(origin[0]+coarse_axis,origin[1]+coarse_axis)
    excluded=window(gx-cg[0],inner[0]-.5*beam,inner[0])*window(gy-cg[1],inner[1]-.5*beam,inner[1])
    far_trace=phases(trace.reshape(-1,2)).reshape(grid,grid,4)*(1-excluded[:,:,None])
    far_spectrum=np.fft.fft2(far_trace,axes=(0,1))
    far_lift=(pbody@far_spectrum.reshape(-1,4)).real
    far_normal=(nbody@far_spectrum.reshape(-1,4)).real
    probe_cells=np.ceil(2*outer/.25).astype(int)
    cells=probe_cells*2**args.local_refinement
    x=np.linspace(cg[0]-outer[0],cg[0]+outer[0],cells[0]+1)
    y=np.linspace(cg[1]-outer[1],cg[1]+outer[1],cells[1]+1)
    xx,yy=np.meshgrid(x,y)
    nodes=np.column_stack((xx.ravel(),yy.ravel(),np.zeros(xx.size)))
    exact_trace=phases(field_at(nodes)).reshape(len(y),len(x),4)
    coarse_trace=tensor_values(far_spectrum,x,y,origin,length)
    taper=window(xx-cg[0],inner[0],outer[0])*window(yy-cg[1],inner[1],outer[1])
    residual=(exact_trace-coarse_trace)*taper[:,:,None]
    residual[[0,-1],:]=0.;residual[:,[0,-1]]=0.
    patches=[CompactP1Trace(x,y,residual[:,:,c]) for c in range(4)]
    probe_x=np.linspace(cg[0]-outer[0],cg[0]+outer[0],probe_cells[0]+1)
    probe_y=np.linspace(cg[1]-outer[1],cg[1]+outer[1],probe_cells[1]+1)
    probe_owner=CompactP1Trace(probe_x,probe_y,np.zeros((len(probe_y),len(probe_x))))
    centroids=probe_owner.panels.centroids
    footprint=ConvexHull(wet.waterline_segments[:,:,:2].reshape(-1,2))
    distance=np.max(centroids[:,:2]@footprint.equations[:,:2].T+footprint.equations[:,2],axis=1)
    candidates=np.flatnonzero((distance>=.03*beam)&(distance<=.5*beam)
        &np.all(abs(centroids[:,:2]-cg[:2])<inner,axis=1))
    if not len(candidates):
        raise ValueError('finite-clearance exterior centroid probes required')
    parent_indices=candidates[np.linspace(0,len(candidates)-1,min(96,len(candidates)),dtype=int)]
    probes=centroids[parent_indices]
    match_distance,selected=cKDTree(patches[0].panels.centroids).query(probes)
    match_allowance=256*np.finfo(float).eps*max(1.,float(abs(probes).max()))
    if np.max(match_distance)>match_allowance:
        raise ValueError('dyadic patch must retain the fixed parent centroid probes')
    angle=2*np.pi*(np.arange(24)+.37)/24
    far_xy=cg[:2]+2*max(outer)*np.column_stack((np.cos(angle),np.sin(angle)))
    targets=np.r_[probes,np.column_stack((far_xy,np.zeros(len(far_xy))))]
    reference_trace=field_at(targets);reference_vertical=field_at(targets,velocity=True)[:,:,2]
    pprobe,nprobe,_=harmonic_maps(targets,np.tile([0.,0.,1.],(len(targets),1)),grid,length/grid,origin)
    trace_base=(pprobe@spectrum.reshape(-1,4)).real
    vertical_base=(nprobe@spectrum.reshape(-1,4)).real
    trace_far=(pprobe@far_spectrum.reshape(-1,4)).real
    vertical_far=(nprobe@far_spectrum.reshape(-1,4)).real
    # Verification-only sampling audit: project a stored finer SAME-domain
    # trace after the SAME exclusion. This does not change the candidate solve.
    fine_index=next(i for i,r in enumerate(grid_rows) if r['length_m']==length and r['grid_size']>grid)
    fine_grid=grid_rows[fine_index]['grid_size']
    with np.load(args.field.with_suffix('.npz')) as z:
        fine_trace=z[f'surface_{fine_index}'].reshape(fine_grid,fine_grid,2)
    fine_axis=np.arange(fine_grid)*length/fine_grid
    fx,fy=np.meshgrid(origin[0]+fine_axis,origin[1]+fine_axis)
    fine_excluded=window(fx-cg[0],inner[0]-.5*beam,inner[0])*window(fy-cg[1],inner[1]-.5*beam,inner[1])
    fine_far=phases(fine_trace.reshape(-1,2)).reshape(fine_grid,fine_grid,4)*(1-fine_excluded[:,:,None])
    fine_spectrum=np.fft.fftshift(np.fft.fft2(fine_far,axes=(0,1)),axes=(0,1))
    first=(fine_grid-grid)//2
    projected=np.fft.ifftshift(fine_spectrum[first:first+grid,first:first+grid],axes=(0,1))*(grid/fine_grid)**2
    projected_trace=amplitude((pprobe@projected.reshape(-1,4)).real)
    source_vertical=2*checked_influence(targets,panels).gradient_single_layer[:,:,2]
    base_vertical=amplitude(vertical_base+source_vertical@base_density)
    # A diagnostic BELOW the interface separates a rough P1 boundary derivative
    # from an inaccurate harmonic field. It does not qualify surface evolution.
    submerged=targets.copy();submerged[:,2]=-.1*beam
    solid=ConvexHull(triangles.reshape(-1,3))
    if np.any(np.all(submerged@solid.equations[:,:3].T+solid.equations[:,3]<=0,axis=1)):
        raise ValueError('subsurface diagnostic targets entered the physical solid')
    reflected=submerged.copy();reflected[:,2]*=-1
    bulk_source=(checked_influence(submerged,panels).gradient_single_layer[:,:,2]
                 +checked_influence(reflected,panels).gradient_single_layer[:,:,2])
    _,nbulk,_=harmonic_maps(submerged,np.tile([0.,0.,1.],(len(targets),1)),grid,length/grid,origin)
    bulk_base=(nbulk@spectrum.reshape(-1,4)).real
    bulk_far=(nbulk@far_spectrum.reshape(-1,4)).real
    reference_bulk=field_at(submerged,velocity=True)[:,:,2]
    base_bulk=amplitude(bulk_base+bulk_source@base_density)
    correction=[];normal_correction=[];vertical_correction=[];trace_correction=[];bulk_correction=[]
    before=perf_counter()
    for patch in patches:
        value,velocity=patch.bulk(panels.centroids)
        correction.append(value);normal_correction.append(np.sum(velocity*panels.normals,axis=1))
        vertical_correction.append(np.r_[patch.surface_velocity(selected)[:,2],patch.outside_velocity(far_xy)[:,2]])
        trace_correction.append(np.r_[patch.values.ravel()[patch.faces[selected]].mean(axis=1),np.zeros(len(far_xy))])
        bulk_correction.append(patch.bulk(submerged)[1][:,2])
    correction=np.array(correction).T;normal_correction=np.array(normal_correction).T
    vertical_correction=np.array(vertical_correction).T;trace_correction=np.array(trace_correction).T
    density=lu_solve(factor,q-far_normal-normal_correction)
    enriched_body=amplitude(far_lift+correction+single@density)
    enriched_vertical=amplitude(vertical_far+vertical_correction+source_vertical@density)
    enriched_bulk=amplitude(bulk_far+np.array(bulk_correction).T+bulk_source@density)
    report.update(qualified=False,scope=__doc__,input_sha256=identity,
        panel_sha256=audit['panel_sha256'],grid=grid,length_m=length,body_panels=len(panels),
        patch_spacing_m=[float(np.diff(x).max()),float(np.diff(y).max())],patch_triangles=len(patches[0].panels),
        patch_extent_m=(2*outer).tolist(),near_probes=len(probes),far_probes=len(far_xy),
        probe_sha256=sha256(targets.tobytes()).hexdigest(),
        probe_matching_max_distance_m=float(np.max(match_distance)),
        coarse_trace_ownership='Smooth near-hull trace detail excluded before coarse sampling; compact local residual restores it on a larger patch',
        local_application_seconds=perf_counter()-before,preparation_plus_application_seconds=perf_counter()-started,
        reference_robin_relative_l2=(norms(reference_vertical-omega**2/GRAVITY*reference_trace,np.ones(len(targets)))
            /norms(reference_vertical,np.ones(len(targets)))).tolist(),
        far_projection_audit=dict(source_grid=fine_grid,target_grid=grid,
            relative_l2=(norms(projected_trace[len(probes):]-reference_trace[len(probes):],np.ones(len(far_xy)))
                /norms(reference_trace[len(probes):],np.ones(len(far_xy)))).tolist(),
            scope='Coarse-far Dirichlet sampling check only; projected coefficients do not enter this body solve'),rows=[])
    for label,value,velocity,surface_trace,bulk in [
        ('coarse global',base_body,base_vertical,amplitude(trace_base),base_bulk),
        ('coarse plus compact local trace',enriched_body,enriched_vertical,amplitude(trace_far+trace_correction),enriched_bulk)]:
        force=1j*WATER_DENSITY*omega*weighted.T@value
        row=dict(method=label,body_potential_relative_l2=(norms(value-reference_body,panels.areas)
            /norms(reference_body,panels.areas)).tolist(),
            principal_added_mass_from_snapshot=(force.imag/omega)[[2,3],[0,1]].tolist(),
            six_dof_force_amplitude_error_real=(force-reference_force).real.tolist(),
            six_dof_force_amplitude_error_imag=(force-reference_force).imag.tolist(),zones=[])
        for name,ids in [('near',np.arange(len(probes))),('far',np.arange(len(probes),len(targets)))]:
            w=np.ones(len(ids))
            row['zones'].append(dict(name=name,
                trace_relative_l2=(norms(surface_trace[ids]-reference_trace[ids],w)/norms(reference_trace[ids],w)).tolist(),
                vertical_velocity_relative_l2=(norms(velocity[ids]-reference_vertical[ids],w)/norms(reference_vertical[ids],w)).tolist(),
                subsurface_depth_m=.1*beam,
                subsurface_vertical_velocity_relative_l2=(norms(bulk[ids]-reference_bulk[ids],w)/norms(reference_bulk[ids],w)).tolist(),
                subsurface_vertical_velocity_max_absolute_error=np.max(abs(bulk[ids]-reference_bulk[ids]),axis=0).tolist(),
                subsurface_vertical_velocity_rms_absolute_error=np.sqrt(np.mean(abs(bulk[ids]-reference_bulk[ids])**2,axis=0)).tolist(),
                reference_subsurface_vertical_velocity_rms=np.sqrt(np.mean(abs(reference_bulk[ids])**2,axis=0)).tolist(),
                vertical_velocity_rms_absolute_error=np.sqrt(np.mean(abs(velocity[ids]-reference_vertical[ids])**2,axis=0)).tolist(),
                reference_vertical_velocity_rms=np.sqrt(np.mean(abs(reference_vertical[ids])**2,axis=0)).tolist(),
                vertical_velocity_max_absolute_error=np.max(abs(velocity[ids]-reference_vertical[ids]),axis=0).tolist()))
        report['rows'].append(row);print(json.dumps(row),flush=True)
    np.savez(args.output.with_suffix('.npz'),targets=targets,base_body=base_body,enriched_body=enriched_body,
        reference_body=reference_body,base_vertical=base_vertical,enriched_vertical=enriched_vertical,
        reference_vertical=reference_vertical,reference_trace=reference_trace,
        local_vertical=amplitude(vertical_correction),patch_x=x,patch_y=y,residual=residual,
        base_bulk=base_bulk,enriched_bulk=enriched_bulk,reference_bulk=reference_bulk,
        base_density=base_density,enriched_density=density)
    report.update(complete=True,arrays_sha256=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest(),
        peak_process_rss_bytes=peak_rss_bytes(),limits=[
            'Frozen external field derivative, not exact Robin/open-water physics; default reference panel integration defect retained',
            'Full reference auxiliary trace, no physical hidden-state ownership or moving geometry',
            'P1 centroid DNO samples are not nodal/edge evolution rates',
            'No nonlinear trajectory, minute accuracy or simulation throughput qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field',type=Path,required=True)
    parser.add_argument('--audit',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--local-refinement',type=int,choices=[0,1,2],default=0)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output/source archive required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    shutil.copyfile(args.field,archive/'parent-field.json')
    shutil.copyfile(args.audit,archive/'parent-audit.json')
    shutil.copyfile(args.audit.with_suffix('.source')/'inputs.json',archive/'inputs.json')
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,
        parent_field_sha256=sha256(args.field.read_bytes()).hexdigest(),parent_audit_sha256=sha256(args.audit.read_bytes()).hexdigest())
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
