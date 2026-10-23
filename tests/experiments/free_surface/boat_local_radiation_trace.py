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
from scipy.spatial import ConvexHull
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
    cells=np.ceil(2*outer/args.spacing).astype(int)
    x=np.linspace(cg[0]-outer[0],cg[0]+outer[0],cells[0]+1)
    y=np.linspace(cg[1]-outer[1],cg[1]+outer[1],cells[1]+1)
    xx,yy=np.meshgrid(x,y)
    nodes=np.column_stack((xx.ravel(),yy.ravel(),np.zeros(xx.size)))
    exact_trace=phases(field_at(nodes)).reshape(len(y),len(x),4)
    coarse_trace=tensor_values(spectrum,x,y,origin,length)
    taper=window(xx-cg[0],inner[0],outer[0])*window(yy-cg[1],inner[1],outer[1])
    residual=(exact_trace-coarse_trace)*taper[:,:,None]
    residual[[0,-1],:]=0.;residual[:,[0,-1]]=0.
    patches=[CompactP1Trace(x,y,residual[:,:,c]) for c in range(4)]
    centroids=patches[0].panels.centroids
    footprint=ConvexHull(wet.waterline_segments[:,:,:2].reshape(-1,2))
    distance=np.max(centroids[:,:2]@footprint.equations[:,:2].T+footprint.equations[:,2],axis=1)
    candidates=np.flatnonzero((distance>=.03*beam)&(distance<=.5*beam)
        &np.all(abs(centroids[:,:2]-cg[:2])<inner,axis=1))
    if not len(candidates):
        raise ValueError('finite-clearance exterior centroid probes required')
    selected=candidates[np.linspace(0,len(candidates)-1,min(96,len(candidates)),dtype=int)]
    probes=centroids[selected]
    angle=2*np.pi*(np.arange(24)+.37)/24
    far_xy=cg[:2]+2*max(outer)*np.column_stack((np.cos(angle),np.sin(angle)))
    targets=np.r_[probes,np.column_stack((far_xy,np.zeros(len(far_xy))))]
    reference_trace=field_at(targets);reference_vertical=field_at(targets,velocity=True)[:,:,2]
    pprobe,nprobe,_=harmonic_maps(targets,np.tile([0.,0.,1.],(len(targets),1)),grid,length/grid,origin)
    trace_base=(pprobe@spectrum.reshape(-1,4)).real
    vertical_base=(nprobe@spectrum.reshape(-1,4)).real
    source_vertical=2*checked_influence(targets,panels).gradient_single_layer[:,:,2]
    base_vertical=amplitude(vertical_base+source_vertical@base_density)
    correction=[];normal_correction=[];vertical_correction=[];trace_correction=[]
    before=perf_counter()
    for patch in patches:
        value,velocity=patch.bulk(panels.centroids)
        correction.append(value);normal_correction.append(np.sum(velocity*panels.normals,axis=1))
        vertical_correction.append(np.r_[patch.surface_velocity(selected)[:,2],patch.outside_velocity(far_xy)[:,2]])
        trace_correction.append(np.r_[patch.values.ravel()[patch.faces[selected]].mean(axis=1),np.zeros(len(far_xy))])
    correction=np.array(correction).T;normal_correction=np.array(normal_correction).T
    vertical_correction=np.array(vertical_correction).T;trace_correction=np.array(trace_correction).T
    density=lu_solve(factor,q-normal-normal_correction)
    enriched_body=amplitude(lift+correction+single@density)
    enriched_vertical=amplitude(vertical_base+vertical_correction+source_vertical@density)
    report.update(qualified=False,scope=__doc__,input_sha256=identity,
        panel_sha256=audit['panel_sha256'],grid=grid,length_m=length,body_panels=len(panels),
        patch_spacing_m=[float(np.diff(x).max()),float(np.diff(y).max())],patch_triangles=len(patches[0].panels),
        patch_extent_m=(2*outer).tolist(),near_probes=len(probes),far_probes=len(far_xy),
        local_application_seconds=perf_counter()-before,preparation_plus_application_seconds=perf_counter()-started,
        reference_robin_relative_l2=(norms(reference_vertical-omega**2/GRAVITY*reference_trace,np.ones(len(targets)))
            /norms(reference_vertical,np.ones(len(targets)))).tolist(),rows=[])
    for label,value,velocity,surface_trace in [
        ('coarse global',base_body,base_vertical,amplitude(trace_base)),
        ('coarse plus compact local trace',enriched_body,enriched_vertical,amplitude(trace_base+trace_correction))]:
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
                vertical_velocity_max_absolute_error=np.max(abs(velocity[ids]-reference_vertical[ids]),axis=0).tolist()))
        report['rows'].append(row);print(json.dumps(row),flush=True)
    np.savez(args.output.with_suffix('.npz'),targets=targets,base_body=base_body,enriched_body=enriched_body,
        reference_body=reference_body,base_vertical=base_vertical,enriched_vertical=enriched_vertical,
        reference_vertical=reference_vertical,reference_trace=reference_trace,
        local_vertical=amplitude(vertical_correction),patch_x=x,patch_y=y,residual=residual)
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
    parser.add_argument('--spacing',type=float,default=.25)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or not np.isfinite(args.spacing) or args.spacing<=0:
        parser.error('fresh output/source archive and positive spacing required')
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
