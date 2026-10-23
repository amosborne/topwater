"""Actual-hull near/far TRACE enrichment against an independent harmonic field.

The scalar is a velocity potential: gradients are velocities and rho*int(phi*n)
are generalized impulses, not forces. The same linear map acts on pressure
potential, but no nonlinear pressure or trajectory is qualified by this test.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.spatial import ConvexHull

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import GRAVITY,WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.benchmark_rankine import peak_rss_bytes
from .body_internal_field import internal_field
from .boat_hos_continuation import window
from .boat_radiation_time import harmonic_maps
from .boat_spatial_work import source_matrix,subdivide
from .compact_p1_trace import CompactP1Trace
from .panel_gradient_reference import checked_influence


def norms(value,weights):
    return np.sqrt(np.sum(weights[:,None]*abs(value)**2,axis=0))


def tensor_values(spectrum,x,y,origin,length):
    """Exact trigonometric interpolation on tensor nodes, no dense point map."""
    size=spectrum.shape[0]
    k=2*np.pi*np.fft.fftfreq(size,d=length/size)
    ex=np.exp(1j*(x-origin[0])[:,None]*k)
    ey=np.exp(1j*(y-origin[1])[:,None]*k)
    return np.stack([(ey@spectrum[:,:,c]@ex.T).real/(size*size)
                     for c in range(spectrum.shape[2])],axis=-1)


def run(inputs,args,report):
    started=perf_counter()
    wet=hydrostatics(inputs.body,inputs.initial_state,density=WATER_DENSITY,gravity=GRAVITY)
    fixture=internal_field(wet.wetted_triangles) # Freeze before discretization changes.
    triangles=wet.wetted_triangles
    for _ in range(args.refinement):
        triangles=subdivide(triangles)
    panels=TrianglePanels(triangles,-1)
    cg=inputs.initial_state.position
    beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    weighted=panels.areas[:,None]*modes
    exact_body,exact_gradient=fixture.field(panels.centroids)
    q=np.einsum('ncj,nj->nc',exact_gradient,panels.normals)
    exact_centroid_impulse=WATER_DENSITY*weighted.T@exact_body
    exact_integrated_impulse=-fixture.pressure_load(panels,cg,WATER_DENSITY)

    def image_field(points):
        potential,gradient=fixture.field(points*[1.,1.,-1.])
        return potential,gradient*[1.,1.,-1.]

    size,length=args.grid,args.length
    origin=cg[:2]-length/2
    axis=np.arange(size)*length/size
    xx,yy=np.meshgrid(origin[0]+axis,origin[1]+axis)
    surface=np.column_stack((xx.ravel(),yy.ravel(),np.zeros(xx.size)))
    trace=fixture.field(surface)[0].reshape(size,size,2)
    spectrum=np.fft.fft2(trace,axes=(0,1))
    pbody,nbody,_=harmonic_maps(panels.centroids,panels.normals,size,length/size,origin)
    coarse_body=(pbody@spectrum.reshape(-1,2)).real
    coarse_normal=(nbody@spectrum.reshape(-1,2)).real
    exact_lift,exact_lift_gradient=image_field(panels.centroids)
    exact_normal=np.einsum('ncj,nj->nc',exact_lift_gradient,panels.normals)
    factor=lu_factor(source_matrix(panels))
    single=assemble_dirichlet_halfspace(panels).single_layer
    oracle_density=lu_solve(factor,q-exact_normal)
    oracle_body=exact_lift+single@oracle_density
    base_density=lu_solve(factor,q-coarse_normal)
    base_body=coarse_body+single@base_density

    inner=np.max(abs(wet.wetted_triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam
    outer=inner+beam
    if np.any(outer>=length/2):
        raise ValueError('compact patch must fit strictly inside the coarse box')
    cells=np.ceil(2*outer/args.spacing).astype(int)
    x=np.linspace(cg[0]-outer[0],cg[0]+outer[0],cells[0]+1)
    y=np.linspace(cg[1]-outer[1],cg[1]+outer[1],cells[1]+1)
    px,py=np.meshgrid(x,y)
    patch_nodes=np.column_stack((px.ravel(),py.ravel(),np.zeros(px.size)))
    exact_trace=fixture.field(patch_nodes)[0].reshape(len(y),len(x),2)
    coarse_trace=tensor_values(spectrum,x,y,origin,length)
    taper=window(px-cg[0],inner[0],outer[0])*window(py-cg[1],inner[1],outer[1])
    residual=(exact_trace-coarse_trace)*taper[:,:,None]
    residual[[0,-1],:]=0.;residual[:,[0,-1]]=0.
    patches=[CompactP1Trace(x,y,residual[:,:,c]) for c in range(2)]
    centroids=patches[0].panels.centroids
    footprint=ConvexHull(wet.waterline_segments[:,:,:2].reshape(-1,2))
    outside_distance=np.max(centroids[:,:2]@footprint.equations[:,:2].T+footprint.equations[:,2],axis=1)
    on_plateau=np.all(abs(centroids[:,:2]-cg[:2])<inner,axis=1)
    candidates=np.flatnonzero((outside_distance>=.03*beam)&(outside_distance<=.5*beam)&on_plateau)
    if not len(candidates):
        raise ValueError('patch must supply finite-clearance exterior probes')
    selected=candidates[np.linspace(0,len(candidates)-1,min(96,len(candidates)),dtype=int)]
    probes=centroids[selected]
    angles=2*np.pi*(np.arange(24)+.37)/24
    far_xy=cg[:2]+2*max(outer)*np.column_stack((np.cos(angles),np.sin(angles)))
    far=np.column_stack((far_xy,np.zeros(len(far_xy))))
    all_probes=np.concatenate((probes,far))
    pprobe,nprobe,_=harmonic_maps(all_probes,np.tile([0.,0.,1.],(len(all_probes),1)),
                                size,length/size,origin)
    coarse_probe=(pprobe@spectrum.reshape(-1,2)).real
    coarse_vertical=(nprobe@spectrum.reshape(-1,2)).real
    image_probe,image_gradient=image_field(all_probes)
    exact_probe,probe_gradient=fixture.field(all_probes)
    exact_vertical=probe_gradient[:,:,2]
    source_vertical=np.empty((len(all_probes),len(panels)))
    for first in range(0,len(all_probes),128):
        last=min(first+128,len(all_probes))
        source_vertical[first:last]=2*checked_influence(all_probes[first:last],panels).gradient_single_layer[:,:,2]
    base_vertical=coarse_vertical+source_vertical@base_density
    oracle_vertical=image_gradient[:,:,2]+source_vertical@oracle_density
    report.update(input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        panel_sha256=sha256(triangles.tobytes()).hexdigest(),body_panels=len(panels),
        grid=size,length_m=length,patch_cells=cells.tolist(),patch_extent_m=(2*outer).tolist(),
        patch_spacing_m=[float(np.diff(x).max()),float(np.diff(y).max())],
        patch_triangles=len(patches[0].panels),near_probes=len(probes),far_probes=len(far),
        fixture=dict(points=fixture.points.tolist(),strengths=fixture.strengths.tolist(),
            center=fixture.center.tolist(),inscribed_radius=fixture.inscribed_radius,
            scalar='velocity potential, unit source strengths in m^3/s'),
        interpretation='Original analytic field is physical-field reference; reflected sources supply exact harmonic lift; same-J exact-lift oracle separates surface error from body discretization',
        potential_quadrature_impulse_error=(exact_centroid_impulse-exact_integrated_impulse).tolist(),
        exact_integrated_generalized_impulse=exact_integrated_impulse.tolist(),
        setup_seconds=perf_counter()-started,rows=[])
    saved=dict(panels=triangles,probe_points=all_probes,exact_body=exact_body,oracle_body=oracle_body,
        exact_vertical=exact_vertical,oracle_vertical=oracle_vertical,base_body=base_body,
        base_vertical=base_vertical,patch_x=x,patch_y=y,residual=residual)
    rows=[('coarse global',base_body,base_vertical,coarse_body,coarse_normal,coarse_probe)]
    before=perf_counter()
    correction=[];correction_normal=[];correction_vertical=[];correction_probe=[]
    for c,patch in enumerate(patches):
        value,velocity=patch.bulk(panels.centroids)
        correction.append(value)
        correction_normal.append(np.einsum('nj,nj->n',velocity,panels.normals))
        correction_vertical.append(np.r_[patch.surface_velocity(selected)[:,2],patch.outside_velocity(far_xy)[:,2]])
        correction_probe.append(np.r_[patch.values.ravel()[patch.faces[selected]].mean(axis=1),np.zeros(len(far))])
    correction=np.array(correction).T;correction_normal=np.array(correction_normal).T
    correction_vertical=np.array(correction_vertical).T;correction_probe=np.array(correction_probe).T
    density=lu_solve(factor,q-coarse_normal-correction_normal)
    body=coarse_body+correction+single@density
    vertical=coarse_vertical+correction_vertical+source_vertical@density
    report['local_application_seconds']=perf_counter()-before
    rows.append(('coarse plus compact local trace',body,vertical,coarse_body+correction,
                 coarse_normal+correction_normal,coarse_probe+correction_probe))
    rows.append(('same body operator exact image lift',oracle_body,oracle_vertical,
                 exact_lift,exact_normal,image_probe))
    for label,body,vertical,lift,normal,surface_trace in rows:
        impulse=WATER_DENSITY*weighted.T@body
        row=dict(method=label,
            body_potential_relative_l2=(norms(body-exact_body,panels.areas)/norms(exact_body,panels.areas)).tolist(),
            body_potential_vs_lift_oracle_relative_l2=(norms(body-oracle_body,panels.areas)/norms(oracle_body,panels.areas)).tolist(),
            body_lift_potential_relative_l2=(norms(lift-exact_lift,panels.areas)/norms(exact_lift,panels.areas)).tolist(),
            body_lift_normal_relative_l2=(norms(normal-exact_normal,panels.areas)/norms(exact_normal,panels.areas)).tolist(),
            generalized_impulse=impulse.tolist(),
            impulse_centroid_field_error=(impulse-exact_centroid_impulse).tolist(),
            impulse_total_error=(impulse-exact_integrated_impulse).tolist(),
            zones=[])
        for name,indices in [('near',np.arange(len(probes))),('far',np.arange(len(probes),len(all_probes)))]:
            w=np.ones(len(indices))
            row['zones'].append(dict(name=name,
                dirichlet_relative_l2=(norms(surface_trace[indices]-exact_probe[indices],w)/norms(exact_probe[indices],w)).tolist(),
                normal_velocity_relative_l2=(norms(vertical[indices]-exact_vertical[indices],w)/norms(exact_vertical[indices],w)).tolist(),
                normal_velocity_vs_lift_oracle_relative_l2=(norms(vertical[indices]-oracle_vertical[indices],w)/norms(oracle_vertical[indices],w)).tolist(),
                normal_velocity_max_absolute_error_m_per_s=np.max(abs(vertical[indices]-exact_vertical[indices]),axis=0).tolist()))
        report['rows'].append(row)
        print(json.dumps(row),flush=True)
    saved.update(enriched_body=rows[1][1],enriched_vertical=rows[1][2],
                 local_vertical=correction_vertical,local_body_potential=correction)
    np.savez(args.output.with_suffix('.npz'),**saved)
    report.update(complete=True,qualified=False,scope=__doc__,peak_process_rss_bytes=peak_rss_bytes(),
        arrays_sha256=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest(),
        limits=['Manufactured body-near-field test, not a supported incident-wave mesh requirement',
            'Full known auxiliary trace, no physical hidden-state transfer or eta/psi evolution',
            'P1 edge/nodal DNO singularities; owned centroid probes are not nodal evolution rates',
            'Periodic coarse plus aperiodic local lift; exact image oracle exposes domain/representation error',
            'No moving nonlinear case, pressure time derivative, full-minute accuracy or throughput qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--grid',type=int,default=129)
    parser.add_argument('--length',type=float,default=128.)
    parser.add_argument('--spacing',type=float,default=.25)
    parser.add_argument('--refinement',type=int,default=1)
    args=parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or args.grid<17 or args.grid%2!=1 or args.refinement not in (0,1,2)
            or not np.isfinite([args.length,args.spacing]).all() or min(args.length,args.spacing)<=0):
        parser.error('fresh output, odd grid>=17, positive length/spacing and refinement0..2 required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs,archive/'inputs.json')
    inputs=SimulationInputs.from_json(args.inputs.read_text())
    report=dict(complete=False,producer_sha256=hashes)
    started=perf_counter()
    try:
        run(inputs,args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error))
        raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
