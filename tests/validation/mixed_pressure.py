"""Bounded calm pressure component and nonflat boundary closure; no TD speed pass."""
from time import perf_counter
STARTED=perf_counter()
import argparse
import json
from pathlib import Path
import shutil
from zipfile import ZipFile
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs,WaveComponents
from topwater.rankine import TrianglePanels
from topwater.bem import assemble_laplace
from topwater.boundary import BoundarySystem
from topwater.hydrostatics import hydrostatics
from topwater.fluid_boundary import closed_boundary
from topwater.pressure import PressureSystem,pressure_acceleration
from .common import prepare,save,read_owned,digest,current_sources
from .mild import initial_body
from .periodic_grid import mean_angles,radius,annulus


def generate(args,c,report,out):
    if (c['format_version']!=1 or c['calm_angular_divisions']!=32 or c['outer_radius_m']!=24.
            or c['depth_m']!=24. or c['radial_offsets_m']!=[0.,.15,.4,1.,2.5,6.,12.]
            or c['density_kg_m3']!=1025. or c['gravity_m_s2']!=9.81
            or c['maximum_elapsed_seconds']!=180 or c['cpu_threads']!=1
            or c['pressure_evaluations_per_simulated_second']!=200
            or c['runtime_limit_wall_seconds_per_simulated_second']!=200.):
        raise ValueError('unsupported prospective pressure controls')
    fields,a=read_owned(args.fields,digest('validation/periodic-grid-controls.json'),'periodic-field-map-candidate')
    current_sources(fields,tuple(fields['producer_sha256']))
    source,s=read_owned(args.sources,digest('validation/periodic-source-controls.json'),'periodic-source-recovery',passed=True)
    if (digest(args.bundle)!=c['field_bundle_sha256'] or fields['source_manifest_sha256']!=digest(args.sources)
            or fields['completed_component_evaluations']!=12 or source['completed_physics_calls']!=84):
        raise ValueError('missing complete identity-matched parent evidence')
    archive=args.output.with_suffix('.source')
    with ZipFile(args.bundle) as z:
        for suffix in ('.json','.npz'):
            if z.read('validation/'+args.fields.stem+suffix)!=args.fields.with_suffix(suffix).read_bytes():
                raise ValueError('field bytes differ from pinned bundle')
        physical=z.read('validation/'+args.fields.stem+'.source/frozen-inputs.json')
    shutil.copyfile(args.bundle,archive/'field-map.zip');shutil.copyfile(c['protocol'],archive/'protocol.md')
    (archive/'frozen-inputs.json').write_bytes(physical)
    mean=SimulationInputs.from_json(physical.decode('utf-8'))
    if np.any(mean.initial_state.linear_velocity) or np.any(mean.initial_state.angular_velocity):raise ValueError('calm probe requires actual rest')
    calm=SimulationInputs(mean.body,WaveComponents([],[],[],[],provenance={'scope':'calm pressure component only'}),mean.initial_state)
    (archive/'calm-inputs.json').write_text(calm.to_json())
    phasors=np.array([s['raw_RAO'][-1,f,f]*s['incident_phasors_m'][f] for f in range(3)])
    displacement=phasors.sum(axis=0).real;velocity=(-1j*s['angular_frequencies_rad_s'][:,None]*phasors).sum(axis=0).real
    mild_state=initial_body(mean.initial_state,displacement,velocity)
    # Physical initial pose/velocity only. No future trajectory or acceleration.
    (archive/'mild-body-state.json').write_text(json.dumps({name:getattr(mild_state,name).tolist() for name in ('position','orientation','linear_velocity','angular_velocity')},indent=2)+'\n')
    tetra=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    tetra_faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
    out.update(warmup_vertices=tetra,warmup_faces=tetra_faces)
    report.update(field_manifest_sha256=digest(args.fields),source_manifest_sha256=digest(args.sources),
        field_bundle_sha256=digest(args.bundle),protocol_sha256=digest(c['protocol']),
        calm_inputs_sha256=digest(archive/'calm-inputs.json'),mild_initial_state_sha256=digest(archive/'mild-body-state.json'),
        frozen_inputs_sha256=digest(archive/'frozen-inputs.json'),
        full_reference_qualified=False,nonlinear_rhs_ready=False,complete_runtime_speed_qualified=False,
        original_raw_passed=fields['original_raw_passed'],scope=__doc__,stages={})
    save(args.output,report,out)
    def budget():
        if perf_counter()-STARTED>c['maximum_elapsed_seconds']:raise TimeoutError('bounded pressure probe exceeded elapsed budget')
    with threadpool_limits(1):
        report['native_thread_pools']=threadpool_info()
        if any(p['num_threads']!=1 for p in threadpool_info()):raise ValueError('one native thread required')
        budget()
        # Load native/JIT paths once. Account separately; do not multiply a cold
        # compilation/cache-load cost by the number of future time steps.
        begin=perf_counter();p=TrianglePanels(tetra[tetra_faces]);op=assemble_laplace(p)
        warm=BoundarySystem(op,dirichlet=np.array([True,False,False,False]),areas=p.areas,length_scale=1.,bounded=True)
        warm.solve(potential=np.zeros(4),normal_velocity=np.zeros(4))
        report['stages']['normal_numerical_warmup_seconds']=perf_counter()-begin;budget();save(args.output,report,out)
        begin=perf_counter();loads=hydrostatics(mean.body,mean.initial_state,density=c['density_kg_m3'],gravity=c['gravity_m_s2'])
        origin=mean.initial_state.position[:2];angles=mean_angles(loads.wetted_triangles,origin,c['calm_angular_divisions'])
        inner=radius(a['mean_hull_planes'],angles,0.,origin)
        vertices,faces,r=annulus(angles,inner,c['outer_radius_m'],c['radial_offsets_m'],origin);n=len(angles)
        calm_boundary=closed_boundary(mean.body,mean.initial_state,vertices,faces,contact=np.arange(n),
            outer=np.arange(len(vertices)-n,len(vertices)),center_xy=origin,depth=c['depth_m'])
        report['stages']['calm_geometry_seconds']=perf_counter()-begin;report['calm_geometry']=calm_boundary.diagnostics
        out.update(calm_vertices=calm_boundary.mesh.vertices,calm_faces=calm_boundary.mesh.faces,
            calm_tags=calm_boundary.tags,calm_surface_vertices=vertices,calm_surface_faces=faces)
        budget();save(args.output,report,out)
        begin=perf_counter();owner=PressureSystem(calm_boundary,calm_boundary.tags=='surface',length_scale=c['outer_radius_m'])
        report['stages']['calm_assembly_factor_seconds']=perf_counter()-begin
        report['laplace']=owner.system.operators.diagnostics;report['factorization']=owner.system.diagnostics
        budget();save(args.output,report,out)
        begin=perf_counter();count=len(owner.panels.areas)
        result=pressure_acceleration(owner,known_potential=np.zeros(count),known_normal_acceleration=np.zeros(count),
            body_speed_squared=np.zeros(owner.body_mask.sum()),density=c['density_kg_m3'],gravity=c['gravity_m_s2'])
        report['stages']['calm_pressure_seconds']=perf_counter()-begin
        out.update(calm_acceleration_world=result.acceleration_world,calm_added_mass=result.added_mass,
            calm_chi=result.acceleration_potential,calm_fluid_load=result.fluid_load,calm_external_load=result.external_load,
            raw_reference_Ainf=s['raw_Ainf'],reference_Ainf_increment_envelope=np.maximum(abs(s['raw_Ainf'][-2]-s['raw_Ainf'][-3]),abs(s['raw_Ainf'][-1]-s['raw_Ainf'][-2])))
        report['pressure']=result.diagnostics
        if np.max(abs(result.acceleration_world))>c['calm_acceleration_absolute_SI']:raise ValueError('exact canonical calm equilibrium failed')
        if result.diagnostics['balance_relative']>c['pressure_body_balance_scaled_tolerance']:raise ValueError('pressure/body balance tolerance failed')
        budget();save(args.output,report,out)
        begin=perf_counter();v=a['scale_1.0_vertices'];f=a['scale_1.0_faces'];n=len(a['scale_1.0_angles_rad'])
        mild_boundary=closed_boundary(mean.body,mild_state,v,f,contact=np.arange(n),
            outer=np.arange(len(v)-n,len(v)),center_xy=origin,depth=c['depth_m'])
        report['stages']['mild_geometry_only_seconds']=perf_counter()-begin;report['mild_geometry']=mild_boundary.diagnostics
        out.update(mild_vertices=mild_boundary.mesh.vertices,mild_faces=mild_boundary.mesh.faces,mild_tags=mild_boundary.tags)
        budget()
    per_evaluation=sum(report['stages'][k] for k in ('calm_geometry_seconds','calm_assembly_factor_seconds','calm_pressure_seconds'))
    projected=report['stages']['normal_numerical_warmup_seconds']+c['pressure_evaluations_per_simulated_second']*per_evaluation
    report.update(complete=True,passed=False,elapsed_seconds=perf_counter()-STARTED,
        pressure_component_seconds_per_evaluation=per_evaluation,pressure_component_projection_wall_seconds_per_simulated_second=projected,
        provisional_midpoint_route_rejected_by_component_cost=projected>c['runtime_limit_wall_seconds_per_simulated_second'],
        cost_scope='component estimate on calm mesh only; no full RHS or mild-runtime cost pass',
        qualification_status='unqualified: derivatives/nonlinear RHS, physical/reference convergence, domain and time-domain comparison remain missing')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,default=Path('validation/mixed-pressure-controls.json'))
    for name in ('fields','sources','bundle','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'mixed-pressure-component');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a)
    print(json.dumps({k:r.get(k) for k in ('complete','passed','stages','calm_geometry','mild_geometry','pressure_component_projection_wall_seconds_per_simulated_second','error')}))
    return 0 if r['complete'] and r['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
