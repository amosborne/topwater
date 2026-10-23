"""Bounded evolving-water falsifier with prescribed radiation/far-field ports.

Physical eta and projected physical psi evolve through the actual weak
kinematic/dynamic equations. Fixed hull, prescribed body motion and prescribed
far chi remain verification restrictions: this is NOT a nonlinear6DOF run.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.integrate import RK23
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import splu

from topwater import SimulationInputs
from topwater.constants import GRAVITY,WATER_DENSITY
from topwater.rankine import TrianglePanels,panel_influence
from tests.benchmark_rankine import peak_rss_bytes
from .boat_benchmark_contract import measurement_environment
from .boat_integration_controls import (RELATIVE_TOLERANCE,SURFACE_HEIGHT_ABSOLUTE_TOLERANCE_M,
    SURFACE_POTENTIAL_ABSOLUTE_TOLERANCE_M2_PER_S)
from .boat_local_radiation_trace import phases,amplitude,far_trace_spectrum
from .boat_local_trace_enrichment import tensor_values
from .boat_regularized_phase_rhs import cut_quadrature,scatter_moments
from .boat_weak_contact_rows import owned,fourier_values
from .compact_p1_weak import WeakPatch,hat_fourier
from .physical_p1_patch import PhysicalP1Patch
from .physical_weak_flux import PhysicalWeakFlux
from .cuda_operator import cp


def run(args,report):
    phase=owned(args.phase);order=phase['orders'][-1]['order']
    # Reusing maps also reuses their hidden-extension and source-basis owner.
    # Matching physical mass alone would not detect a changed auxiliary fill.
    for name in ('tests/experiments/free_surface/physical_p1_patch.py',
                 'tests/experiments/free_surface/compact_p1_trace.py',
                 'tests/experiments/free_surface/compact_p1_maps.py',
                 'tests/experiments/free_surface/panel_gradient_reference.py',
                 'src/topwater/rankine.py','src/topwater/_rankine_kernel.py',
                 'src/topwater/_rankine_extended.py','src/topwater/numerics.py','src/topwater/constants.py'):
        path=Path(name)
        if sha256(path.read_bytes()).hexdigest()!=phase['producer_sha256'].get(str(path)):
            raise ValueError('cached phase map owner changed; rebuild phase operators: '+name)
    contact_path=args.phase.with_suffix('.source')/'contact.json'
    if sha256(contact_path.read_bytes()).hexdigest()!=phase['parent_contact_sha256']:
        raise ValueError('phase/contact ownership changed')
    contact=json.loads(contact_path.read_text())
    paths={name:Path(record['path']) for name,record in contact['parents'].items()}
    parents={}
    for name,path in paths.items():
        if sha256(path.read_bytes()).hexdigest()!=contact['parents'][name]['sha256']:
            raise ValueError('contact parent ownership changed: '+name)
        parents[name]=owned(path)
    inputs=SimulationInputs.from_json((paths['audit'].with_suffix('.source')/'inputs.json').read_text())
    if sha256(inputs.to_json().encode()).hexdigest()!=phase['input_sha256']:
        raise ValueError('physical input ownership changed')
    with np.load(args.phase.with_suffix('.npz')) as z:
        x=z['x'].copy();y=z['y'].copy();polygon=z['polygon'].copy()
        pmap=z['pmap'].copy();nmap=z['nmap'].copy();even=z['even'].copy();neumann=z['neumann'].copy()
        far_body=z['far_body'].copy();far_normal=z['far_normal'].copy();cross=z[f'q{order}_cross'].copy()
        psi_mom=z[f'q{order}_psi_mom'].copy();initial_c=z[f'q{order}_coefficients'].copy()
        initial_sigma=z[f'q{order}_sigma'].copy();old_force=z[f'q{order}_force'].copy()
        old_active=z['active_nodes'].copy();old_mass=(z['mass_data'].copy(),z['mass_indices'].copy(),z['mass_indptr'].copy())
    patch=PhysicalP1Patch(x,y,polygon);count=len(patch.active_nodes)
    np.testing.assert_array_equal(patch.active_nodes,old_active)
    for current,previous in zip((patch.mass.data,patch.mass.indices,patch.mass.indptr),old_mass,strict=True):
        np.testing.assert_array_equal(current,previous)
    with np.load(paths['audit'].with_suffix('.npz')) as z:
        triangles=z['triangles'].copy();external_sources=z['sources_0'].copy()
    if sha256(triangles.tobytes()).hexdigest()!=phase['panel_sha256']:
        raise ValueError('physical panel ownership changed')
    panels=TrianglePanels(triangles,-1);cg=inputs.initial_state.position
    field=parents['field'];snapshot=parents['snapshot'];grid=snapshot['grid'];length=snapshot['length_m']
    omega=2*np.pi*field['frequency_hz'];origin=cg[:2]-length/2
    beam=float(np.ptp(inputs.body.mesh.vertices[:,1]));inner=np.max(abs(triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam
    amplitudes=np.array([.01,.005,.01,.005])

    def rotate(values,time):
        cosine,sine=np.cos(omega*time),np.sin(omega*time)
        return np.column_stack((values[:,:2]*cosine+values[:,2:]*sine,
            values[:,2:]*cosine-values[:,:2]*sine))*amplitudes

    def rate(values): return np.column_stack((omega*values[:,2:],-omega*values[:,:2]))

    started=perf_counter()
    mass_factor=splu(patch.mass.tocsc());inverse_cross=mass_factor.solve(cross)
    schur=neumann-nmap@inverse_cross;factor=lu_factor(schur)
    psi_bar=mass_factor.solve(psi_mom)
    far_bar=psi_bar-initial_c-inverse_cross@initial_sigma
    eta_bar=phases(1j*omega/GRAVITY*amplitude(psi_bar))
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    q=np.column_stack((modes[:,[2,3]],np.zeros((len(panels),2))))
    load=WATER_DENSITY*(panels.areas[:,None]*modes).T
    pressure_local=load@pmap;pressure_source=load@even-pressure_local@inverse_cross
    pressure_far=load@far_body
    with np.load(paths['field'].with_suffix('.npz')) as z:
        trace=z['surface_0'].reshape(grid,grid,2);reference_force=z['reference_force'].copy()
    external_density=phases(-external_sources/(-1j*omega))
    axis=np.arange(grid)*length/grid;gx,gy=np.meshgrid(origin[0]+axis,origin[1]+axis)
    gp=np.column_stack((gx.ravel(),gy.ravel(),np.zeros(grid*grid)))
    correction=np.empty((len(gp),4))
    for first in range(0,len(gp),128):
        correction[first:first+128]=2*panel_influence(gp[first:first+128],panels).single_layer@external_density
    trace-=amplitude(correction).reshape(grid,grid,2)
    far=far_trace_spectrum(trace,origin,length,cg[:2],inner,beam)
    kx,ky=np.meshgrid(2*np.pi*np.fft.fftfreq(grid,d=length/grid),2*np.pi*np.fft.fftfreq(grid,d=length/grid))
    far_vertical=far*np.hypot(kx,ky)[:,:,None]
    hx,hy=float(x[1]-x[0]),float(y[1]-y[0])
    far_flux=np.zeros((len(patch.nodes),4))
    far_flux.reshape(len(y),len(x),4)[1:-1,1:-1]=tensor_values(
        far_vertical*hat_fourier(hx,hy,kx,ky)[:,:,None],x[1:-1],y[1:-1],origin,length)
    for points,nodes,weights in cut_quadrature(patch,order,hidden=True):
        scatter_moments(far_flux,nodes,-weights,fourier_values(far_vertical,points[:,:2],origin,length))
    far_flux=far_flux[patch.active_nodes]
    kernel=WeakPatch((len(y)-2,len(x)-2),hx,hy)
    report['kernel_diagnostics']=kernel.diagnostics
    report['constructing_physical_flux']=True
    print(json.dumps(dict(constructing_physical_flux=True,source_triangles=len(patch.faces))),flush=True)
    flux=PhysicalWeakFlux(patch,polygon,kernel,order=args.cut_order)
    report.pop('constructing_physical_flux',None)
    report.update(physical_flux_diagnostics=flux.diagnostics,setup_seconds=perf_counter()-started,
        input_sha256=phase['input_sha256'],panel_sha256=phase['panel_sha256'],phase_order=order,
        active_water_coefficients=count,parallel_phase_cases=4,
        numerical_contract=dict(integrator='SciPy RK23',rtol=RELATIVE_TOLERANCE,
            eta_atol_m=SURFACE_HEIGHT_ABSOLUTE_TOLERANCE_M,
            projected_physical_psi_atol_m2_per_s=SURFACE_POTENTIAL_ABSOLUTE_TOLERANCE_M2_PER_S,
            real_state_count=8*count,error_norm='RMS over ALL real eta and projected physical-psi coefficients in all four cases',
            body_state_count=0,prescribed_velocity_amplitudes=amplitudes.tolist(),first_step=None))
    audit_coefficients=np.column_stack((initial_c,np.random.default_rng(0).standard_normal((count,3)),np.ones(count)))
    saved=dict(x=x,y=y,active_nodes=patch.active_nodes,
        physical_contact_rows=flux.contact_rows,physical_contact_active_rows=flux.affected_rows,
        weak_full_kernel=kernel.kernel,weak_full_shape=np.array(kernel.shape),weak_fft_shape=np.array(kernel.fft_shape),
        weak_spacing=np.array(kernel.spacing),
        operator_audit_coefficients=audit_coefficients,operator_audit_loads=flux.apply(audit_coefficients),
        far_flux=far_flux,initial_eta=eta_bar,initial_projected_psi=psi_bar)
    mass_area=float(patch.mass.sum())

    def rms(values): return np.sqrt(np.sum(values*(patch.mass@values),axis=0)/mass_area)

    def current_coefficients(time,psi):
        direct=psi-rotate(far_bar,time)
        sigma=lu_solve(factor,rotate(q,time)-rotate(far_normal,time)-nmap@direct)
        return direct-inverse_cross@sigma

    def pressure(time,eta):
        direct=-GRAVITY*eta-rotate(rate(far_bar),time)
        sigma_dot=lu_solve(factor,rotate(rate(q),time)-rotate(rate(far_normal),time)-nmap@direct)
        return -pressure_local@direct-pressure_source@sigma_dot-rotate(rate(pressure_far),time)

    initial=np.stack((rotate(eta_bar,0.),rotate(psi_bar,0.)))
    np.testing.assert_allclose(pressure(0.,initial[0]),phases(old_force)*amplitudes,rtol=1e-9,atol=1e-10)
    initial_eta_rate=mass_factor.solve(flux.apply(current_coefficients(0.,initial[1]))+rotate(far_flux,0.))
    report['initial_kinematic_rhs_vs_harmonic_derivative_rms_m_per_s']=rms(initial_eta_rate-rotate(rate(eta_bar),0.)).tolist()
    samples=np.linspace(0.,args.duration,11)
    reports=[];histories=[]
    for max_step in (args.max_step,args.max_step/2):
        evaluations=0;last_evaluation_time=0.
        def rhs(time,vector):
            nonlocal evaluations,last_evaluation_time
            evaluations+=1;last_evaluation_time=float(time)
            if evaluations>args.maximum_rhs: raise RuntimeError('bounded evolving-water RHS budget exhausted')
            eta,psi=vector.reshape(2,count,4)
            coefficient=current_coefficients(time,psi)
            eta_rate=mass_factor.solve(flux.apply(coefficient)+rotate(far_flux,time))
            return np.stack((eta_rate,-GRAVITY*eta)).ravel()
        started=perf_counter()
        atol=np.r_[np.full(count*4,SURFACE_HEIGHT_ABSOLUTE_TOLERANCE_M),
                   np.full(count*4,SURFACE_POTENTIAL_ABSOLUTE_TOLERANCE_M2_PER_S)]
        solver=RK23(rhs,0.,initial.ravel(),args.duration,rtol=RELATIVE_TOLERANCE,atol=atol,max_step=max_step)
        values=[initial.copy()];next_sample=1;steps=0;step_sizes=[]
        while solver.status=='running':
            message=solver.step()
            if solver.status=='failed': raise RuntimeError('RK23 failed: '+str(message))
            steps+=1;step_sizes.append(float(solver.step_size));dense=solver.dense_output()
            while next_sample<len(samples) and samples[next_sample]<=solver.t:
                value=dense(samples[next_sample]).reshape(2,count,4)
                if not np.isfinite(value).all(): raise RuntimeError('nonfinite evolving-water state')
                values.append(value);next_sample+=1
            if steps%10==0:
                print(json.dumps(dict(max_step=max_step,accepted_steps=steps,rhs=evaluations,time=float(solver.t))),flush=True)
        integration_seconds=perf_counter()-started
        if len(values)!=len(samples): raise RuntimeError('missing common-time water output')
        states=np.array(values);started=perf_counter()
        forces=np.array([pressure(t,value[0]) for t,value in zip(samples,states,strict=True)])
        force_seconds=perf_counter()-started
        reference_forces=np.array([rotate(phases(reference_force),t) for t in samples])
        force_error=forces-reference_forces
        complex_force_error=force_error[:,:,:2]+1j*force_error[:,:,2:]
        eta_error=np.array([rms(value[0]-rotate(eta_bar,t)) for t,value in zip(samples,states,strict=True)])
        row=dict(max_step_seconds=max_step,achieved_duration_seconds=float(solver.t),accepted_steps=steps,
            rhs_evaluations=evaluations,last_evaluation_time=last_evaluation_time,
            integration_seconds=integration_seconds,pressure_output_seconds=force_seconds,
            integration_only_wall_per_sim_second=integration_seconds/args.duration,
            minimum_accepted_step=min(step_sizes),maximum_accepted_step=max(step_sizes),
            eta_max_coefficient_m=float(abs(states[:,0]).max()),
            physical_mass_eta_rms_error_vs_external_harmonic_max=eta_error.max(axis=0).tolist(),
            six_force_phase_amplitude_error_max=abs(complex_force_error).max(axis=0).tolist(),
            principal_force_phase_relative_error_max=(abs(complex_force_error).max(axis=0)[[2,3],[0,1]]
                /(abs(reference_force)[[2,3],[0,1]]*amplitudes[:2])).tolist())
        reports.append(row);histories.append((states,forces))
        saved[f'run_{len(reports)}_states']=states;saved[f'run_{len(reports)}_forces']=forces
        saved['times']=samples
        report['runs']=reports
        np.savez(args.output.with_suffix('.npz'),**saved)
        report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
        print(json.dumps(row),flush=True)
    difference=histories[0][0]-histories[1][0]
    report['step_halving_eta_physical_rms_max']=np.max([rms(v[0]) for v in difference],axis=0).tolist()
    report['step_halving_force_absolute_max']=abs(histories[0][1]-histories[1][1]).max(axis=0).tolist()
    report.update(complete=True,qualified=False,scope=__doc__,peak_process_rss_bytes=peak_rss_bytes(),
        limits=['Fixed geometry and prescribed body velocity/far chi; body pressure does not advance six body states',
            'Both eta and projected physical psi evolve; no external force/eta reset is applied during integration',
            'Short-time comparison with a finite-panel harmonic reference includes spatial/Robin/far-closure error, not pure temporal error',
            'Local RK tolerance acceptance is not body-motion or water-field accuracy qualification',
            'Integration-only timing is NOT full nonlinear6DOF simulation throughput; source preparation/output are separately reported',
            'No moving geometry, global energy certificate,60second0.1degree/1cm accuracy or baseline promotion'])
    report['measurement_environment']=measurement_environment(cp,'unknown',
        'User released competing workload; no independent CPU contention measurement. Numerical jobs serialized.',
        'Common-time states and pressures; artifacts written after each integration, no speed ratio between different max-step contracts')
    report['measurement_environment']['timing_scope']='Fixed-hull verification setup, RK integration and pressure-output timings reported separately; total includes result writes, not imports/environment launcher'
    report['measurement_environment']['startup_treatment']='Fresh process; first-use kernel/GPU context setup included in setup; uncontrolled persistent caches; steady RHS uses batched GPU FFT/host transfers and CPU sparse/body maps'


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('phase','output'): parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--duration',type=float,default=.1)
    parser.add_argument('--max-step',type=float,default=.01)
    parser.add_argument('--cut-order',type=int,choices=[4,8],default=4)
    parser.add_argument('--maximum-rhs',type=int,default=2000)
    args=parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists() or not 0<args.duration<=1
            or not 0<args.max_step<=args.duration or args.maximum_rhs<1):
        parser.error('fresh archive, duration(0,1], positive max-step<=duration and positive RHS budget required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parent_phase_sha256=sha256(args.phase.read_bytes()).hexdigest())
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
