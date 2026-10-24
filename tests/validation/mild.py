"""Array-only periodic motion candidate; no Topwater runtime or transient solver."""
import argparse
import json
from pathlib import Path
import shutil
from zipfile import ZipFile
import numpy as np
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import InitialState,SimulationInputs,WaveComponents
from topwater.frames import rotation_matrix
from .common import (prepare,save,read_owned,digest,current_sources,REFERENCE,
    PHYSICAL,VALIDATORS,physical_identity)


def histories(phasors,omega,times):
    """World generalized displacement and its exact first/second derivatives."""
    p=np.asarray(phasors);w=np.asarray(omega);t=np.asarray(times)
    if (p.ndim!=3 or p.shape[1:]!=(len(w),6) or t.ndim!=1
            or not all(np.isfinite(x).all() for x in (p,w,t)) or np.any(w<=0)):
        raise ValueError('finite level/wave/six-mode phasors and positive frequencies required')
    carrier=np.exp(-1j*w[:,None]*t[None,:])
    return tuple(np.einsum('lwj,wt->ltj',p*(-1j*w[None,:,None])**order,carrier).real
                 for order in range(3))


def spatial_envelope(phasors,omega,step,limits,fraction,floor):
    p=np.asarray(phasors)
    if len(p)<4 or not np.isfinite(p).all():
        raise ValueError('four complete complex motion levels required')
    previous=abs(p[-2]-p[-3]);latest=abs(p[-1]-p[-2])
    component=np.maximum(previous,latest)
    spatial=component.sum(axis=0)
    # Every time is within step/2 of a sample; the exact harmonic derivative
    # gives a conservative reference-only interpolation/sampling bound.
    sampling=(abs(p[-1])*omega[:,None]).sum(axis=0)*step/2
    allowance=np.asarray(limits)*fraction
    trend=bool(np.all(latest<=previous+floor[None,:]))
    return dict(spatial=spatial,sampling=sampling,component=component,
        velocity=(component*omega[:,None]).sum(axis=0),
        acceleration=(component*omega[:,None]**2).sum(axis=0),
        gate=dict(passed=bool(np.all(spatial+sampling<=allowance) and trend),
            trend_passed=trend,maximum_budget_ratio=float(np.max((spatial+sampling)/allowance))))


def initial_body(mean,displacement,velocity):
    """Physical pose and body angular velocity from world rotation coordinates."""
    theta=np.asarray(displacement[3:]);rate=np.asarray(velocity[3:]);angle=np.linalg.norm(theta)
    x,y,z=theta;cross=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
    if angle<1e-5:
        a=.5-angle**2/24+angle**4/720
        b=1/6-angle**2/120+angle**4/5040
    else:
        a=(1-np.cos(angle))/angle**2;b=(angle-np.sin(angle))/angle**3
    world_omega=(np.eye(3)+a*cross+b*(cross@cross))@rate
    world_R=Rotation.from_rotvec(theta).as_matrix()@rotation_matrix(mean.orientation)
    return InitialState(position=mean.position+displacement[:3],
        orientation=Rotation.from_matrix(world_R).as_quat(scalar_first=True),
        linear_velocity=velocity[:3],angular_velocity=world_R.T@world_omega)


def informative_signal(motion,uncertainty,limits):
    signal=np.max(abs(motion-motion[0]),axis=0)
    lower=np.maximum(0.,signal-2*uncertainty)
    return signal,lower,bool(np.any(lower[:3]>limits[:3]) and np.any(lower[3:]>limits[3:]))


def load_case(path):
    c=json.loads(Path(path).read_text());s=c['sea'];q=c['comparison'];init=c['initialization']
    arrays=[np.asarray(s[key],float) for key in ('amplitudes_m','frequencies_hz','propagation_headings_degrees','phases_rad')]
    if (c['format_version']!=1 or any(a.shape!=(3,) or not np.isfinite(a).all() for a in arrays)
            or np.any(arrays[0]<=0) or np.any(arrays[1]<=0)
            or s['time_convention']!='Re(amplitude*exp(-i*omega*t))'
            or s['phase_origin_world_xy_m']!=[0.,0.]
            or init['rest_start_reference'] is not False
            or init['kind']!='periodic steady state with prehistory extending before t=0'
            or q['duration_s']<1 or not 0<q['sample_step_s']<=.005
            or q['translation_error_limit_m']!=.01 or q['rotation_error_limit_degrees']!=.1
            or q['rotation_error_chart']!='theta=Log(R*R_equilibrium^T), components in mean world axes; compare the same chart in both models'
            or q['reference_budget_fraction']!=.1
            or q['require_nontrivial_translation_and_rotation'] is not True
            or q['reference_roundoff_floor_m']!=1e-10 or q['reference_roundoff_floor_rad']!=1e-10
            or q['runtime_limit_wall_seconds_per_simulated_second']!=200):
        raise ValueError('unsupported prospective mild periodic case/acceptance')
    for name in ('duration_s','sample_step_s'):
        if not np.isfinite(q[name]):raise ValueError('finite reference duration and step required')
    for name in ('selected_frequency_indices','selected_heading_indices'):
        indices=c['raw_reference'][name]
        if indices!=[0,1,2] or any(type(i)!=int for i in indices):
            raise ValueError('exact declared discrete raw response coverage required')
    if c['raw_reference']['required_levels']!=[1,2,3,4]:
        raise ValueError('all four raw mesh levels required')
    return c


def raw_inputs(args,report):
    contract=json.loads(Path('validation/canonical-boat.json').read_text())
    raw,arrays=read_owned(args.raw,digest('validation/canonical-boat.json'),'reference')
    frozen,data=read_owned(args.frozen,raw['contract_sha256'],'freeze',passed=True)
    current_sources(raw,REFERENCE);current_sources(frozen,PHYSICAL+VALIDATORS)
    if (physical_identity(frozen,data)!=frozen['physical_identity_sha256']
            or raw['frozen_physical_identity_sha256']!=frozen['physical_identity_sha256']
            or raw['completed_physics_calls']!=132 or len(raw['levels'])!=4
            or not all(row['complete'] for row in raw['levels'])):
        raise ValueError('incomplete or identity-mismatched raw physical reference')
    if (digest(args.frozen.with_suffix('.source')/'inputs.json')!=frozen['input_sha256']
            or frozen['input_sha256']!=raw['input_sha256']):
        raise ValueError('frozen input JSON identity differs')
    shapes=dict(Ainf=(4,6,6),A=(4,3,6,6),B=(4,3,6,6),
        **{key:(4,3,3,6) for key in ('diffraction','FK','excitation','RAO')})
    for name,shape in shapes.items():
        if arrays[name].shape!=shape or not np.isfinite(arrays[name]).all():
            raise ValueError('incomplete raw reference quantity '+name)
    np.testing.assert_array_equal(arrays['excitation'],arrays['diffraction']+arrays['FK'])
    for level in range(4):
        np.testing.assert_array_equal(arrays[f'triangles_{level+1}'],data[f'reference_triangles_{level+1}'])
        for f,w in enumerate(2*np.pi*np.asarray(contract['cases']['frequencies_hz'])):
            Z=data['hydrostatic_stiffness']-w*w*(data['dry_mass']+arrays['A'][level,f])-1j*w*arrays['B'][level,f]
            np.testing.assert_allclose(arrays['RAO'][level,f],np.linalg.solve(Z,arrays['excitation'][level,f].T).T,rtol=1e-13,atol=1e-13)
    report.update(raw_manifest_sha256=digest(args.raw),frozen_manifest_sha256=digest(args.frozen),
        frozen_physical_identity_sha256=frozen['physical_identity_sha256'],
        original_raw_passed=raw['passed'],original_raw_gates=raw['gates'])
    return contract,raw,arrays,frozen,data


@threadpool_limits.wrap(limits=1)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,default=Path('validation/mild-periodic-case.json'))
    for name in ('raw','raw-bundle','frozen','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();case=load_case(args.case)
    report=prepare(args.output,args.case,'periodic-motion-candidate');out={}
    report['native_thread_pools']=threadpool_info()
    save(args.output,report,out)
    try:
        if digest(args.raw_bundle)!=case['raw_reference']['bundle_sha256']:
            raise ValueError('raw bundle differs from the prospectively pinned data')
        with ZipFile(args.raw_bundle) as bundle:
            for suffix in ('.json','.npz'):
                if bundle.read('validation/'+case['raw_reference']['study_id']+suffix)!=args.raw.with_suffix(suffix).read_bytes():
                    raise ValueError('raw file bytes differ from the pinned bundle')
        shutil.copyfile(args.raw_bundle,args.output.with_suffix('.source')/'raw-reference.zip')
        contract,raw,r,frozen,data=raw_inputs(args,report)
        s=case['sea'];q=case['comparison'];spec=case['raw_reference']
        f=np.asarray(spec['selected_frequency_indices']);h=np.asarray(spec['selected_heading_indices'])
        np.testing.assert_array_equal(np.asarray(contract['cases']['frequencies_hz'])[f],s['frequencies_hz'])
        np.testing.assert_allclose(np.asarray(contract['cases']['propagation_headings_rad'])[h],np.deg2rad(s['propagation_headings_degrees']),rtol=0,atol=1e-15)
        omega=2*np.pi*np.asarray(s['frequencies_hz']);amplitude=np.asarray(s['amplitudes_m']);phase=np.asarray(s['phases_rad'])
        phasors=r['RAO'][:,f,h,:]*amplitude[None,:,None]*np.exp(1j*phase[None,:,None])
        duration,step=q['duration_s'],q['sample_step_s'];count=round(duration/step)
        if abs(count*step-duration)>1e-12:raise ValueError('duration must be an integer multiple of sample step')
        times=np.linspace(0,duration,count+1);motion,velocity,acceleration=histories(phasors,omega,times)
        refined=histories(phasors,omega,np.linspace(0,duration,2*count+1))[0]
        np.testing.assert_allclose(motion,refined[:,::2],rtol=0,atol=1e-14)
        limits=np.array([q['translation_error_limit_m']]*3+[np.deg2rad(q['rotation_error_limit_degrees'])]*3)
        floor=np.array([q['reference_roundoff_floor_m']]*3+[q['reference_roundoff_floor_rad']]*3)
        uncertainty=spatial_envelope(phasors,omega,step,limits,q['reference_budget_fraction'],floor)
        signal,lower,informative=informative_signal(motion[-1],uncertainty['spatial']+uncertainty['sampling'],limits)
        mean=SimulationInputs.from_json(args.frozen.with_suffix('.source').joinpath('inputs.json').read_text(encoding='utf-8'))
        waves=WaveComponents(amplitude,omega,np.deg2rad(s['propagation_headings_degrees']),phase,
            provenance=dict(case_id=case['case_id'],kind=s['kind'],initialization=case['initialization']['kind']))
        body=initial_body(mean.initial_state,motion[-1,0],velocity[-1,0])
        physical=SimulationInputs(mean.body,waves,body).to_json()
        target=args.output.with_suffix('.source')/'body-wave-inputs.json';target.write_text(physical,encoding='utf-8',newline='\n')
        out.update(times_s=times,angular_frequencies_rad_s=omega,component_phasors=phasors,
            displacement=motion,velocity=velocity,acceleration=acceleration,
            spatial_uncertainty=uncertainty['spatial'],reference_sampling_bound=uncertainty['sampling'],
            velocity_spatial_uncertainty=uncertainty['velocity'],acceleration_spatial_uncertainty=uncertainty['acceleration'])
        for name in ('Ainf','A','B','excitation'):
            out['raw_'+name]=r[name]
            out['raw_mesh_increment_envelope_'+name]=np.maximum(abs(r[name][-2]-r[name][-3]),abs(r[name][-1]-r[name][-2]))
        report.update(complete=True,passed=False,case_id=case['case_id'],body_wave_inputs_sha256=digest(target),
            equivalent_component_Hm0_m=float(4*np.sqrt(np.sum(amplitude**2)/2)),
            motion_spatial_gate=uncertainty['gate'],case_informative=informative,
            constant_initial_pose_error=signal.tolist(),signal_to_user_tolerance=(signal/limits).tolist(),
            constant_initial_pose_error_lower_bound=lower.tolist(),
            motion_coordinate_convention=q['rotation_error_chart'],
            initial_body_state={name:getattr(body,name).tolist() for name in ('position','orientation','linear_velocity','angular_velocity')},
            reference_budget=(limits*q['reference_budget_fraction']).tolist(),
            qualification_status='unqualified: periodic fluid initial fields, field uncertainty, compatibility and nonlinear-limit checks are missing',
            required_readiness=dict(motion_spatial=uncertainty['gate']['passed'],informative_case=informative,
                periodic_fluid_fields=False,field_uncertainty=False,nonlinear_compatibility=False,linear_limit=False),
            runtime_status='unsupported; no Topwater TD simulation executed',
            scope='Periodic body-motion spatial candidate only; original raw coefficient decision unchanged')
    except (Exception,KeyboardInterrupt) as error:
        report.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,report,out)
    print(json.dumps(dict(output=str(args.output),complete=report['complete'],passed=report['passed'],motion_gate=report.get('motion_spatial_gate'))))
    return 0 if report['complete'] and report['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
