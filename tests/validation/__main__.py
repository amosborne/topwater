"""Executable canonical validation sequence. Reference generation is explicit.

python -m tests.validation freeze --output artifacts/validation/frozen.json
python -m tests.validation reference --frozen ... --output ...
python -m tests.validation compare-inertia --frozen ... --reference ... --output ...
python -m tests.validation compare-waves --output ...  # currently FAILS unsupported
python -m tests.validation qualify --frozen ... --reference ... --inertia ... --waves ... --output ...
"""
import argparse
from pathlib import Path
from time import perf_counter
import numpy as np
from threadpoolctl import threadpool_limits
from .common import (PHYSICAL,CANDIDATE,REFERENCE,VALIDATORS,case_envelope,physical_identity,digest,hashes,load_contract,prepare,
    save,read_owned,current_sources,study_preflight,triangles_at,entrywise)


def freeze(args,contract,report,arrays):
    from tests import reference_boat
    from topwater import SimulationInputs,WaveComponents
    from topwater.frames import rotation_matrix
    p=contract['physical_problem']
    if (p['density_kg_m3']!=reference_boat.WATER_DENSITY
            or p['gravity_m_s2']!=reference_boat.GRAVITY):
        raise ValueError('contract constants differ from actual canonical equilibrium fixture')
    body=reference_boat.reference_body();pose,hydro=reference_boat.equilibrium()
    if np.any(pose.linear_velocity) or np.any(pose.angular_velocity):
        raise ValueError('canonical reference must freeze a calm zero-velocity equilibrium')
    inputs=SimulationInputs(body,WaveComponents([],[],[],[]),pose)
    text=inputs.to_json();input_path=args.output.with_suffix('.source')/'inputs.json'
    input_path.write_text(text)
    R=rotation_matrix(pose.orientation)
    dry=np.zeros((6,6));dry[:3,:3]=body.mass*np.eye(3);dry[3:,3:]=R@body.inertia_cg@R.T
    length=float(np.ptp(body.mesh.vertices[np.unique(body.mesh.faces[body.mesh.physical_faces]),0]))
    beam=float(np.ptp(body.mesh.vertices[np.unique(body.mesh.faces[body.mesh.physical_faces]),1]))
    report['incident_case_envelope']=case_envelope(contract,beam)
    arrays.update(wetted_triangles=hydro.wetted_triangles,world_cg=pose.position,
        orientation=pose.orientation,dry_mass=dry,hydrostatic_stiffness=hydro.stiffness,
        body_mass=np.array(body.mass),hull_length=np.array(length),hull_beam=np.array(beam))
    report.update(input_sha256=digest(input_path),physical_source_sha256=hashes(PHYSICAL),
        implementation=contract['implementation'],physical_problem=p,
        physical_freeze_only=True,qualification_scope='Not physical solver qualification',
        equilibrium_force_N=hydro.force.tolist(),equilibrium_moment_Nm=hydro.moment.tolist(),
        wet_panel_count=len(hydro.wetted_triangles),complete=False)
    report['physical_identity_sha256']=physical_identity(report,arrays)
    save(args.output,report,arrays)
    report['reference_study_plan']=study_preflight(len(hydro.wetted_triangles),
        contract['reference_numerics']['body_subdivision_levels'],
        contract['reference_numerics']['maximum_estimated_workspace_bytes'])
    report.update(complete=True,passed=True)


def frozen_input(args,contract,report):
    frozen,data=read_owned(args.frozen,report['contract_sha256'],'freeze',passed=True)
    current_sources(frozen,PHYSICAL+VALIDATORS)
    case_envelope(contract,float(data['hull_beam']))
    if digest(Path(args.frozen).with_suffix('.source')/'inputs.json')!=frozen['input_sha256']:
        raise ValueError('frozen physical input JSON changed')
    if physical_identity(frozen,data)!=frozen['physical_identity_sha256']:
        raise ValueError('frozen stable physical identity differs from its actual arrays/inputs')
    report.update(input_sha256=frozen['input_sha256'],
        frozen_physical_identity_sha256=frozen['physical_identity_sha256'],
        frozen_manifest_sha256_provenance=digest(args.frozen))
    return frozen,data


def qualified_reference(args,contract,report,frozen,frozen_data):
    from .reference import require_qualified
    reference,data=read_owned(args.reference,report['contract_sha256'],'reference',passed=True)
    if (reference['input_sha256']!=frozen['input_sha256']
            or reference['frozen_physical_identity_sha256']!=frozen['physical_identity_sha256']):
        raise ValueError('qualified reference belongs to different frozen inputs')
    current_sources(reference,REFERENCE)
    require_qualified(reference,data,contract,frozen_data)
    report['reference_sha256']=digest(args.reference)
    return reference,data


def compare_inertia(args,contract,report,arrays):
    from topwater.inertia import instantaneous_added_mass
    from topwater.rankine import TrianglePanels
    from topwater.numerics import LaplaceControls
    frozen,data=frozen_input(args,contract,report)
    reference,ref=qualified_reference(args,contract,report,frozen,data)
    spec=contract['candidate_numerics'];resource=contract['reference_numerics']
    levels=spec['body_subdivision_levels']
    report['study_plan']=study_preflight(len(data['wetted_triangles']),levels,resource['maximum_estimated_workspace_bytes'])
    report['candidate_source_sha256']=hashes(CANDIDATE)
    controls=LaplaceControls(quadrature_rtol=spec['quadrature_rtol'],quadrature_atol=spec['quadrature_atol'],
        max_refinement=spec['max_refinement'],target_block_size=spec['target_block_size'],
        max_operator_bytes=resource['maximum_estimated_workspace_bytes'])
    rho=contract['physical_problem']['density_kg_m3'];started=perf_counter()
    scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:])
    values=[];report['levels']=[]
    for level in levels:
        if perf_counter()-started>resource['maximum_elapsed_seconds']:
            raise TimeoutError('candidate elapsed resource ceiling exhausted between physics calls')
        triangles=triangles_at(data['wetted_triangles'],level)
        with threadpool_limits(1):
            result=instantaneous_added_mass(TrianglePanels(triangles,-1),data['world_cg'],density=rho,controls=controls)
        if not np.isfinite(result.added_mass).all():
            raise FloatingPointError('nonfinite candidate all-six response')
        arrays[f'level{level}_Ainf']=result.added_mass
        values.append(result.added_mass/scale)
        report['levels'].append(dict(level=level,panels=len(triangles),diagnostics=result.diagnostics,
            comparison=entrywise(values[-1],ref['normalized_Ainf'][-1],ref['uncertainty_Ainf'],
                contract['acceptance']['candidate_absolute_normalized'],contract['acceptance']['candidate_relative'])))
        save(args.output,report,arrays)
    if perf_counter()-started>resource['maximum_elapsed_seconds']:
        raise TimeoutError('candidate elapsed ceiling exceeded after final solve; calls are not interruptible')
    acceptance=contract['acceptance']
    # The declared comparison tolerance is relative to the finest C, not Cprevious.
    allowance=acceptance['candidate_absolute_normalized']+acceptance['candidate_relative']*abs(values[-1])
    mesh=dict(passed=bool(np.all(abs(values[-1]-values[-2])<=allowance)),
        maximum_tolerance_ratio=float(np.max(abs(values[-1]-values[-2])/allowance)))
    report.update(complete=True,passed=mesh['passed'] and report['levels'][-1]['comparison']['passed'],
        gate='all_six_added_mass',candidate_mesh_gate=mesh,
        qualification_scope=acceptance['scope'])
    arrays['normalized_Ainf']=np.stack(values)


def compare_waves(args,contract,report,arrays):
    # No alternate solver, supplied force table or optional adapter can turn
    # the selected instantaneous-only core into a time-domain wave path.
    report.update(complete=True,passed=False,gate='small_amplitude_time_domain_wave_responses',
        status='unsupported',reason='Selected P0 Rankine core has no coupled time-domain wave implementation',
        required_contract_revision=contract['time_domain_wave_controls'])


def qualify(args,contract,report,arrays):
    frozen,frozen_data=frozen_input(args,contract,report)
    rows={}
    for name,kind,gate in (('reference','reference','reference'),('inertia','compare-inertia','all_six_added_mass'),
                           ('waves','compare-waves','small_amplitude_time_domain_wave_responses')):
        try:
            item,_=read_owned(getattr(args,name),report['contract_sha256'],kind,passed=True)
            if item.get('input_sha256')!=frozen['input_sha256']:
                raise ValueError('gate frozen physical input identity differs')
            if item.get('frozen_physical_identity_sha256')!=frozen['physical_identity_sha256']:
                raise ValueError('gate stable physical identity differs')
            if kind=='reference':
                qualified_reference(args,contract,report,frozen,frozen_data)
            elif kind=='compare-inertia' and item.get('candidate_source_sha256')!=hashes(CANDIDATE):
                raise ValueError('candidate implementation source changed after comparison')
            if kind=='compare-inertia':
                current_sources(item,CANDIDATE+VALIDATORS)
            if kind=='compare-inertia' and item.get('reference_sha256')!=digest(args.reference):
                raise ValueError('candidate comparison uses a different reference artifact')
            rows[gate]=dict(passed=True,artifact_sha256=digest(getattr(args,name)))
        except (OSError,ValueError,KeyError) as error:
            rows[gate]=dict(passed=False,reason=str(error))
    rows['small_amplitude_time_domain_wave_responses']=dict(passed=False,status='unsupported',
        reason='Selected implementation has no TD-wave path; a supplied passed artifact cannot override this capability')
    report.update(complete=True,passed=False,gates=rows,qualification_status='unqualified',
        reason='All required physical gates must pass; missing/unsupported cases never pass')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=('freeze','reference','compare-inertia','compare-waves','qualify'))
    parser.add_argument('--contract',type=Path,default=Path('validation/canonical-boat.json'))
    for name in ('frozen','reference','inertia','waves'):
        parser.add_argument('--'+name,type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    required={'freeze':(),'reference':('frozen',),'compare-inertia':('frozen','reference'),
        'compare-waves':(),'qualify':('frozen','reference','inertia','waves')}[args.command]
    if any(getattr(args,name) is None for name in required):
        parser.error('missing command inputs: '+','.join(required))
    contract=load_contract(args.contract)
    report=prepare(args.output,args.contract,args.command);arrays={};started=perf_counter()
    try:
        if args.command=='reference':
            from .reference import generate
            frozen,data=frozen_input(args,contract,report)
            study_preflight(len(data['wetted_triangles']),contract['reference_numerics']['body_subdivision_levels'],
                contract['reference_numerics']['maximum_estimated_workspace_bytes'])
            generate(args,contract,frozen,data,report,arrays)
        else:
            {'freeze':freeze,'compare-inertia':compare_inertia,'compare-waves':compare_waves,'qualify':qualify}[args.command](
                args,contract,report,arrays)
    except (Exception,KeyboardInterrupt) as error:
        report.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    report['elapsed_seconds']=perf_counter()-started
    save(args.output,report,arrays)
    print(__import__('json').dumps(dict(output=str(args.output),complete=report['complete'],passed=report['passed'])))
    return 0 if report['complete'] and report['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
