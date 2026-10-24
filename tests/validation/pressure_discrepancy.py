"""Separate coarse body and explicit free-surface/domain effects; no promotion."""
import argparse
from importlib import import_module
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from topwater import SimulationInputs
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from topwater.inertia import instantaneous_added_mass
from .common import prepare,save,read_owned,digest,current_sources,entrywise
from .reference_mesh import farthest_first,verify_surface


def generate(args,c,r,a):
    if (c['format_version']!=1 or c['maximum_elapsed_seconds']!=60 or c['cpu_threads']!=1
            or c['required_capytaine_calls']!=6 or c['qualification_allowed'] is not False
            or c['diagnostic_relative']!=.01 or c['diagnostic_absolute_normalized']!=1e-5
            or c['density_kg_m3']!=1025. or c['gravity_m_s2']!=9.81):
        raise ValueError('unsupported discrepancy attribution controls')
    parent,p=read_owned(args.probe,digest('validation/mixed-pressure-controls.json'),'mixed-pressure-component')
    current_sources(parent,tuple(parent['producer_sha256']))
    if (parent['full_reference_qualified'] is not False or parent['nonlinear_rhs_ready'] is not False
            or parent['complete_runtime_speed_qualified'] is not False
            or parent['scope']!='Bounded calm pressure component and nonflat boundary closure; no TD speed pass.'):
        raise ValueError('parent scope or unqualified status differs')
    archive=args.output.with_suffix('.source');source=args.probe.with_suffix('.source')
    if digest(source/'calm-inputs.json')!=parent['calm_inputs_sha256']:raise ValueError('probe physical input differs')
    physical=(source/'calm-inputs.json').read_bytes();(archive/'calm-inputs.json').write_bytes(physical)
    if digest(source/'protocol.md')!=parent['protocol_sha256']:raise ValueError('parent protocol differs')
    shutil.copyfile(source/'contract.json',archive/'probe-controls.json')
    shutil.copyfile(source/'protocol.md',archive/'probe-protocol.md')
    shutil.copyfile(source/'field-map.zip',archive/'field-map.zip');shutil.copyfile(c['protocol'],archive/'protocol.md')
    if digest(archive/'field-map.zip')!=parent['field_bundle_sha256']:raise ValueError('parent evidence bundle differs')
    for suffix in ('.json','.npz'):shutil.copyfile(args.probe.with_suffix(suffix),archive/('probe'+suffix))
    mean=SimulationInputs.from_json(physical.decode());R=rotation_matrix(mean.initial_state.orientation)
    mass=np.zeros((6,6));mass[:3,:3]=mean.body.mass*np.eye(3);mass[3:,3:]=R@mean.body.inertia_cg@R.T
    scale=np.sqrt(mass.diagonal()[:,None]*mass.diagonal()[None,:])
    triangles=farthest_first(p['calm_vertices'][p['calm_faces'][p['calm_tags']=='body']][:,::-1])
    a.update(body_triangles=triangles,world_cg=mean.initial_state.position,dry_mass=mass,
        closed_added_mass=p['calm_added_mass'],fine_reference_Ainf=p['raw_reference_Ainf'],
        fine_reference_uncertainty=p['reference_Ainf_increment_envelope'],entry_scale=scale)
    r.update(probe_manifest_sha256=digest(args.probe),probe_arrays_sha256=digest(args.probe.with_suffix('.npz')),
        calm_inputs_sha256=digest(archive/'calm-inputs.json'),field_bundle_sha256=parent['field_bundle_sha256'],
        protocol_sha256=digest(c['protocol']),completed_capytaine_calls=0,full_reference_qualified=False,
        baseline_qualified=False,scope=__doc__)
    save(args.output,r,a);started=perf_counter()
    def budget():
        if perf_counter()-started>c['maximum_elapsed_seconds']:raise TimeoutError('discrepancy diagnostic exceeded time budget')
    original=hydrostatics(mean.body,mean.initial_state,density=c['density_kg_m3'],gravity=c['gravity_m_s2'])
    r['physical_surface_proof']=verify_surface(original.wetted_triangles,triangles)
    save(args.output,r,a)
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles,MODES
        r['native_thread_pools']=threadpool_info()
        if cpt.__version__!='3.0.0' or any(x['num_threads']!=1 for x in threadpool_info()):raise ValueError('pinned Capytaine and one thread required')
        budget();begin=perf_counter()
        top=instantaneous_added_mass(TrianglePanels(triangles,-1),mean.initial_state.position,density=c['density_kg_m3'])
        a['halfspace_added_mass']=top.added_mass;r['halfspace_seconds']=perf_counter()-begin
        r['halfspace_diagnostics']=top.diagnostics;budget();save(args.output,r,a)
        body=body_from_triangles(triangles,mean.initial_state.position,mean.body.mass)
        green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None)
        r['green_settings']=green.exportable_settings;r['capytaine_version']=cpt.__version__
        solver=cpt.BEMSolver(green_function=green);coarse=np.zeros((6,6));a['same_body_capytaine_Ainf']=coarse
        for j,mode in enumerate(MODES):
            budget();result=solver.solve(cpt.RadiationProblem(body=body,radiating_dof=mode,omega=np.inf,
                rho=c['density_kg_m3'],g=c['gravity_m_s2'],water_depth=np.inf,forward_speed=0.))
            coarse[:,j]=[result.added_mass[m] for m in MODES]
            if any(result.radiation_damping[m]!=0 for m in MODES):raise ValueError('nonzero infinite-frequency damping')
            r['completed_capytaine_calls']+=1;budget();save(args.output,r,a)
    comparisons={}
    for name,x,y,u in (
            ('closed_vs_halfspace',a['closed_added_mass'],a['halfspace_added_mass'],np.zeros((6,6))),
            ('halfspace_vs_same_body_capytaine',a['halfspace_added_mass'],coarse,np.zeros((6,6))),
            ('same_body_capytaine_vs_fine',coarse,a['fine_reference_Ainf'][-1],a['fine_reference_uncertainty']),
            ('halfspace_vs_fine',a['halfspace_added_mass'],a['fine_reference_Ainf'][-1],a['fine_reference_uncertainty'])):
        a[name+'_normalized_difference']=(x-y)/scale
        comparisons[name]=entrywise(x/scale,y/scale,u/scale,c['diagnostic_absolute_normalized'],c['diagnostic_relative'])
    budget();r.update(complete=True,passed=False,diagnostic_comparisons=comparisons,elapsed_seconds=perf_counter()-started,
        qualification_status='unqualified; same-body and finite-domain comparisons diagnose error only, not convergence')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--controls',type=Path,default=Path('validation/pressure-discrepancy-controls.json'))
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();c=json.loads(args.controls.read_text());r=prepare(args.output,args.controls,'pressure-discrepancy-diagnostic');a={};save(args.output,r,a)
    try:generate(args,c,r,a)
    except (Exception,KeyboardInterrupt) as error:r.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,r,a);print(json.dumps({k:r.get(k) for k in ('complete','passed','completed_capytaine_calls','diagnostic_comparisons','error')}))
    return 0 if r['complete'] and r['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
