"""Offline periodic indirect sources; not a compatible nonlinear fluid restart."""
import argparse
from importlib import import_module
import json
from pathlib import Path
import shutil
from time import perf_counter
from zipfile import ZipFile
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from .common import prepare,save,digest,reference_plan
from .mild import load_case,raw_inputs


def combine(basis,zeta,xi):
    """Radiation columns followed by unit-elevation diffraction; no incident field."""
    basis=np.asarray(basis);xi=np.asarray(xi)
    if (basis.ndim!=2 or basis.shape[1]!=7 or xi.shape!=(6,)
            or not np.isfinite(basis).all() or not np.isfinite(xi).all()
            or not np.isfinite(zeta)):
        raise ValueError('finite six radiation columns and one diffraction column required')
    return basis[:,:6]@xi+basis[:,6]*zeta


def integrate_pressure(triangles,center,pressure):
    """Body-outward geometry; pressure force acts inward on the body."""
    t=np.asarray(triangles);p=np.asarray(pressure)
    area_normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])/2
    load=np.concatenate((area_normal,np.cross(t.mean(axis=1)-center,area_normal)),axis=1)
    return -load.T@p


def require_close(actual,expected,scale,relative=1e-11,absolute=1e-12):
    actual,expected=np.asarray(actual),np.asarray(expected)
    if actual.shape!=expected.shape or not np.isfinite(actual).all() or not np.isfinite(expected).all():
        raise ValueError('missing or nonfinite reproduced reference')
    ratio=abs(actual-expected)/(absolute*np.asarray(scale)+relative*abs(expected))
    if np.any(ratio>1):
        raise ValueError(f'reference reproduction failed: maximum ratio {np.max(ratio)}')
    return float(np.max(ratio))


def generate(args,controls,case,report,out):
    contract,raw,r,frozen,data=raw_inputs(args,report)
    spec=contract['reference_numerics'];physical=contract['physical_problem']
    if (controls['format_version']!=1 or controls['required_physics_calls']!=84 or controls['cpu_threads']!=1
            or controls['maximum_elapsed_seconds']!=1800
            or controls['maximum_estimated_workspace_bytes']!=spec['maximum_estimated_workspace_bytes']
            or controls['reproduction_relative']!=1e-11
            or controls['reproduction_absolute_normalized']!=1e-12
            or controls['trace_relative']!=1e-12 or controls['trace_absolute_normalized']!=1e-12):
        raise ValueError('unsupported prospective source-recovery controls')
    report['study_plan']=reference_plan(data,spec)
    report.update(completed_physics_calls=0,levels=[],checks=[],scope=__doc__,
        full_reference_qualified=False,nonlinear_restart_ready=False,
        original_raw_passed=raw['passed'])
    archive=args.output.with_suffix('.source')
    shutil.copyfile(controls['case'],archive/'mild-case.json')
    shutil.copyfile(controls['protocol'],archive/'protocol.md')
    if digest(args.raw_bundle)!=case['raw_reference']['bundle_sha256']:
        raise ValueError('raw bundle differs from frozen mild case')
    with ZipFile(args.raw_bundle) as bundle:
        for suffix in ('.json','.npz'):
            if bundle.read('validation/'+case['raw_reference']['study_id']+suffix)!=args.raw.with_suffix(suffix).read_bytes():
                raise ValueError('raw artifact differs from pinned bundle')
    shutil.copyfile(args.raw_bundle,archive/'raw-reference.zip')
    report.update(case_sha256=digest(controls['case']),protocol_sha256=digest(controls['protocol']),
                  raw_bundle_sha256=digest(args.raw_bundle))
    # This save precedes native Green construction and all numerical solves.
    save(args.output,report,out)
    started=perf_counter()
    import_module('capytaine.green_functions.Delhommeau_float64')
    with threadpool_limits(1):
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles,MODES
        from capytaine.bem.airy_waves import airy_waves_potential,airy_waves_velocity,froude_krylov_force
        report['native_thread_pools']=threadpool_info()
        if cpt.__version__!=spec['capytaine_version'] or any(p['num_threads']!=1 for p in threadpool_info()):
            raise ValueError('pinned Capytaine and single native thread required')
        report['capytaine_version']=cpt.__version__
        save(args.output,report,out)
        def budget():
            if perf_counter()-started>controls['maximum_elapsed_seconds']:
                raise TimeoutError('source-recovery elapsed budget exhausted')
        budget()
        green=cpt.Delhommeau(tabulation_cache_dir=None)
        budget()
        report['green_settings']=green.exportable_settings
        if green.exportable_settings!=raw['green_settings']['finite']:
            raise ValueError('finite Green settings differ from pinned raw study')
        rho,g=physical['density_kg_m3'],physical['gravity_m_s2']
        report.update(density_kg_m3=rho,gravity_m_s2=g,water_depth='infinite',forward_speed_m_s=0.,
            time_convention='Re(amplitude*exp(-i*omega*t))',
            basis_columns=list(MODES)+['Diffraction'],
            basis_normalization='Unit displacement (m or rad) radiation; unit elevation (m) diffraction',
            source_convention='Capytaine 3.0 indirect density with archived Delhommeau settings; not Rankine density',
            trace_units=dict(potential='m^2/s per unit basis amplitude',pressure='Pa per unit basis amplitude',
                normal_velocity='m/s per unit basis amplitude',
                combined_potential='m^2/s',combined_pressure='Pa',combined_normal_velocity='m/s'),
            mesh_normals='body outward',runtime_use='External reference only; no future loads/trajectories in Topwater')
        omega=2*np.pi*np.asarray(case['sea']['frequencies_hz'])
        headings=np.deg2rad(case['sea']['propagation_headings_degrees'])
        zeta=np.asarray(case['sea']['amplitudes_m'])*np.exp(1j*np.asarray(case['sea']['phases_rad']))
        mass_scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:])
        L=float(data['hull_length']);force_scale=rho*g*np.array([L**2]*3+[L**3]*3)
        out.update(angular_frequencies_rad_s=omega,headings_rad=headings,incident_phasors_m=zeta,
                   world_cg=data['world_cg'],dry_mass=data['dry_mass'],hydrostatic_stiffness=data['hydrostatic_stiffness'],
                   raw_Ainf=r['Ainf'],raw_A=r['A'],raw_B=r['B'],raw_excitation=r['excitation'],raw_RAO=r['RAO'])
        def check(name,actual,expected,scale,trace=False):
            ratio=require_close(actual,expected,scale,1e-12 if trace else 1e-11,1e-12)
            report['checks'].append(dict(name=name,maximum_tolerance_ratio=ratio))
        for index,level in enumerate(case['raw_reference']['required_levels']):
            triangles=data[f'reference_triangles_{level}'];n=len(triangles)
            body=body_from_triangles(triangles,data['world_cg'],float(data['body_mass']))
            solver=cpt.BEMSolver(green_function=green)
            prefix=f'level{level}_';out[prefix+'triangles']=triangles
            row=dict(level=level,panels=n,complete=False);report['levels'].append(row)
            fields={name:np.zeros((3,n,7),complex) for name in
                    ('sources','potential','pressure','normal_velocity','source_potential','source_normal_velocity')}
            for name,value in fields.items():out[prefix+name]=value
            combined={name:np.zeros((3,n),complex) for name in ('scattered_sources','total_hull_potential','total_hull_pressure','total_hull_normal_velocity')}
            for name,value in combined.items():out[prefix+name]=value
            modes=np.stack([body.dofs[m] for m in MODES],axis=2)
            mode_normal=np.einsum('nic,ni->nc',modes,body.mesh.faces_normals)
            for f,w in enumerate(omega):
                common=dict(body=body,omega=w,rho=rho,g=g,water_depth=np.inf,forward_speed=0.)
                for column in range(7):
                    budget()
                    problem=(cpt.RadiationProblem(**common,radiating_dof=MODES[column]) if column<6
                             else cpt.DiffractionProblem(**common,wave_direction=headings[f]))
                    result=solver.solve(problem,method='indirect',keep_details=True)
                    report['completed_physics_calls']+=1
                    budget()
                    for name in ('sources','potential','pressure'):
                        value=np.asarray(getattr(result,name))
                        if value.shape!=(n,) or not np.isfinite(value).all():
                            raise ValueError('incomplete source basis '+name)
                        fields[name][f,:,column]=value
                    fields['normal_velocity'][f,:,column]=problem.boundary_condition
                    label=f'L{level}/f{f}/c{column}'
                    check(label+'/pressure',result.pressure,1j*rho*w*result.potential,rho*w,True)
                    force=np.array([result.forces[m] for m in MODES])
                    check(label+'/integral',integrate_pressure(triangles,data['world_cg'],result.pressure),force,force_scale,True)
                    if column<6:
                        check(label+'/A',np.array([result.added_mass[m] for m in MODES]),r['A'][index,f,:,column],mass_scale[:,column])
                        check(label+'/B',np.array([result.radiation_damping[m] for m in MODES]),r['B'][index,f,:,column],w*mass_scale[:,column])
                    else:
                        check(label+'/diffraction',force,r['diffraction'][index,f,f],force_scale)
                        fk=np.array([froude_krylov_force(problem)[m] for m in MODES])
                        check(label+'/FK',fk,r['FK'][index,f,f],force_scale)
                    save(args.output,report,out)
                    print(f'source recovery: {report["completed_physics_calls"]}/84',flush=True)
                # Same pinned indirect operators as Capytaine's solve; cached for
                # these immutable mesh/frequency inputs. Check produced fields,
                # not merely a recombination of prescribed boundary conditions.
                S,K=solver.engine.build_matrices(body.mesh,body.mesh,
                    free_surface=problem.free_surface,water_depth=problem.water_depth,
                    wavenumber=problem.wavenumber,adjoint_double_layer=True,
                    diagonal_term_in_double_layer=True)
                fields['source_potential'][f]=S@fields['sources'][f]
                fields['source_normal_velocity'][f]=K@fields['sources'][f]
                del S,K
                check(f'L{level}/f{f}/source_potential',fields['source_potential'][f],fields['potential'][f],1.,True)
                check(f'L{level}/f{f}/source_normal',fields['source_normal_velocity'][f],fields['normal_velocity'][f],1.,True)
                budget()
                xi=r['RAO'][index,f,f]*zeta[f]
                incident=airy_waves_potential(body.mesh.faces_centers,problem)*zeta[f]
                incident_normal=np.einsum('ni,ni->n',airy_waves_velocity(body.mesh.faces_centers,problem),body.mesh.faces_normals)*zeta[f]
                combined['scattered_sources'][f]=combine(fields['sources'][f],zeta[f],xi)
                combined['total_hull_potential'][f]=combine(fields['source_potential'][f],zeta[f],xi)+incident
                combined['total_hull_pressure'][f]=1j*rho*w*combined['total_hull_potential'][f]
                combined['total_hull_normal_velocity'][f]=combine(fields['source_normal_velocity'][f],zeta[f],xi)+incident_normal
                check(f'L{level}/f{f}/normal',combined['total_hull_normal_velocity'][f],-1j*w*(mode_normal@xi),1.,True)
                expected=(w*w*r['A'][index,f]+1j*w*r['B'][index,f])@xi+zeta[f]*r['excitation'][index,f,f]
                check(f'L{level}/f{f}/total_load',integrate_pressure(triangles,data['world_cg'],combined['total_hull_pressure'][f]),expected,force_scale)
            row['complete']=True;save(args.output,report,out)
        if report['completed_physics_calls']!=84:raise ValueError('missing source solves')
        budget()
        report.update(complete=True,passed=True,elapsed_seconds=perf_counter()-started,
            qualification_scope='Source recovery/reproduction only; no field spatial or nonlinear restart qualification')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--controls',type=Path,default=Path('validation/periodic-source-controls.json'))
    for name in ('raw','raw-bundle','frozen','output'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();controls=json.loads(args.controls.read_text());case=load_case(controls['case'])
    report=prepare(args.output,args.controls,'periodic-source-recovery');out={}
    save(args.output,report,out)
    try:generate(args,controls,case,report,out)
    except (Exception,KeyboardInterrupt) as error:
        report.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    save(args.output,report,out)
    print(json.dumps({k:report.get(k) for k in ('complete','passed','completed_physics_calls','error')}))
    return 0 if report['complete'] and report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
