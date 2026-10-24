"""Declared offline Capytaine study; never a runtime force model."""
from time import perf_counter
import numpy as np
from threadpoolctl import threadpool_limits
from .common import REFERENCE,current_sources,hashes,mesh_quality,save,study_preflight,triangles_at


def normalized(raw,mass,length,rho,gravity,omega):
    scale=np.sqrt(mass.diagonal()[:,None]*mass.diagonal()[None,:])
    force=rho*gravity*np.array([length**2]*3+[length**3]*3)
    return dict(Ainf=raw['Ainf']/scale,A=raw['A']/scale,
        B=raw['B']/omega[None,:,None,None]/scale,
        diffraction=raw['diffraction']/force,FK=raw['FK']/force,
        excitation=raw['excitation']/force,
        RAO=raw['RAO']*np.array([1.,1.,1.,length,length,length]))


def require_qualified(report,arrays,contract,frozen):
    levels=contract['reference_numerics']['body_subdivision_levels']
    nf,nh=len(contract['cases']['frequencies_hz']),len(contract['cases']['propagation_headings_rad'])
    shapes=dict(Ainf=(len(levels),6,6),A=(len(levels),nf,6,6),B=(len(levels),nf,6,6),
        **{name:(len(levels),nf,nh,6) for name in ('diffraction','FK','excitation','RAO')})
    if ([row['level'] for row in report['levels']]!=levels
            or not all(row['complete'] is True for row in report['levels'])
            or report['capytaine_version']!=contract['reference_numerics']['capytaine_version']):
        raise ValueError('missing reference levels, cases or pinned version')
    for name,shape in shapes.items():
        if arrays[name].shape!=shape or not np.isfinite(arrays[name]).all():
            raise ValueError('missing complete raw reference entries: '+name)
    for row,level in zip(report['levels'],levels,strict=True):
        expected=triangles_at(frozen['wetted_triangles'],level)
        if row['panels']!=len(expected):
            raise ValueError('reference panel count differs from declared nested physical mesh')
        np.testing.assert_array_equal(expected,arrays[f'triangles_{level}'])
    p=contract['physical_problem'];omega=2*np.pi*np.asarray(contract['cases']['frequencies_hz'])
    values=normalized(arrays,frozen['dry_mass'],float(frozen['hull_length']),p['density_kg_m3'],p['gravity_m_s2'],omega)
    for name,value in values.items():
        uncertainty,gate=mesh_quality(value,contract['acceptance'])
        if not gate['passed'] or report['gates'][name]['passed'] is not True:
            raise ValueError('unresolved reference mesh uncertainty: '+name)
        np.testing.assert_array_equal(value,arrays['normalized_'+name])
        np.testing.assert_array_equal(uncertainty,arrays['uncertainty_'+name])


def generate(args,contract,frozen,data,report,arrays):
    started=perf_counter()
    import capytaine as cpt
    from tests.capytaine_reference import body_from_triangles,MODES
    from capytaine.bem.airy_waves import froude_krylov_force
    spec=contract['reference_numerics'];physical=contract['physical_problem']
    if cpt.__version__!=spec['capytaine_version']:
        raise ValueError('required pinned Capytaine dependency/version unavailable')
    levels=spec['body_subdivision_levels']
    plan=study_preflight(len(data['wetted_triangles']),levels,spec['maximum_estimated_workspace_bytes'])
    rho,g=physical['density_kg_m3'],physical['gravity_m_s2']
    omega=2*np.pi*np.asarray(contract['cases']['frequencies_hz'])
    headings=contract['cases']['propagation_headings_rad']

    def budget():
        if perf_counter()-started>spec['maximum_elapsed_seconds']:
            raise TimeoutError('reference study elapsed budget exhausted between physics calls; no skipped-case pass')

    report.update(reference_source_sha256=hashes(REFERENCE),capytaine_version=cpt.__version__,
        study_plan=plan,levels=[],completed_physics_calls=0,
        reference_scope='Raw all36 Ainf and A/B, Cartesian complex excitation/RAO; measured mesh uncertainty, no Richardson')
    raw={name:[] for name in ('Ainf','A','B','diffraction','FK','excitation','RAO')}
    inf_green=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None)
    finite_green=cpt.Delhommeau(tabulation_cache_dir=None)
    report['green_settings']=dict(infinite=inf_green.exportable_settings,finite=finite_green.exportable_settings)
    for level in levels:
        budget();triangles=triangles_at(data['wetted_triangles'],level)
        body=body_from_triangles(triangles,data['world_cg'],float(data['body_mass']))
        row=dict(level=level,panels=len(triangles),complete=False)
        report['levels'].append(row);arrays[f'triangles_{level}']=triangles
        Ainf=np.zeros((6,6));A=np.zeros((len(omega),6,6));B=np.zeros_like(A)
        D=np.zeros((len(omega),len(headings),6),complex);FK=np.zeros_like(D)
        arrays.update({f'level{level}_Ainf':Ainf,f'level{level}_A':A,f'level{level}_B':B,
                       f'level{level}_diffraction':D,f'level{level}_FK':FK})
        save(args.output,report,arrays)
        with threadpool_limits(1):
            solver=cpt.BEMSolver(green_function=inf_green)
            for column,mode in enumerate(MODES):
                budget();result=solver.solve(cpt.RadiationProblem(body=body,radiating_dof=mode,
                    omega=np.inf,rho=rho,g=g,water_depth=np.inf))
                Ainf[:,column]=[result.added_mass[name] for name in MODES]
                if any(result.radiation_damping[name]!=0 for name in MODES):
                    raise ValueError('infinite-frequency reference returned nonzero damping')
                if not np.isfinite(Ainf[:,column]).all():
                    raise FloatingPointError('nonfinite raw Ainf column')
                report['completed_physics_calls']+=1;save(args.output,report,arrays)
            solver=cpt.BEMSolver(green_function=finite_green)
            for f,w in enumerate(omega):
                common=dict(body=body,omega=w,rho=rho,g=g,water_depth=np.inf)
                for column,mode in enumerate(MODES):
                    budget();result=solver.solve(cpt.RadiationProblem(**common,radiating_dof=mode))
                    A[f,:,column]=[result.added_mass[name] for name in MODES]
                    B[f,:,column]=[result.radiation_damping[name] for name in MODES]
                    if not np.isfinite([A[f,:,column],B[f,:,column]]).all():
                        raise FloatingPointError('nonfinite finite-frequency radiation column')
                    report['completed_physics_calls']+=1;save(args.output,report,arrays)
                for h,heading in enumerate(headings):
                    budget();problem=cpt.DiffractionProblem(**common,wave_direction=heading)
                    result=solver.solve(problem);incident=froude_krylov_force(problem)
                    D[f,h]=[result.forces[name] for name in MODES];FK[f,h]=[incident[name] for name in MODES]
                    if not np.isfinite([D[f,h],FK[f,h]]).all():
                        raise FloatingPointError('nonfinite complex excitation')
                    report['completed_physics_calls']+=1;save(args.output,report,arrays)
        X=D+FK;rao=np.empty_like(X)
        for f,w in enumerate(omega):
            impedance=-w*w*(data['dry_mass']+A[f])-1j*w*B[f]+data['hydrostatic_stiffness']
            rao[f]=np.linalg.solve(impedance,X[f].T).T
        for name,value in dict(Ainf=Ainf,A=A,B=B,diffraction=D,FK=FK,excitation=X,RAO=rao).items():
            if not np.isfinite(value).all():
                raise FloatingPointError('nonfinite coupled reference '+name)
            raw[name].append(value.copy())
        arrays[f'level{level}_RAO']=rao;arrays[f'level{level}_excitation']=X
        row.update(complete=True,raw_Ainf_reciprocity_max=float(abs(Ainf-Ainf.T).max()))
        save(args.output,report,arrays);budget()
    raw={name:np.stack(values) for name,values in raw.items()}
    norm=normalized(raw,data['dry_mass'],float(data['hull_length']),rho,g,omega)
    gates={}
    for name,value in norm.items():
        uncertainty,gates[name]=mesh_quality(value,contract['acceptance'])
        arrays['normalized_'+name]=value;arrays['uncertainty_'+name]=uncertainty
    arrays.update(raw);arrays['linear_motion_amplitude']=raw['RAO']*physical['incident_elevation_amplitude_m']
    report.update(complete=True,passed=all(row['passed'] for row in gates.values()),
        gates=gates,elapsed_seconds=perf_counter()-started,
        uncertainty_scope='Componentwise envelope of last two normalized increments; empirical, not a certified bound',
        qualification_scope=contract['acceptance']['scope'])
    save(args.output,report,arrays)
