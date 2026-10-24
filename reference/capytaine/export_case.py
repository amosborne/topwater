"""Minimum-mesh native Capytaine exports. Incident forcing belongs to the consumer."""
import argparse,json,shutil,subprocess,time,hashlib,sys,os
from pathlib import Path
from importlib.metadata import version,distribution,distributions
for thread_variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[thread_variable]='1'
import numpy as np
import xarray as xr
from scipy.spatial.transform import Rotation
from scipy.io import savemat
from threadpoolctl import threadpool_limits

MODES=['Surge','Sway','Heave','Roll','Pitch','Yaw']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def triangles_at(triangles,level):
    t=np.asarray(triangles)
    for _ in range(level):
        a,b,c=t.transpose(1,0,2);ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
        t=np.concatenate([np.stack(v,axis=1) for v in ((a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca))])
    # Preserve winding, put maximal vertex radius first (Capytaine3 triangle convention).
    starts=np.argmax(np.linalg.norm(t-t.mean(axis=1)[:,None,:],axis=2),axis=1)
    return np.array([np.roll(row,-int(start),axis=0) for row,start in zip(t,starts)])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--case',type=Path,default=Path('reference/capytaine/canonical-boat.json'))
    parser.add_argument('--contract',type=Path,default=Path('reference/capytaine/interface-contract.json'))
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    contract=json.loads(args.contract.read_text());controls=contract['producer_controls']
    if contract['producer_review']['status']!='ready_for_coordinator_execution_allocation':raise ValueError('Interface not reviewed')
    if controls['mesh_refinement']!=0:raise ValueError('Only minimum mesh is authorized')
    if controls['cpu_threads']!=1:raise ValueError('This minimum-case resource declaration permits one thread')
    if sha(args.case)!=contract['physical_input_sha256']:raise ValueError('Physical input identity differs from contract')
    started=time.perf_counter();out=args.output;out.mkdir(parents=True,exist_ok=False);src=out/'source';src.mkdir()
    c=json.loads(args.case.read_text());report=dict(complete=False,scope='Native Capytaine export only, no incident-forcing or motion claim',command=[sys.executable,*sys.argv],
        current_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),case_sha256=sha(args.case),mesh_level=0,contract_sha256=sha(args.contract),
        producer_sha256={},versions={name:version(name) for name in ('capytaine','numpy','scipy','xarray','threadpoolctl','pandas')},
        environment_packages={d.metadata['Name']:d.version for d in distributions()})
    for p in (Path(__file__),args.case,args.contract):
        shutil.copyfile(p,src/p.name);report['producer_sha256'][p.name]=sha(p)
    def save(): (out/'manifest.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    def budget():
        if time.perf_counter()-started>controls['maximum_elapsed_seconds']:raise TimeoutError('declared native export time exceeded')
    save()
    try:
        import capytaine as cpt
        from capytaine.io.legacy import export_hydrostatics
        from capytaine.green_functions.delhommeau import Delhommeau
        if version('capytaine')!=contract['capytaine_version']:raise ValueError('Capytaine version differs from reviewed contract')
        report['dependency_sha256']={}
        for name in report['versions']:
            package=distribution(name)
            report['dependency_sha256'][name]={str(p):sha(package.locate_file(p)) for p in package.files
                if str(p).endswith(('.py','.pyd','.so','.dll'))}
        t=triangles_at(c['minimum_wetted_triangles_world_m'],0)
        vertices,indices=np.unique(t.reshape(-1,3),axis=0,return_inverse=True)
        mesh=cpt.Mesh(vertices=vertices,faces=indices.reshape(-1,3),auto_clean=False,auto_check=False)
        body=cpt.FloatingBody(mesh=mesh,dofs=cpt.rigid_body_dofs(rotation_center=c['world_cg_m']),
            center_of_mass=c['world_cg_m'],mass=c['body_mass_kg'],name='canonical_boat')
        rho=c['water']['density_kg_m3'];g=c['water']['gravity_m_s2']
        body.inertia_matrix=xr.DataArray(c['dry_mass_matrix_world'],dims=['influenced_dof','radiating_dof'],coords=dict(influenced_dof=MODES,radiating_dof=MODES))
        body.hydrostatic_stiffness=body.compute_hydrostatic_stiffness(rho=rho,g=g)
        frequencies=np.asarray(controls['omega_rad_s']);headings=np.deg2rad(controls['propagation_headings_degrees'])
        report.update(panels=len(t),expected_calls=6+(6+len(headings))*len(frequencies),completed_calls=0,omega=frequencies.tolist(),native_thread_count=1)
        save();results=[]
        with threadpool_limits(1):
            green=Delhommeau(tabulation_cache_dir=None)
            solver=cpt.BEMSolver(green_function=green);Ainf=np.empty((6,6))
            for j,mode in enumerate(MODES):
                budget();r=solver.solve(cpt.RadiationProblem(body=body,omega=np.inf,radiating_dof=mode,rho=rho,g=g,water_depth=np.inf),keep_details=False)
                Ainf[:,j]=[r.added_mass[row] for row in MODES];report['completed_calls']+=1
            for i,w in enumerate(frequencies):
                for mode in MODES:
                    budget();r=solver.solve(cpt.RadiationProblem(body=body,omega=w,radiating_dof=mode,rho=rho,g=g,water_depth=np.inf),keep_details=False)
                    results.append(r);report['completed_calls']+=1
                for heading in headings:
                    budget();r=solver.solve(cpt.DiffractionProblem(body=body,omega=w,wave_direction=heading,rho=rho,g=g,water_depth=np.inf),keep_details=False)
                    results.append(r);report['completed_calls']+=1
                if (i+1)%20==0:save();print(f'Capytaine: {i+1}/{len(frequencies)} frequencies',flush=True)
            data=cpt.assemble_dataset(results,hydrostatics=True).sel(influenced_dof=MODES,radiating_dof=MODES)
            for key in ('added_mass','radiation_damping','diffraction_force','Froude_Krylov_force','excitation_force','inertia_matrix','hydrostatic_stiffness'):
                if key not in data or not np.isfinite(data[key]).all():raise ValueError('incomplete native '+key)
            if not np.isfinite(Ainf).all():raise ValueError('incomplete native Ainf')
            with xr.set_options(netcdf_engine_order=['scipy']):
                cpt.export_dataset(out/'canonical.nc',data)
            export_hydrostatics(str(out),body)
        np.savez(out/'native.npz',triangles=t,added_mass_infinite=Ainf,dry_mass=np.array(c['dry_mass_matrix_world']),
            hydrostatic_stiffness=body.hydrostatic_stiffness.sel(influenced_dof=MODES,radiating_dof=MODES).values,
            world_cg=np.array(c['world_cg_m']),displaced_volume=np.array(body.disp_volume),center_of_buoyancy=body.center_of_buoyancy)
        savemat(out/'native.mat',dict(Ainf=Ainf,M=np.array(c['dry_mass_matrix_world']),
            K=body.hydrostatic_stiffness.sel(influenced_dof=MODES,radiating_dof=MODES).values,cg=np.array(c['world_cg_m']),
            volume=body.disp_volume,cb=body.center_of_buoyancy,rho=rho,g=g))
        hull=c['minimum_hull'];v=np.asarray(hull['vertices_material_m']);q=hull['equilibrium_quaternion_wxyz']
        R=Rotation.from_quat([*q[1:],q[0]]).as_matrix();v=(v-np.asarray(hull['cg_material_m']))@R.T+np.asarray(c['world_cg_m'])
        with (out/'boat.stl').open('w') as f:
            f.write('solid canonical\n')
            for i,face in enumerate(hull['triangular_faces']):
                if i in hull['nonphysical_closure_faces']:continue
                tri=v[face];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
                f.write('facet normal '+' '.join(map(str,n))+'\nouter loop\n')
                for p in tri:f.write('vertex '+' '.join(map(str,p))+'\n')
                f.write('endloop\nendfacet\n')
            f.write('endsolid canonical\n')
        budget();report.update(complete=True,elapsed_seconds=time.perf_counter()-started,displaced_volume_m3=float(body.disp_volume),
            coefficient_thresholds_blocking=False,output_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'})
    except BaseException as error:report.update(error=dict(type=type(error).__name__,message=str(error)),elapsed_seconds=time.perf_counter()-started)
    save();print(json.dumps({k:report.get(k) for k in ('complete','panels','completed_calls','elapsed_seconds','error')}));return 0 if report['complete'] else 1


if __name__=='__main__':raise SystemExit(main())
