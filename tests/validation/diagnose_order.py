"""Bounded same-surface vertex-order diagnostic, never qualification evidence."""
import argparse
from importlib import import_module
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from threadpoolctl import threadpool_limits,threadpool_info
from .common import prepare,load_contract,save,digest
from .__main__ import frozen_input


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--contract',type=Path,default=Path('validation/canonical-boat.json'))
    parser.add_argument('--frozen',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();contract=load_contract(args.contract)
    report=prepare(args.output,args.contract,'vertex-order-diagnostic');arrays={}
    started=perf_counter()
    try:
        frozen,data=frozen_input(args,contract,report)
        import capytaine as cpt
        from tests.capytaine_reference import body_from_triangles,MODES
        import_module('capytaine.green_functions.Delhommeau_float64')
        package=Path(cpt.__file__).parent
        report['inspected_dependency_sources']={}
        for name in ('meshes/geometry.py','meshes/meshes.py','green_functions/delhommeau.py'):
            path=package/name;target=args.output.with_suffix('.source')/'dependency-sources'/name
            target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
            report['inspected_dependency_sources'][name]=digest(path)
        level=contract['reference_numerics']['body_subdivision_levels'][0]
        triangles=data[f'reference_triangles_{level}']
        centers=triangles.mean(axis=1)
        first=np.argmax(np.linalg.norm(triangles-centers[:,None],axis=2),axis=1)
        farthest=np.take_along_axis(triangles,((first[:,None]+np.arange(3))%3)[:,:,None],axis=1)
        variants=dict(original=triangles,cyclic1=np.roll(triangles,-1,axis=1),
            cyclic2=np.roll(triangles,-2,axis=1),farthest_first=farthest)
        arrays['original_triangles']=triangles
        report.update(panels=len(triangles),level=level,capytaine_version=cpt.__version__,
            scope='All36 Ainf and A/B at 0.1Hz; only vertex cyclic order varies; no qualification',
            variants=[],completed_physics_calls=0)
        scale=np.sqrt(data['dry_mass'].diagonal()[:,None]*data['dry_mass'].diagonal()[None,:])
        omega=2*np.pi*.1;rho=contract['physical_problem']['density_kg_m3'];g=contract['physical_problem']['gravity_m_s2']
        with threadpool_limits(1):
            report['native_thread_pools']=threadpool_info();save(args.output,report,arrays)
            inf=cpt.Delhommeau(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None)
            finite=cpt.Delhommeau(tabulation_cache_dir=None)
            for name,t in variants.items():
                body=body_from_triangles(t,data['world_cg'],float(data['body_mass']))
                np.testing.assert_allclose(body.mesh.faces_centers,centers,rtol=0,atol=1e-14)
                if name=='original':
                    arrays['normals']=body.mesh.faces_normals
                    arrays['areas']=body.mesh.faces_areas
                else:
                    np.testing.assert_allclose(body.mesh.faces_normals,arrays['normals'],rtol=1e-10,atol=1e-12)
                    np.testing.assert_allclose(body.mesh.faces_areas,arrays['areas'],rtol=1e-10,atol=1e-12)
                arrays[name+'_radii']=body.mesh.faces_radiuses
                arrays[name+'_triangles']=t
                report['variants'].append(dict(name=name,radius_min=float(body.mesh.faces_radiuses.min()),
                    radius_max=float(body.mesh.faces_radiuses.max())))
                for label,w,green in (('Ainf',np.inf,inf),('A',omega,finite)):
                    A=np.zeros((6,6));B=np.zeros((6,6));solver=cpt.BEMSolver(green_function=green)
                    for j,mode in enumerate(MODES):
                        if perf_counter()-started>contract['reference_numerics']['maximum_elapsed_seconds']:
                            raise TimeoutError('diagnostic elapsed limit; no skipped-call success')
                        result=solver.solve(cpt.RadiationProblem(body=body,radiating_dof=mode,
                            omega=w,rho=rho,g=g,water_depth=np.inf))
                        A[:,j]=[result.added_mass[row] for row in MODES]
                        B[:,j]=[result.radiation_damping[row] for row in MODES]
                        report['completed_physics_calls']+=1
                    arrays[name+'_'+label]=A/scale
                    if label=='A':arrays[name+'_B']=B/omega/scale
                    save(args.output,report,arrays)
        comparisons={}
        for quantity in ('Ainf','A','B'):
            reference=arrays['farthest_first_'+quantity]
            allowance=1e-6+.001*abs(reference)
            comparisons[quantity]={name:float(np.max(abs(arrays[name+'_'+quantity]-reference)/allowance))
                for name in variants if name!='farthest_first'}
        report.update(complete=True,passed=all(value<=1 for row in comparisons.values() for value in row.values()),
            comparisons_maximum_reference_budget_ratio=comparisons,
            qualification_status='unqualified; diagnostic only')
    except (Exception,KeyboardInterrupt) as error:
        report.update(complete=False,passed=False,error=dict(type=type(error).__name__,message=str(error)))
    report['elapsed_seconds']=perf_counter()-started;save(args.output,report,arrays)
    print(__import__('json').dumps(dict(complete=report['complete'],passed=report['passed'],output=str(args.output))))
    return 0 if report['complete'] and report['passed'] else 1


if __name__=='__main__':
    raise SystemExit(main())
