"""Check the EXTERNAL reference's near-surface Green evaluation independently.

Keep panels fixed and distinguish tabulation, pole extraction and a new solve.
These settings are verification controls; no Topwater runtime switch is added.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.constants import GRAVITY,WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels,panel_influence
from tests.capytaine_reference import body_from_triangles,MODES
from .boat_image_plane import exterior_probes
from .boat_spatial_work import subdivide


def run(inputs,output,report):
    import capytaine as cpt
    wet=hydrostatics(inputs.body,inputs.initial_state,density=WATER_DENSITY,gravity=GRAVITY)
    triangles=subdivide(wet.wetted_triangles)
    body=body_from_triangles(triangles,inputs.initial_state.position,inputs.body.mass)
    footprint=ConvexHull(wet.waterline_segments[:,:,:2].reshape(-1,2))
    beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    probes=exterior_probes(wet.waterline_segments,footprint,beam)
    omega=2*np.pi*.1
    panels=TrianglePanels(triangles,1)
    integrated=panel_influence(probes,panels).single_layer
    radius=np.linalg.norm(probes[:,None,:]-body.mesh.faces_centers[None,:,:],axis=2)
    centroid=body.mesh.faces_areas[None,:]/(4*np.pi*radius)
    # Capytaine's negative-Rankine convention, converted to our positive S:
    # dz(Phi)-kappa*Phi = 2*kappa*(S_exact-S_centroid)*sigma_velocity.
    defect_map=2*omega**2/GRAVITY*(integrated-centroid)
    problems=[cpt.RadiationProblem(body=body,radiating_dof=mode,omega=omega,
        rho=WATER_DENSITY,g=GRAVITY,water_depth=np.inf) for mode in ['Heave','Roll']]
    cache=Path('artifacts/capytaine-surface-audit-cache');cache.mkdir(exist_ok=True)
    variants=[('default tabulated low_freq',dict(tabulation_cache_dir=str(cache))),
              ('untabulated low_freq',dict(tabulation_nr=0,tabulation_nz=0,tabulation_cache_dir=None)),
              ('untabulated extracted Rankine',dict(tabulation_nr=0,tabulation_nz=0,
                  tabulation_cache_dir=None,gf_singularities='low_freq_with_rankine_part'))]
    baseline=None;saved=dict(triangles=triangles,probes=probes);report['rows']=[]
    with threadpool_limits(1):
        for index,(name,settings) in enumerate(variants):
            before=perf_counter()
            green=cpt.Delhommeau(**settings)
            rankine,_=green.evaluate(probes,body.mesh,free_surface=np.inf,
                                     early_dot_product=False,diagonal_term_in_double_layer=False)
            capytaine_defect_map=2*omega**2/GRAVITY*(-rankine-centroid)
            solver=cpt.BEMSolver(green_function=green)
            results=[solver.solve(problem,n_threads=1) for problem in problems]
            if baseline is None:
                baseline=results
            # Existing sources isolate Green evaluation from re-solving the
            # boundary equation with that different numerical Green function.
            for source_name,solutions in [('fixed baseline sources',baseline),('new matching solve',results)]:
                value=np.column_stack([solver.compute_potential(probes,r)/(-1j*omega) for r in solutions])
                velocity=np.column_stack([solver.compute_velocity(probes,r)[:,2]/(-1j*omega) for r in solutions])
                finite=[]
                for step in [1e-3,5e-4]:
                    down=probes.copy();down[:,2]-=step
                    lower=probes.copy();lower[:,2]-=2*step
                    a=np.column_stack([solver.compute_potential(down,r)/(-1j*omega) for r in solutions])
                    b=np.column_stack([solver.compute_potential(lower,r)/(-1j*omega) for r in solutions])
                    finite.append((3*value-4*a+b)/(2*step))
                robin=omega**2/GRAVITY*value
                predicted=defect_map@np.column_stack([r.sources/(-1j*omega) for r in solutions])
                predicted_native=capytaine_defect_map@np.column_stack([r.sources/(-1j*omega) for r in solutions])
                if settings.get('gf_singularities')=='low_freq_with_rankine_part':
                    predicted=np.zeros_like(predicted)
                    predicted_native=np.zeros_like(predicted_native)
                norm=np.linalg.norm(robin,axis=0)
                row=dict(variant=name,sources=source_name,
                    green_settings=green.exportable_settings,
                    robin_max_absolute_difference=np.max(abs(velocity-robin),axis=0).tolist(),
                    robin_relative_l2=(np.linalg.norm(velocity-robin,axis=0)/norm).tolist(),
                    analytic_robin_defect_prediction_max_error=np.max(abs(velocity-robin-predicted),axis=0).tolist(),
                    native_rankine_defect_prediction_max_error=np.max(abs(velocity-robin-predicted_native),axis=0).tolist(),
                    finite_difference_relative_l2=[(np.linalg.norm(f-velocity,axis=0)/norm).tolist() for f in finite],
                    finite_difference_refinement_change=(np.linalg.norm(finite[1]-finite[0],axis=0)/norm).tolist(),
                    source_solution_added_mass=np.array([[r.added_mass[mode] for r in solutions] for mode in MODES]).tolist(),
                    source_solution_damping=np.array([[r.radiation_damping[mode] for r in solutions] for mode in MODES]).tolist())
                report['rows'].append(row)
                key=f'{index}_{source_name[:5]}'
                saved[key+'_potential']=value;saved[key+'_vertical']=velocity
                saved[key+'_predicted_defect']=predicted
                saved[key+'_finite']=np.array(finite)
                print(json.dumps(row),flush=True)
            saved[f'sources_{index}']=np.column_stack([r.sources for r in results])
            saved[f'body_{index}']=np.column_stack([r.potential/(-1j*omega) for r in results])
            report['rows'][-1]['variant_seconds']=perf_counter()-before
            np.savez(output.with_suffix('.npz'),**saved)
    report.update(complete=True,qualified=False,scope=__doc__,
        body_panels=len(triangles),panel_sha256=sha256(triangles.tobytes()).hexdigest(),
        probes=len(probes),capytaine_version=cpt.__version__,
        arrays_sha256=sha256(output.with_suffix('.npz').read_bytes()).hexdigest(),
        limits=['Reference audit only; does not validate Topwater or a trajectory',
            'Analytic prediction tests finite-panel wave-centroid versus exact Rankine integration, not continuum convergence',
            'Two one-sided finite-difference steps test the same evaluated field, not an independent continuum solution',
            'New matching solves change the external numerical reference explicitly'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output/source archive required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs,archive/'inputs.json')
    inputs=SimulationInputs.from_json(args.inputs.read_text())
    report=dict(complete=False,producer_sha256=hashes,input_sha256=sha256(inputs.to_json().encode()).hexdigest())
    started=perf_counter()
    try:
        run(inputs,args.output,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error))
        raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
