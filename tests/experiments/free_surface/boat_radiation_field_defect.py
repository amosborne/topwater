"""External harmonic snapshots isolate spatial error before time integration.

Capytaine supplies unit-velocity complex amplitudes. Re/Im are two independent
real time phases passed through the existing REAL flat time-domain BVP. No
frequency-domain solve or external load model is added to Topwater execution.
Full auxiliary data and physical-exterior harmonic fill are separate controls.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from numpy.testing import assert_allclose
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.hydrostatics import hydrostatics
from topwater.rankine import TrianglePanels
from tests.capytaine_reference import body_from_triangles, MODES
from tests.benchmark_rankine import peak_rss_bytes
from .boat_radiation_time import Radiation, harmonic_maps
from .boat_spatial_work import subdivide


def amplitude_norm(value, weights):
    return np.sqrt(np.sum(weights[:,None]*abs(value)**2,axis=0))


def run(inputs,reference,grids,refinement,output,report):
    import capytaine as cpt
    triangles=hydrostatics(inputs.body,inputs.initial_state).wetted_triangles
    for _ in range(refinement):
        triangles=subdivide(triangles)
    panels=TrianglePanels(triangles,-1)
    panel_hash=sha256(triangles.tobytes()).hexdigest()
    omega=2*np.pi*.1
    matching=[row['reference'] for row in reference['rows']
        if row['panel_sha256']==panel_hash
        and row['reference']['angular_frequency']!='infinite'
        and np.isclose(row['reference']['angular_frequency'],omega,rtol=0.,atol=1e-14)]
    if len(matching)!=1:
        raise ValueError('one exact-panel finite-frequency reference required')
    expected=matching[0]
    if (expected['density']!=WATER_DENSITY or expected['gravity']!=GRAVITY
            or expected['water_depth']!='infinite'
            or expected['radiating_modes']!=['Heave','Roll']
            or expected['influenced_modes']!=list(MODES)):
        raise ValueError('reference physical/force conventions differ')
    body=body_from_triangles(triangles,inputs.initial_state.position,inputs.body.mass)
    green=cpt.Delhommeau(tabulation_cache_dir=None)
    solutions=[]; added=[]; damping=[]
    before=perf_counter()
    with threadpool_limits(1):
        solver=cpt.BEMSolver(green_function=green)
        for mode in ['Heave','Roll']:
            result=solver.solve(cpt.RadiationProblem(body=body,radiating_dof=mode,
                omega=omega,rho=WATER_DENSITY,g=GRAVITY,water_depth=np.inf),n_threads=1)
            solutions.append(result)
            added.append([result.added_mass[name] for name in MODES])
            damping.append([result.radiation_damping[name] for name in MODES])
    added=np.array(added).T;damping=np.array(damping).T
    assert_allclose(added,expected['added_mass'],atol=1e-8,rtol=1e-10)
    assert_allclose(damping,expected['radiation_damping'],atol=1e-8,rtol=1e-10)
    reference_body=np.column_stack([r.potential/(-1j*omega) for r in solutions])
    modes=np.column_stack((panels.normals,
        np.cross(panels.centroids-inputs.initial_state.position,panels.normals)))
    weights=panels.areas[:,None]*modes
    reference_force=-damping+1j*omega*added
    assert_allclose(1j*WATER_DENSITY*omega*weights.T@reference_body,
                    reference_force,atol=1e-8,rtol=1e-10)
    report.update(panel_sha256=panel_hash,body_panels=len(panels),
        frequency_hz=.1,capytaine_version=cpt.__version__,
        capytaine_solve_seconds=perf_counter()-before,
        capytaine_added_mass=added.tolist(),capytaine_damping=damping.tolist(),rows=[])
    saved=dict(reference_body=reference_body,reference_force=reference_force)
    for case,(length,size) in enumerate(grids):
        before=perf_counter()
        system=Radiation(inputs,size,length,refinement,0.)
        cp=system.cp
        x,y=np.meshgrid(system.origin[0]+np.arange(size)*system.spacing,
                        system.origin[1]+np.arange(size)*system.spacing)
        surface=np.column_stack((x.ravel(),y.ravel(),np.zeros(size*size)))
        reference_surface=np.empty((size*size,2),dtype=complex)
        with threadpool_limits(1):
            for column,result in enumerate(solutions):
                for first in range(0,len(surface),256):
                    last=min(first+256,len(surface))
                    reference_surface[first:last,column]=solver.compute_potential(surface[first:last],result)/(-1j*omega)
        distance=np.linalg.norm(surface[system.exterior,:2]-inputs.initial_state.position[:2],axis=1)
        near=system.exterior[np.argsort(distance)[:12]]
        far=system.exterior[np.linspace(0,len(system.exterior)-1,12,dtype=int)]
        probes=np.unique(np.r_[near,far])
        reference_vertical=np.empty((len(probes),2),dtype=complex)
        with threadpool_limits(1):
            for column,result in enumerate(solutions):
                reference_vertical[:,column]=solver.compute_velocity(surface[probes],result)[:,2]/(-1j*omega)
        robin=omega**2/GRAVITY*reference_surface[probes]
        assert_allclose(reference_vertical,robin,atol=1e-8,rtol=1e-5)
        robin_check=dict(points=len(probes),
            max_absolute_velocity_difference=np.max(abs(reference_vertical-robin),axis=0).tolist(),
            relative_l2=(np.linalg.norm(reference_vertical-robin,axis=0)/np.linalg.norm(robin,axis=0)).tolist())
        cp.cuda.get_current_stream().synchronize()
        preparation=perf_counter()-before
        potential,_,_=harmonic_maps(panels.centroids,panels.normals,size,system.spacing,system.origin)
        potential=cp.asarray(potential)
        single=cp.asarray(assemble_dirichlet_halfspace(panels).single_layer)
        q=cp.asarray(modes[:,[2,3]])
        # Linear methods take real data. The two pairs of columns below are
        # the t=0 and t=period/4 phases of Re(Phi exp(-i*omega*t)).
        raw=np.column_stack((reference_surface.real,reference_surface.imag))
        full=cp.asarray(raw)
        filled=system.extend(full[system.device_exterior]).reshape(size*size,4)
        for label,trace in [('full reference auxiliary trace',full),('physical exterior harmonic fill',filled)]:
            started=perf_counter()
            spectrum=cp.fft.fft2(trace.reshape(size,size,4),axes=(0,1)).reshape(-1,4)
            right=cp.column_stack((q,cp.zeros_like(q)))-system.apply_columns(system.normal_map,spectrum).real
            density=system.device_lu_solve(system.factor,right)
            body_value=(potential@spectrum).real+single@density
            flat_flux=cp.fft.ifft2((spectrum.reshape(size,size,4)*system.wave[:,:,None]),
                axes=(0,1)).real.reshape(-1,4)[system.device_exterior]
            source_flux=system.apply_columns(system.flux_map,density)
            flux=flat_flux+source_flux
            computed_body=cp.asnumpy(body_value[:,:2]+1j*body_value[:,2:])
            computed_flux=cp.asnumpy(flux[:,:2]+1j*flux[:,2:])
            computed_lift=cp.asnumpy(flat_flux[:,:2]+1j*flat_flux[:,2:])
            computed_source=cp.asnumpy(source_flux[:,:2]+1j*source_flux[:,2:])
            exact_flux=omega**2/GRAVITY*reference_surface[system.exterior]
            force=1j*WATER_DENSITY*omega*weights.T@computed_body
            # Do not divide cross loads by their nearly-zero true amplitudes.
            # Each driven column is normalized by its principal reactive load.
            force_scale=abs(reference_force[[2,3],[0,1]])
            field_error=amplitude_norm(computed_body-reference_body,panels.areas)/amplitude_norm(reference_body,panels.areas)
            zones=[]
            for low,high in [(0.,5.),(5.,15.),(15.,30.)]:
                selected=(distance>=low)&(distance<high)
                if np.any(selected):
                    error=computed_flux[selected]-exact_flux[selected]
                    zones.append(dict(radius_m=[low,high],nodes=int(selected.sum()),
                        relative_l2=(amplitude_norm(error,np.ones(selected.sum()))
                            /amplitude_norm(exact_flux[selected],np.ones(selected.sum()))).tolist(),
                        lift_max_velocity_amplitude=np.max(abs(computed_lift[selected]),axis=0).tolist(),
                        source_max_velocity_amplitude=np.max(abs(computed_source[selected]),axis=0).tolist(),
                        exact_max_velocity_amplitude=np.max(abs(exact_flux[selected]),axis=0).tolist(),
                        max_absolute_velocity_error_m_per_s=np.max(abs(error),axis=0).tolist()))
            row=dict(length_m=length,grid_size=size,spacing_m=system.spacing,
                trace_control=label,body_potential_relative_l2=field_error.tolist(),
                capytaine_vertical_robin_check=robin_check,
                normalized_six_dof_load_error=(abs(force-reference_force)/force_scale).tolist(),
                principal_added_mass_from_snapshot=(force.imag/omega)[[2,3],[0,1]].tolist(),
                principal_damping_from_snapshot=(-force.real)[[2,3],[0,1]].tolist(),
                near_surface_flux=zones,preparation_seconds=preparation,
                application_seconds=perf_counter()-started)
            report['rows'].append(row)
            print(json.dumps(row),flush=True)
            saved[f'body_{case}_{label[:4]}']=computed_body
            saved[f'flux_{case}_{label[:4]}']=computed_flux
        saved[f'surface_{case}']=reference_surface
        saved[f'exterior_{case}']=system.exterior
        del system,potential,single,full,filled,spectrum,body_value,flux,density
        cp.get_default_memory_pool().free_all_blocks()
        np.savez(output.with_suffix('.npz'),**saved)
    report.update(complete=True,qualified=False,peak_process_rss_bytes=peak_rss_bytes(),
        arrays_sha256=sha256(output.with_suffix('.npz').read_bytes()).hexdigest(),
        scope=__doc__,limits=['Two real time phases of an external harmonic field, no Topwater frequency-domain runtime',
            'No trajectory or nonlinear moving geometry; both phase amplitudes are needed for errors',
            'Periodic outer data and physical hidden fill both differ from the infinite-domain reference',
            'Surface norms are sampled near-field defects, not certified continuum error bounds'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--reference',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--grid',type=int,nargs='+',default=[129,257])
    parser.add_argument('--length',type=float,nargs='+',default=[128.,128.])
    parser.add_argument('--refinement',type=int,default=1)
    args=parser.parse_args()
    if (args.output.exists() or args.output.with_suffix('.source').exists()
            or len(args.grid)!=len(args.length) or not args.grid or args.refinement not in (0,1,2)
            or any(n<17 or n%2!=1 for n in args.grid)
            or not np.isfinite(args.length).all() or min(args.length)<=0):
        parser.error('fresh output, paired positive lengths/odd grids and refinement0..2 required')
    inputs=SimulationInputs.from_json(args.inputs.read_text())
    reference=json.loads(args.reference.read_text())
    if sha256(inputs.to_json().encode()).hexdigest()!=reference['input_sha256']:
        parser.error('reference and input identity differ')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.inputs,archive/'inputs.json')
    shutil.copyfile(args.reference,archive/'reference.json')
    report=dict(complete=False,producer_sha256=hashes,
        input_sha256=sha256(inputs.to_json().encode()).hexdigest(),
        reference_sha256=sha256(args.reference.read_bytes()).hexdigest())
    started=perf_counter()
    try:
        run(inputs,reference,list(zip(args.length,args.grid)),args.refinement,args.output,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error))
        raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
