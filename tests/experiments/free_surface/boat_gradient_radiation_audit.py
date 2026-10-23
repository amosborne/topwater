"""Replace the bulk velocity evaluator at an archived physical radiation state.

This is an implementation parity gate, not new Capytaine agreement. The old
linear discretization's spatial errors remain present in both calculations.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter

import numpy as np
from numpy.testing import assert_allclose

from topwater import SimulationInputs
from topwater.bem import assemble_dirichlet_halfspace
from topwater.constants import GRAVITY, WATER_DENSITY
from .boat_radiation_time import Radiation, drive, harmonic_maps
from .cuda_source_potential import paired_velocity


def audit(parent_path, output, report):
    parent=json.loads(parent_path.read_text())
    archive=parent_path.with_suffix('.source')
    for relative, expected in parent['producer_sha256'].items():
        if sha256((archive/relative).read_bytes()).hexdigest()!=expected:
            raise ValueError('archived parent source hash mismatch: '+relative)
    arrays=parent_path.with_suffix('.npz')
    if sha256(arrays.read_bytes()).hexdigest()!=parent['arrays_sha256']:
        raise ValueError('archived parent arrays hash mismatch')
    inputs=SimulationInputs.from_json((archive/'inputs.json').read_text())
    if sha256(inputs.to_json().encode()).hexdigest()!=parent['input_sha256']:
        raise ValueError('archived parent input hash mismatch')
    # Execute the OLD producer graph in a fresh process, isolated from the
    # editable installation. Empty package markers prevent namespace fallback.
    legacy=output.with_suffix('.source')/'legacy-reference'
    shutil.copytree(archive,legacy)
    for relative in ['tests','tests/experiments','tests/experiments/free_surface']:
        (legacy/relative/'__init__.py').touch()
    reference_arrays=output.with_name(output.stem+'-legacy-rhs.npz')
    helper=Path(__file__).with_name('boat_archived_radiation_rhs.py').resolve()
    with output.with_name(output.stem+'-legacy.log').open('w') as log:
        subprocess.run([sys.executable,str(helper),str(legacy.resolve()),
                        str(parent_path.resolve()),str(reference_arrays.resolve())],
                       cwd=legacy,stdout=log,stderr=subprocess.STDOUT,check=True)
    system=Radiation(inputs,parent['grid_size'],parent['box_length_m'],
                     parent.get('body_refinement',1),parent['sponge_strength_per_second'])
    cp=system.cp
    if sha256(system.panels.triangles.tobytes()).hexdigest()!=parent['panel_sha256']:
        raise ValueError('current and archived body panels differ')
    if parent['driven_modes']!=['Heave','Roll']:
        raise ValueError('archived drive modes differ')
    assert_allclose(system.amplitudes,parent['velocity_amplitudes'],atol=0.,rtol=0.)
    with np.load(arrays) as saved:
        assert np.array_equal(saved['exterior'],system.exterior)
        assert np.array_equal(saved['hidden'],system.hidden)
        state=cp.asarray(saved['final_state'])
        archived_force=saved['forces'][-1].copy()
        time=float(saved['times'][-1])
    omega=2*np.pi*parent['frequency_hz']; ramp=parent['ramp_seconds']
    old_rate,old_force,_,_=system.rhs(time,state,omega,ramp)
    assert_allclose(cp.asnumpy(old_force),archived_force,atol=1e-8,rtol=1e-10)
    with np.load(reference_arrays) as reference:
        assert np.array_equal(reference['panels'],system.panels.triangles)
        assert np.array_equal(reference['exterior'],system.exterior)
        assert np.array_equal(reference['hidden'],system.hidden)
        legacy_rate=reference['rate'].copy()
        assert_allclose(cp.asnumpy(old_rate),legacy_rate,atol=1e-11,rtol=1e-9)
        assert_allclose(cp.asnumpy(old_force),reference['force'],atol=1e-8,rtol=1e-10)
    eta,psi=state
    spectrum=cp.fft.fft2(system.extend(psi),axes=(0,1)).reshape(-1,2)
    speed,acceleration=drive(time,omega,ramp)
    density=system.device_lu_solve(system.factor,
        system.velocity*speed-system.apply_columns(system.normal_map,spectrum).real)
    row,col=np.divmod(system.exterior,system.size)
    surface=np.column_stack((system.origin[0]+col*system.spacing,
                             system.origin[1]+row*system.spacing,np.zeros(len(row))))
    velocities=[]; diagnostics=[]
    for mode in range(2):
        result=paired_velocity(surface,system.panels,0.,density[:,mode])
        velocities.append(result.values); diagnostics.append(result.diagnostics)
    velocity=cp.stack(velocities,axis=-1)
    old_source=system.apply_columns(system.flux_map,density)
    assert_allclose(cp.asnumpy(velocity[:,:2,:]),0.,atol=1e-14,rtol=0.)
    assert_allclose(cp.asnumpy(velocity[:,2,:]),cp.asnumpy(old_source),atol=1e-11,rtol=1e-9)
    new_rate=old_rate.copy()
    new_rate[0]+=velocity[:,2,:]-old_source
    assert_allclose(cp.asnumpy(new_rate),cp.asnumpy(old_rate),atol=1e-11,rtol=1e-9)
    # Independent direct pressure-potential solve rather than the old adjoint
    # force contraction. With these FIXED linear operators only, its source
    # density also equals the motion density's derivative. The sponge term
    # retains the archived artificial pressure in the damping zone.
    psi_rate=-GRAVITY*eta-system.sponge[:,None]*psi
    rate_spectrum=cp.fft.fft2(system.extend(psi_rate),axes=(0,1)).reshape(-1,2)
    pressure_density=system.device_lu_solve(system.factor,
        system.velocity*acceleration-system.apply_columns(system.normal_map,rate_spectrum).real)
    potential,_,_=harmonic_maps(system.panels.centroids,system.panels.normals,
        system.size,system.spacing,system.origin)
    single=assemble_dirichlet_halfspace(system.panels).single_layer
    phi_rate=(cp.asarray(potential)@rate_spectrum).real+cp.asarray(single)@pressure_density
    weights=cp.asarray(system.panels.areas[:,None]*system.modes)
    direct_force=-WATER_DENSITY*weights.T@phi_rate
    assert_allclose(cp.asnumpy(direct_force),cp.asnumpy(old_force),atol=1e-8,rtol=1e-10)
    cp.cuda.get_current_stream().synchronize()
    report.update(complete=True,qualified=False,parent=str(parent_path),
        parent_sha256=sha256(parent_path.read_bytes()).hexdigest(),
        parent_arrays_sha256=parent['arrays_sha256'],input_sha256=parent['input_sha256'],
        panel_sha256=parent['panel_sha256'],grid_size=system.size,body_panels=len(system.panels),
        archived_time_seconds=time,peak_eta_m=cp.asnumpy(cp.max(abs(eta),axis=0)).tolist(),
        max_source_velocity_difference=float(cp.max(abs(velocity[:,2,:]-old_source)).item()),
        max_rhs_difference=float(cp.max(abs(new_rate-old_rate)).item()),
        max_archived_rhs_difference=float(cp.max(abs(new_rate-cp.asarray(legacy_rate))).item()),
        isolated_reference_arrays_sha256=sha256(reference_arrays.read_bytes()).hexdigest(),
        max_archived_load_difference=float(cp.max(abs(old_force-cp.asarray(archived_force))).item()),
        max_direct_pressure_load_difference=float(cp.max(abs(direct_force-old_force)).item()),
        six_dof_loads=cp.asnumpy(direct_force).tolist(),gradient_diagnostics=diagnostics,
        limits=['Linear fixed hull and flat surface only; both physical exterior states are nonzero',
                'Sponge pressure retained for parity; p=0 is not imposed in the damping zone',
                'Archived endpoint is near a velocity peak: prescribed acceleration is nearly zero',
                'Parity preserves the archived discretization, including its Capytaine discrepancy',
                'No evolving nonlinear surface, moving ownership or throughput qualification'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists():
        parser.error('fresh output and source archive required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
                   *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),
                   *Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source,target);hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes)
    started=perf_counter()
    try:
        audit(args.parent,args.output,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error))
        raise
    finally:
        report['audit_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
