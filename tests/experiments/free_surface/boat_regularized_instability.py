"""Bounded homogeneous-mode and independent-row audit of evolving water.

This diagnoses the recorded linear instability; it does not alter the operator,
time-step controls, eigenvalues or boundary closure.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import LinearOperator,eigs,ArpackNoConvergence,splu
from scipy.fft import rfft2,irfft2

from topwater import SimulationInputs
from topwater.constants import GRAVITY
from tests.laplace_reference import triangle_samples
from .physical_p1_patch import PhysicalP1Patch
from .plane_p1_field import PlaneP1Field
from .plane_single_layer import plane_single_layer
from .cut_hat_geometry import cut_hat


def owned(path):
    report=json.loads(path.read_text())
    if not report['complete'] or sha256(path.with_suffix('.npz').read_bytes()).hexdigest()!=report['arrays_sha256']:
        raise ValueError('complete owned arrays required: '+str(path))
    return report


def run(args,report):
    water=owned(args.water);phase=owned(args.phase)
    if sha256(args.phase.read_bytes()).hexdigest()!=water['parent_phase_sha256']:
        raise ValueError('water/phase ownership changed')
    for name in ('tests/experiments/free_surface/physical_p1_patch.py',
                 'tests/experiments/free_surface/compact_p1_trace.py'):
        path=Path(name)
        if sha256(path.read_bytes()).hexdigest()!=water['producer_sha256'].get(str(path)):
            raise ValueError('recorded extension/source geometry owner changed: '+name)
    contact=json.loads((args.phase.with_suffix('.source')/'contact.json').read_text())
    if sha256((args.phase.with_suffix('.source')/'contact.json').read_bytes()).hexdigest()!=phase['parent_contact_sha256']:
        raise ValueError('phase contact owner changed')
    audit_path=Path(contact['parents']['audit']['path']);audit=owned(audit_path)
    if sha256(audit_path.read_bytes()).hexdigest()!=contact['parents']['audit']['sha256']:
        raise ValueError('audit input changed')
    inputs=SimulationInputs.from_json((audit_path.with_suffix('.source')/'inputs.json').read_text())
    if sha256(inputs.to_json().encode()).hexdigest()!=phase['input_sha256']:
        raise ValueError('physical inputs changed')
    with np.load(args.water.with_suffix('.npz')) as z:
        x=z['x'].copy();y=z['y'].copy();kernel=z['weak_full_kernel'].copy()
        interior_shape=tuple(int(v) for v in z['weak_full_shape'])
        fft_shape=tuple(int(v) for v in z['weak_fft_shape'])
        contact_rows=z['physical_contact_rows'].copy();contact_indices=z['physical_contact_active_rows'].copy()
        audit_coefficients=z['operator_audit_coefficients'].copy();audit_loads=z['operator_audit_loads'].copy()
        times=z['times'].copy();history=z['run_2_states'].copy()
    with np.load(args.phase.with_suffix('.npz')) as z:
        polygon=z['polygon'].copy();n=z['nmap'].copy();j=z['neumann'].copy()
        cross=z[f"q{phase['orders'][-1]['order']}_cross"].copy()
        active=z['active_nodes'].copy()
        archived_mass=(z['mass_data'].copy(),z['mass_indices'].copy(),z['mass_indptr'].copy())
    patch=PhysicalP1Patch(x,y,polygon);count=len(active)
    np.testing.assert_array_equal(patch.active_nodes,active)
    for current,previous in zip((patch.mass.data,patch.mass.indices,patch.mass.indptr),archived_mass,strict=True):
        np.testing.assert_array_equal(current,previous)
    if interior_shape!=(len(y)-2,len(x)-2) or contact_rows.shape!=(len(contact_indices),count):
        raise ValueError('recorded flux operator dimensions differ')
    position=np.full(len(patch.nodes),-1,dtype=int);position[patch.interior_nodes]=np.arange(len(patch.interior_nodes))
    active_positions=position[active]
    embed=np.zeros(fft_shape);dy,dx=np.indices(kernel.shape)
    embed[(dy-interior_shape[0]+1)%fft_shape[0],(dx-interior_shape[1]+1)%fft_shape[1]]=kernel
    spectrum=rfft2(embed)

    def flux_real(values):
        vector=values.ndim==1
        if vector: values=values[:,None]
        expanded=patch.expand(values)[patch.interior_nodes].reshape(*interior_shape,-1)
        padded=np.zeros((*fft_shape,values.shape[1]));padded[:interior_shape[0],:interior_shape[1]]=expanded
        result=irfft2(rfft2(padded,axes=(0,1))*spectrum[:,:,None],s=fft_shape,axes=(0,1))[
            :interior_shape[0],:interior_shape[1]].reshape(-1,values.shape[1])[active_positions].copy()
        result[contact_indices]=contact_rows@values
        return result[:,0] if vector else result

    def flux(values):
        if np.iscomplexobj(values): return flux_real(values.real)+1j*flux_real(values.imag)
        return flux_real(values)

    cpu_audit=flux(audit_coefficients)
    np.testing.assert_allclose(cpu_audit,audit_loads,rtol=1e-11,atol=1e-12)
    report['cpu_reference_vs_production_action_max']=float(abs(cpu_audit-audit_loads).max())
    if args.prior_instability is not None:
        prior=owned(args.prior_instability)
        if prior['parents']['phase']!=sha256(args.phase.read_bytes()).hexdigest():
            raise ValueError('prior mode has different phase/source-space ownership')
        with np.load(args.prior_instability.with_suffix('.npz')) as z: prior_mode=z['mode_0_chi'].copy()
        corrected=flux(prior_mode)
        report['prior_growing_mode_contact_rows']=[dict(active_row=r['active_row'],
            previous_real=r['assembled_real'],corrected_real=float(corrected[r['active_row']].real),
            corrected_imag=float(corrected[r['active_row']].imag),
            independent_q16_real=r['orders'][-1]['direct_real'],
            corrected_absolute_difference=float(abs(corrected[r['active_row']]-(r['orders'][-1]['direct_real']+1j*r['orders'][-1]['direct_imag']))))
            for r in prior['independent_rows']]
    mass=patch.mass;mf=splu(mass.tocsc());yc=mf.solve(cross);schur=j-n@yc;sf=lu_factor(schur)
    length=float(np.ptp(inputs.body.mesh.vertices[:,1]));velocity=np.sqrt(GRAVITY*length)
    scaling=np.sqrt(mass.diagonal());evaluations=0

    def real_action(value,count_work=True):
        nonlocal evaluations
        if count_work:
            evaluations+=1
            if evaluations>args.maximum_matvecs: raise RuntimeError('homogeneous eigen matvec budget exhausted')
        eta=value[:count]/scaling;psi=velocity*value[count:]/scaling
        coefficient=psi+yc@lu_solve(sf,n@psi)
        return np.r_[scaling*mf.solve(flux(coefficient)),-GRAVITY/velocity*value[:count]]

    def action(value):
        if np.iscomplexobj(value): return real_action(value.real,False)+1j*real_action(value.imag,False)
        return real_action(value,False)

    operator=LinearOperator((2*count,2*count),matvec=real_action,dtype=np.float64)
    print(json.dumps(dict(eigen_search=True,dimension=2*count,maximum_matvecs=args.maximum_matvecs)),flush=True)
    try:
        eigenvalues,vectors=eigs(operator,k=4,which='LR',ncv=40,tol=1e-8,maxiter=250,
            v0=np.sin(np.arange(2*count)*.371)+np.cos(np.arange(2*count)*.613))
        status='converged'
    except ArpackNoConvergence as error:
        eigenvalues,vectors=error.eigenvalues,error.eigenvectors;status='partial convergence'
    if not len(eigenvalues): raise RuntimeError('no resolved homogeneous modes')
    order=np.argsort(eigenvalues.real)[::-1];eigenvalues=eigenvalues[order];vectors=vectors[:,order]
    parent_contact=np.intersect1d(patch.parent_id,patch.hidden_parent_id)
    contact_piece=np.isin(patch.parent_id,parent_contact)
    mass_area=float(mass.sum());records=[];saved={}
    for index,(value,vector) in enumerate(zip(eigenvalues,vectors.T,strict=True)):
        eta=vector[:count]/scaling;psi=velocity*vector[count:]/scaling
        rms=np.sqrt(np.vdot(eta,mass@eta).real/mass_area)
        eta/=rms;psi/=rms
        sigma=-lu_solve(sf,n@psi);coefficient=psi-yc@sigma
        weak_kin=flux(coefficient)-value*(mass@eta)
        weak_dyn=-GRAVITY*(mass@eta)-value*(mass@psi)
        denominator=np.linalg.norm(value*(mass@eta))+np.linalg.norm(flux(coefficient))
        expanded=patch.expand(eta)
        field=np.einsum('tvi,ti->tv',patch.basis_values,expanded[patch.faces[patch.parent_id]])
        energy=patch.cut_areas/12*(np.sum(abs(field)**2,axis=1)+abs(field.sum(axis=1))**2)
        base=flux(psi);feedback=flux(yc@lu_solve(sf,n@psi));normalizer=np.vdot(psi,mass@psi)
        base_quotient=np.vdot(psi,base)/normalizer;feedback_quotient=np.vdot(psi,feedback)/normalizer
        row=dict(real_growth_per_s=float(value.real),oscillation_rad_per_s=float(value.imag),
            inverted_operator_residual=float(np.linalg.norm(action(vector)-value*vector)/
                (np.linalg.norm(action(vector))+abs(value)*np.linalg.norm(vector))),
            weak_kinematic_relative_residual=float(np.linalg.norm(weak_kin)/denominator),
            weak_dynamic_relative_residual=float(np.linalg.norm(weak_dyn)/
                (GRAVITY*np.linalg.norm(mass@eta)+abs(value)*np.linalg.norm(mass@psi))),
            body_normal_residual_max=float(abs(n@coefficient+j@sigma).max()),
            physical_trace_residual_max=float(abs(mass@coefficient+cross@sigma-mass@psi).max()),
            physical_contact_parent_eta_energy_fraction=float(energy[contact_piece].sum()/energy.sum()),
            physical_contact_parent_area_fraction=float(patch.cut_areas[contact_piece].sum()/patch.cut_areas.sum()),
            base_rayleigh_real=float(base_quotient.real),base_rayleigh_imag=float(base_quotient.imag),
            body_feedback_rayleigh_real=float(feedback_quotient.real),body_feedback_rayleigh_imag=float(feedback_quotient.imag),
            expected_total_rayleigh_real=float((-value*value/GRAVITY).real),
            expected_total_rayleigh_imag=float((-value*value/GRAVITY).imag))
        records.append(row);saved[f'mode_{index}_eta']=eta;saved[f'mode_{index}_psi']=psi
        saved[f'mode_{index}_chi']=coefficient;saved[f'mode_{index}_sigma']=sigma
    report.update(input_sha256=phase['input_sha256'],panel_sha256=phase['panel_sha256'],
        eigensolver_status=status,eigen_matvecs=evaluations,maximum_matvecs=args.maximum_matvecs,
        coordinate_scaling='diagonal sqrt physical mass; psi also divided by sqrt(g*beam); similarity only, no operator repair',
        modes=records,independent_rows=[])
    # Independently integrate selected rows with per-hat cuts and higher
    # target orders; this does not reuse the global contact-row assembly.
    coefficient=saved['mode_0_chi'];eta=saved['mode_0_eta'];mode_load=flux(coefficient)
    full_to_active=np.full(len(patch.nodes),-1,dtype=int);full_to_active[active]=np.arange(count)
    touched=full_to_active[patch.faces[patch.hidden_parent_id]];affected=np.unique(touched[touched>=0])
    work=abs(np.conj(eta[affected])*mode_load[affected])
    selected=affected[np.argsort(work)[-args.row_count:][::-1]]
    full=patch.expand(coefficient).reshape(len(y),len(x))
    trace=PlaneP1Field(x,y,np.stack((full.real,full.imag),axis=-1))
    pairs=0
    def potential(points):
        nonlocal pairs
        result=np.empty((len(points),2,2))
        for first in range(0,len(points),64):
            pxy=points[first:first+64];pairs+=len(pxy)*len(trace.panels)
            if pairs>args.maximum_source_pairs: raise RuntimeError('independent-row source-pair budget exhausted')
            scalar,_=plane_single_layer(pxy,trace.panels)
            result[first:first+len(pxy)]=(scalar@trace.gradient.reshape(len(trace.panels),4)).reshape(-1,2,2)
        return result
    hx,hy=float(x[1]-x[0]),float(y[1]-y[0])
    for row in selected:
        geometry=cut_hat(patch.nodes[active[row]],hx,hy,polygon);record=dict(active_row=int(row),
            full_node=int(active[row]),center=patch.nodes[active[row]].tolist(),
            physical_hat_mass=geometry['physical_hat_mass'],physical_mass_diagonal=float(mass[row,row]),
            assembled_real=float(mode_load[row].real),assembled_imag=float(mode_load[row].imag),orders=[])
        for q in (8,16):
            print(json.dumps(dict(independent_row=int(row),order=q,start=True)),flush=True)
            part=geometry['exterior'];bary,weights=triangle_samples(np.eye(3)[None],q);bary=bary[0]
            ep=np.einsum('qj,tji->tqi',bary,part['triangles']).reshape(-1,2)
            a=potential(ep);w=(part['areas'][:,None]*weights).ravel()
            value=2*np.einsum('m,mi,mic->c',w,np.repeat(part['gradients'],len(weights),axis=0),a)
            line=geometry['waterline'];nodes,weights=leggauss(q);t=(nodes+1)/2;weights=weights/2
            if len(line['segments']):
                lp=(line['segments'][:,:1]+t[None,:,None]*(line['segments'][:,1:]-line['segments'][:,:1])).reshape(-1,2)
                lv=((1-t)*line['values'][:,:1]+t*line['values'][:,1:]).ravel()
                lw=(np.linalg.norm(line['segments'][:,1]-line['segments'][:,0],axis=1)[:,None]*weights).ravel()*lv
                value+=2*np.einsum('m,mi,mic->c',lw,np.repeat(line['hole_outward_normals'],q,axis=0),potential(lp))
            direct=value[0]+1j*value[1]
            record['orders'].append(dict(order=q,direct_real=float(direct.real),direct_imag=float(direct.imag),
                absolute_difference=float(abs(direct-mode_load[row])),
                relative_difference=float(abs(direct-mode_load[row])/abs(direct)) if direct!=0 else None))
        report['independent_rows'].append(record)
    eta_norm=np.array([np.sqrt(np.sum(v[0]*(mass@v[0]),axis=0)/mass_area) for v in history])
    amplitude_norm=np.hypot(eta_norm[:,0],eta_norm[:,2]);use=amplitude_norm[-5:]>0
    report.update(observed_late_eta_growth_per_s=float(np.polyfit(times[-5:][use],np.log(amplitude_norm[-5:][use]),1)[0]),
        independent_row_source_pairs=pairs,complete=True,qualified=False,scope=__doc__,
        limits=['Computed eigenmodes and direct-row comparisons are diagnostics, not a global energy proof',
            'Fixed far chi is an unresolved nonphysical closure; localization alone cannot identify a single cause',
            'No eigenvalue/row/symmetry repair or integration-tolerance change is applied',
            'No qualified body-motion trajectory, nonlinear accuracy or realtime result'])
    np.savez(args.output.with_suffix('.npz'),eigenvalues=eigenvalues,**saved)
    report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('water','phase','output'): parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--maximum-matvecs',type=int,default=1000)
    parser.add_argument('--maximum-source-pairs',type=int,default=500_000_000)
    parser.add_argument('--row-count',type=int,choices=[1,2,3],default=2)
    parser.add_argument('--prior-instability',type=Path,help='Independent recorded counterexample on the same owned phase space')
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists(): parser.error('fresh output/source archive required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest() for name in ('phase','water')})
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
