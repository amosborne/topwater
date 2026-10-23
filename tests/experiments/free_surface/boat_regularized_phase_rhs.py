"""Physical-psi weak coupling and pressure at two external radiation phases.

This is a fixed-hull instantaneous RHS falsifier with prescribed far-field
data. Physical eta moments come from the independent external potential.
It does NOT advance eta with the kinematic equation or qualify a trajectory.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import lu_factor,lu_solve
from scipy.sparse.linalg import splu
from threadpoolctl import threadpool_limits

from topwater import SimulationInputs
from topwater.constants import GRAVITY,WATER_DENSITY
from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels,panel_influence,single_layer_trace
from tests.benchmark_rankine import peak_rss_bytes
from tests.capytaine_reference import body_from_triangles
from tests.laplace_reference import triangle_samples
from .boat_local_radiation_trace import phases,amplitude,far_trace_spectrum
from .boat_local_trace_enrichment import tensor_values
from .boat_radiation_time import harmonic_maps
from .boat_spatial_work import source_matrix
from .boat_weak_contact_rows import owned,triangle_rule,fourier_values
from .compact_p1_maps import body_maps
from .compact_p1_weak import hat_fourier
from .cut_hat_geometry import cut_hat
from .physical_p1_patch import PhysicalP1Patch
from .plane_p1_field import PlaneP1Field
from .panel_gradient_reference import checked_influence


def cut_quadrature(patch,order,hidden=False):
    triangles=patch.hidden_triangles if hidden else patch.cut_triangles
    basis=patch.hidden_basis_values if hidden else patch.basis_values
    parents=patch.hidden_parent_id if hidden else patch.parent_id
    areas=patch.hidden_areas if hidden else patch.cut_areas
    bary,weights=triangle_samples(np.eye(3)[None],order);bary=bary[0]
    for first in range(0,len(triangles),32):
        last=min(first+32,len(triangles));tri=triangles[first:last]
        xy=np.einsum('qj,tji->tqi',bary,tri)
        values=np.einsum('qj,tji->tqi',bary,basis[first:last])
        weighted=areas[first:last,None,None]*weights[None,:,None]*values
        points=np.column_stack((xy.reshape(-1,2),np.zeros(len(tri)*len(weights))))
        nodes=np.repeat(patch.faces[parents[first:last]],len(weights),axis=0)
        yield points,nodes,weighted.reshape(-1,3)


def scatter_moments(result,nodes,weights,field):
    for local in range(3):
        np.add.at(result,nodes[:,local],weights[:,local,None]*field)


def run(args,report):
    import capytaine as cpt
    contact=owned(args.contact);paths={k:Path(v['path']) for k,v in contact['parents'].items()}
    metadata={}
    for name,path in paths.items():
        if sha256(path.read_bytes()).hexdigest()!=contact['parents'][name]['sha256']:
            raise ValueError('changed owned contact input: '+name)
        metadata[name]=owned(path)
    snapshot,field,audit=(metadata[k] for k in ('snapshot','field','audit'))
    inputs=SimulationInputs.from_json((paths['audit'].with_suffix('.source')/'inputs.json').read_text())
    identity=sha256(inputs.to_json().encode()).hexdigest()
    if any(r['input_sha256']!=identity for r in (contact,snapshot,field,audit)):
        raise ValueError('physical input ownership differs')
    with np.load(paths['snapshot'].with_suffix('.npz')) as z:
        x=z['patch_x'].copy();y=z['patch_y'].copy()
    stride=2**args.coarsening
    if (len(x)-1)%stride or (len(y)-1)%stride: raise ValueError('nested source axes required')
    x=x[::stride];y=y[::stride]
    with np.load(args.contact.with_suffix('.npz')) as z: polygon=z['polygon'].copy()
    patch=PhysicalP1Patch(x,y,polygon)
    with np.load(paths['audit'].with_suffix('.npz')) as z:
        triangles=z['triangles'].copy();external_sources=z['sources_0'].copy();reference_body=z['body_0'].copy()
    panel_hash=sha256(triangles.tobytes()).hexdigest()
    if any(r['panel_sha256']!=panel_hash for r in (contact,snapshot,field,audit)):
        raise ValueError('physical body panel ownership differs')
    panels=TrianglePanels(triangles,-1);nb=len(panels);nn=len(patch.nodes)
    if 8*nb*nn*4>LaplaceControls().max_operator_bytes:
        raise MemoryError('phase-gate retained cross maps exceed operator byte budget')
    query_count=sum((len(patch.cut_triangles)+len(patch.hidden_triangles))*q*q for q in args.orders)
    if query_count>args.maximum_quadrature_points:
        raise RuntimeError('physical quadrature query budget exhausted before construction')
    omega=2*np.pi*field['frequency_hz'];cg=inputs.initial_state.position
    grid=snapshot['grid'];length=snapshot['length_m'];origin=cg[:2]-length/2
    beam=float(np.ptp(inputs.body.mesh.vertices[:,1]))
    inner=np.max(abs(triangles[:,:,:2]-cg[:2]),axis=(0,1))+.5*beam
    with np.load(paths['field'].with_suffix('.npz')) as z:
        trace=z['surface_0'].reshape(grid,grid,2);reference_force=z['reference_force'].copy()
    external_density=phases(-external_sources/(-1j*omega))
    body=body_from_triangles(triangles,cg,inputs.body.mass)
    settings=dict(audit['rows'][0]['green_settings']);settings.pop('green_function')
    cache=Path('artifacts/capytaine-surface-audit-cache');cache.mkdir(exist_ok=True)
    solver=cpt.BEMSolver(green_function=cpt.Delhommeau(**settings,tabulation_cache_dir=str(cache)))
    external=[SimpleNamespace(body=body,sources=external_sources[:,c],free_surface=0.,water_depth=np.inf,
        encounter_wavenumber=omega**2/GRAVITY,forward_speed=0.) for c in range(2)]

    def reference(points,velocity=False):
        values=[];function=solver.compute_velocity if velocity else solver.compute_potential
        with threadpool_limits(1):
            for result in external:
                column=np.concatenate([function(points[first:first+256],result)
                    for first in range(0,len(points),256)])
                values.append((column[:,2] if velocity else column)/(-1j*omega))
        return np.column_stack(values)

    def source_trace(points):
        result=np.empty((len(points),nb))
        for first in range(0,len(points),128):
            result[first:first+128]=2*panel_influence(points[first:first+128],panels).single_layer
        return result

    def rate(values): return np.column_stack((omega*values[:,2:],-omega*values[:,:2]))

    started=perf_counter()
    axis=np.arange(grid)*length/grid;gx,gy=np.meshgrid(origin[0]+axis,origin[1]+axis)
    gp=np.column_stack((gx.ravel(),gy.ravel(),np.zeros(grid*grid)))
    trace=trace-amplitude(source_trace(gp)@external_density).reshape(grid,grid,2)
    far=far_trace_spectrum(trace,origin,length,cg[:2],inner,beam)
    pfar,nfar,_=harmonic_maps(panels.centroids,panels.normals,grid,length/grid,origin)
    far_body=(pfar@far.reshape(-1,4)).real;far_normal=(nfar@far.reshape(-1,4)).real
    full_p,full_n=body_maps(patch,panels.centroids,panels.normals)
    pmap=patch.condense_maps(full_p);nmap=patch.condense_maps(full_n)
    even=np.empty((nb,nb));neumann=source_matrix(panels)
    for first in range(0,nb,128):
        last=min(first+128,nb);indices=np.arange(first,last)
        image=checked_influence(panels.centroids[first:last]*[1.,1.,-1.],panels)
        even[first:last]=single_layer_trace(panels,indices).single_layer+image.single_layer
        neumann[first:last]+=2*np.einsum('mni,mi->mn',image.gradient_single_layer,
            panels.normals[first:last]*[1.,1.,-1.])
    nfactor=lu_factor(neumann)
    modes=np.column_stack((panels.normals,np.cross(panels.centroids-cg,panels.normals)))
    q=np.column_stack((modes[:,[2,3]],np.zeros((nb,2))))
    weighted_modes=panels.areas[:,None]*modes
    hx,hy=float(x[1]-x[0]),float(y[1]-y[0])
    kx,ky=np.meshgrid(2*np.pi*np.fft.fftfreq(grid,d=length/grid),2*np.pi*np.fft.fftfreq(grid,d=length/grid))
    full_far=np.zeros((nn,4))
    full_far.reshape(len(y),len(x),4)[1:-1,1:-1]=tensor_values(
        far*hat_fourier(hx,hy,kx,ky)[:,:,None],x[1:-1],y[1:-1],origin,length)
    mass_factor=splu(patch.mass.tocsc())
    report.update(input_sha256=identity,panel_sha256=panel_hash,active_physical_coefficients=len(patch.active_nodes),
        covered_auxiliary_coefficients=len(patch.covered_nodes),source_triangles=len(patch.faces),
        physical_mass_diagonal_min=float(patch.mass.diagonal().min()),
        physical_mass_diagonal_max=float(patch.mass.diagonal().max()),
        preparation_seconds=perf_counter()-started,quadrature_points=query_count,orders=[])
    saved=dict(x=x,y=y,polygon=polygon,active_nodes=patch.active_nodes,covered_nodes=patch.covered_nodes,
        mass_data=patch.mass.data,mass_indices=patch.mass.indices,mass_indptr=patch.mass.indptr,
        pmap=pmap,nmap=nmap,even=even,neumann=neumann,far_body=far_body,far_normal=far_normal)
    for order in args.orders:
        print(json.dumps(dict(order=order,start=True)),flush=True);started=perf_counter()
        cross=np.zeros((nn,nb));psi_mom=np.zeros((nn,4));far_mom=full_far.copy()
        for points,nodes,weights in cut_quadrature(patch,order):
            scatter_moments(cross,nodes,weights,source_trace(points))
            scatter_moments(psi_mom,nodes,weights,phases(reference(points)))
        for points,nodes,weights in cut_quadrature(patch,order,hidden=True):
            scatter_moments(far_mom,nodes,-weights,fourier_values(far,points[:,:2],origin,length))
        cross=cross[patch.active_nodes];psi_mom=psi_mom[patch.active_nodes];far_mom=far_mom[patch.active_nodes]
        inverse_cross=mass_factor.solve(cross)
        schur=neumann-nmap@inverse_cross;condition=float(np.linalg.cond(schur))
        if not np.isfinite(condition): raise RuntimeError('singular source Schur complement; no diagonal repair')
        factor=lu_factor(schur)

        def solve(physical_psi_mom,far_trace_mom,rigid,far_n):
            direct=mass_factor.solve(physical_psi_mom-far_trace_mom)
            sigma=lu_solve(factor,rigid-far_n-nmap@direct)
            return direct-inverse_cross@sigma,sigma

        coefficients,sigma=solve(psi_mom,far_mom,q,far_normal)
        eta_mom=phases(1j*omega/GRAVITY*amplitude(psi_mom))
        derivative,sigma_dot=solve(-GRAVITY*eta_mom,rate(far_mom),rate(q),rate(far_normal))
        body_phi=far_body+pmap@coefficients+even@sigma
        body_phi_dot=rate(far_body)+pmap@derivative+even@sigma_dot
        force=-WATER_DENSITY*weighted_modes.T@amplitude(body_phi_dot)
        no_cross=mass_factor.solve(-GRAVITY*eta_mom-rate(far_mom))
        sigma_dot_no_cross=lu_solve(nfactor,rate(q)-rate(far_normal)-nmap@no_cross)
        force_no_cross=-WATER_DENSITY*weighted_modes.T@amplitude(rate(far_body)+pmap@no_cross+even@sigma_dot_no_cross)
        force_no_source=-WATER_DENSITY*weighted_modes.T@amplitude(rate(far_body)+pmap@derivative)
        state_residual=patch.mass@coefficients+cross@sigma+far_mom-psi_mom
        dynamic_residual=patch.mass@derivative+cross@sigma_dot+rate(far_mom)+GRAVITY*eta_mom
        normal_residual=nmap@derivative+neumann@sigma_dot+rate(far_normal)-rate(q)
        body_error=amplitude(body_phi)-reference_body
        entry=dict(order=order,construction_and_reference_seconds=perf_counter()-started,
            source_schur_condition_number=condition,
            body_potential_relative_l2=np.sqrt(np.sum(panels.areas[:,None]*abs(body_error)**2,axis=0)
                /np.sum(panels.areas[:,None]*abs(reference_body)**2,axis=0)).tolist(),
            six_dof_force_absolute_error=abs(force-reference_force).tolist(),
            principal_force_relative_error=(abs(force-reference_force)[[2,3],[0,1]]/abs(reference_force)[[2,3],[0,1]]).tolist(),
            omission_cross_principal_force_relative_error=(abs(force_no_cross-reference_force)[[2,3],[0,1]]/abs(reference_force)[[2,3],[0,1]]).tolist(),
            omission_source_pressure_principal_force_relative_error=(abs(force_no_source-reference_force)[[2,3],[0,1]]/abs(reference_force)[[2,3],[0,1]]).tolist(),
            state_weak_trace_residual_max=float(abs(state_residual).max()),
            dynamic_weak_trace_residual_max=float(abs(dynamic_residual).max()),
            derivative_normal_residual_max=float(abs(normal_residual).max()),
            derivative_vs_rotated_state_max=float(abs(derivative-rate(coefficients)).max()),
            sigma_derivative_vs_rotated_state_max=float(abs(sigma_dot-rate(sigma)).max()))
        # Independent physical field samples on the original straight cut hat.
        base_case=contact['cases'][0]
        geometry=cut_hat(base_case['centre'],*contact_test_spacing(contact),polygon)
        ep,ew,ev,eg,_=triangle_rule(geometry['exterior'],8)
        full=patch.expand(coefficients).reshape(len(y),len(x),4)
        field_owner=PlaneP1Field(x,y,full)
        owners=field_owner.point_owners(ep)
        local_values=full.reshape(-1,4)[field_owner.faces[owners,0]]+np.einsum('mi,mic->mc',
            ep-field_owner.panels.triangles[owners,0,:2],field_owner.gradient[owners])
        points=np.column_stack((ep,np.zeros(len(ep))))
        reconstructed=amplitude(local_values+fourier_values(far,ep,origin,length)+source_trace(points)@sigma)
        reference_trace=reference(points)
        entry['straight_physical_trace_weighted_rms_error']=np.sqrt((ew*ev)@abs(reconstructed-reference_trace)**2
            /geometry['physical_hat_mass']).tolist()
        a=field_owner.vector_potential(ep,owners)
        local_flux=2*np.einsum('m,mi,mic->c',ew,eg,a)
        line=geometry['waterline'];nodes,weights=leggauss(8);t=(nodes+1)/2;weights=weights/2
        lp=(line['segments'][:,:1]+t[None,:,None]*(line['segments'][:,1:]-line['segments'][:,:1])).reshape(-1,2)
        lv=((1-t)*line['values'][:,:1]+t*line['values'][:,1:]).ravel()
        lw=(np.linalg.norm(line['segments'][:,1]-line['segments'][:,0],axis=1)[:,None]*weights).ravel()*lv
        a=field_owner.vector_potential(lp,field_owner.point_owners(lp))
        local_flux+=2*np.einsum('m,mi,mic->c',lw,np.repeat(line['hole_outward_normals'],8,axis=0),a)
        candidate_flux=amplitude((local_flux+(ew*ev)@fourier_values(far*np.hypot(kx,ky)[:,:,None],ep,origin,length))[None])[0]
        actual_flux=(ew*ev)@reference(points,velocity=True)
        robin_flux=omega**2/GRAVITY*((ew*ev)@reference_trace)
        entry['straight_actual_velocity_moment_absolute_error']=abs(candidate_flux-actual_flux).tolist()
        entry['reference_robin_moment_defect']=abs(actual_flux-robin_flux).tolist()
        entry['straight_candidate_vs_robin_moment_error']=abs(candidate_flux-robin_flux).tolist()
        report['orders'].append(entry)
        saved[f'q{order}_cross']=cross;saved[f'q{order}_psi_mom']=psi_mom
        saved[f'q{order}_coefficients']=coefficients;saved[f'q{order}_sigma']=sigma
        saved[f'q{order}_derivative']=derivative;saved[f'q{order}_sigma_dot']=sigma_dot
        saved[f'q{order}_force']=force;saved[f'q{order}_body_phi']=body_phi
        np.savez(args.output.with_suffix('.npz'),**saved)
        report['arrays_sha256']=sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
        print(json.dumps(entry),flush=True)
    report.update(complete=True,qualified=False,scope=__doc__,peak_process_rss_bytes=peak_rss_bytes(),
        limits=['Both physical pressure terms are tested, but eta is externally supplied; kinematic evolution does not generate the state',
            'Far-field chi and its derivative are prescribed reference data; no open-water boundary qualification',
            'All imposed harmonic rates are rotations of state: phase agreement is algebraic; independent fields/forces supply physical evidence',
            'Float cut geometry and algebraic hidden-chi extension are verification references, not moving topology certificates',
            'No nonlinear moving geometry, complete energy proof, full6DOF trajectory,60second accuracy or simulation throughput qualification'])


def contact_test_spacing(contact):
    with np.load(Path(contact['parents']['snapshot']['path']).with_suffix('.npz')) as z:
        return float(z['patch_x'][1]-z['patch_x'][0]),float(z['patch_y'][1]-z['patch_y'][0])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('contact','output'): parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--orders',type=int,nargs='+',default=[4,8])
    parser.add_argument('--coarsening',type=int,choices=[0,1,2],default=1)
    parser.add_argument('--maximum-quadrature-points',type=int,default=2_000_000)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or not args.orders or min(args.orders)<2 or max(args.orders)>8:
        parser.error('fresh archive and integration orders2..8 required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),
        *Path(__file__).parent.glob('*.py'),*Path(__file__).parent.glob('*.cu'),*Path(__file__).parent.glob('*.cuh')]:
        relative=source.resolve().relative_to(Path.cwd())
        target=archive/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    shutil.copyfile(args.contact,archive/'contact.json')
    report=dict(complete=False,producer_sha256=hashes,parent_contact_sha256=sha256(args.contact.read_bytes()).hexdigest())
    started=perf_counter()
    try: run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
