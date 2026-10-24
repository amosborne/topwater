"""Existing closed-boundary P1/P0 energy on an exact odd-mirrored wet hull."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from scipy.linalg import cholesky, cho_solve, eigvalsh, solve_triangular
from topwater import HullMesh
from topwater.constants import WATER_DENSITY
from topwater.rankine import TrianglePanels
from topwater.surface import body_load
from tests.laplace_reference import triangle_samples
from tests.benchmark_rankine import peak_rss_bytes
from .regularized_snapshot import owned
from .operator import assemble_energy
from .boat_spatial_work import subdivide


def checkpoint(args, report, arrays):
    np.savez(args.output.with_suffix('.npz'), **arrays)
    report['arrays_sha256'] = sha256(args.output.with_suffix('.npz').read_bytes()).hexdigest()
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False))


def ratio(a, b):
    return float(a/b) if b else (0. if a == 0 else None)


def mirrored_mesh(triangles):
    """Reflect/reverse; only exactly equal coordinates share vertices."""
    if np.any(triangles[:, :, 2] > 0) or np.any(triangles[:, :, 2].mean(axis=1) >= 0):
        raise ValueError('physical input facets must be wetted strictly below the plane in their interiors')
    full = np.concatenate((triangles, triangles[:, ::-1]*[1., 1., -1.]))
    vertices, inverse = np.unique(full.reshape(-1, 3), axis=0, return_inverse=True)
    mesh = HullMesh(vertices, inverse.reshape(-1, 3))
    lower = np.flatnonzero(vertices[:, 2] < 0)
    lookup = {tuple(point): i for i, point in enumerate(vertices)}
    upper = np.array([lookup[tuple(vertices[i]*[1., 1., -1.])] for i in lower])
    odd = np.zeros((len(vertices), len(lower)))
    odd[lower, np.arange(len(lower))] = 1.
    odd[upper, np.arange(len(lower))] = -1.
    if not len(lower) or not np.any(vertices[:, 2] == 0):
        raise ValueError('this surface-piercing verification requires off-plane and exact seam vertices')
    np.testing.assert_array_equal(mesh.vertices[mesh.faces[:len(triangles)]], triangles)
    return mesh, odd


def refine_once(triangles):
    """Existing midpoint owner and reviewed arithmetic geometry preflight.

    The same plane/area/moment checks were reviewed for
    boat_odd_body_refinement; skinny normalized-normal differences are not a
    geometry acceptance condition. No vertex changes beyond midpoint arithmetic.
    """
    old=TrianglePanels(triangles,-1)
    fine=subdivide(triangles)
    panels=TrianglePanels(fine,-1)
    parent=np.tile(np.arange(len(old)),4)
    np.testing.assert_allclose(panels.areas,np.tile(old.areas/4,4),rtol=1e-12,atol=1e-14)
    original_offset=triangles-triangles[:,:1]
    original_defect=np.max(abs(np.einsum('tvi,ti->tv',original_offset,old.geometric_normals)),axis=1)
    offset=fine-triangles[parent,:1]
    defect=abs(np.einsum('tvi,ti->tv',offset,old.geometric_normals[parent]))
    allowance=original_defect[parent,None]+64*np.finfo(float).eps*np.sum(
        abs(fine)+abs(triangles[parent,:1]),axis=2)
    if np.any(defect>allowance):
        raise ValueError('midpoint children leave original parent planes beyond arithmetic allowance')
    if np.any(np.sum(panels.normals*old.normals[parent],axis=1)<=0):
        raise ValueError('midpoint subdivision reverses physical winding')
    first=(panels.areas[:,None]*panels.centroids).reshape(4,len(old),3).sum(axis=0)
    original_first=old.areas[:,None]*old.centroids
    first_allowance=(abs(panels.areas.reshape(4,len(old))-old.areas[None,:]/4)[:,:,None]
        *abs(panels.centroids.reshape(4,len(old),3))).sum(axis=0)
    first_allowance+=64*np.finfo(float).eps*old.areas[:,None]*np.max(abs(triangles),axis=1)
    if np.any(abs(first-original_first)>first_allowance):
        raise ValueError('midpoint subdivision changes original first moments beyond arithmetic allowance')
    return fine,dict(max_parent_plane_defect_m=float(defect.max()),
        max_parent_plane_arithmetic_allowance_m=float(allowance.max()),
        area_difference=float(panels.areas.sum()-old.areas.sum()),
        first_moment_difference=(first-original_first).sum(axis=0).tolist(),
        max_normal_difference_report_only=float(abs(panels.normals-old.normals[parent]).max()),
        scope='Reviewed floating arithmetic plane/area/first-moment preflight, not a certified interval')


def run(args, report):
    family, reference = owned(args.family), owned(args.reference)
    if family['parents']['reference'] != sha256(args.reference.read_bytes()).hexdigest():
        raise ValueError('cached physical six-port reference ownership differs')
    for name in ('water','phase','mode'):
        if family['parents'][name] != reference['parents'][name]:
            raise ValueError('cached physical parents differ')
    for name in ('src/topwater/rankine.py','src/topwater/inputs.py','src/topwater/surface.py'):
        relative=Path(name)
        if sha256(relative.read_bytes()).hexdigest()!=family['producer_sha256'].get(str(relative)):
            raise ValueError('physical geometry/load owner differs: '+name)
    with np.load(args.family.with_suffix('.npz')) as z:
        triangles, center = z['triangles'].copy(), z['world_cg'].copy()
    with np.load(args.reference.with_suffix('.npz')) as z:
        cached={name:z[name].copy() for name in ('ritz_added_mass','direct_added_mass','capytaine_added_mass','dry_mass')}
        np.testing.assert_array_equal(triangles,z['triangles'])
        np.testing.assert_array_equal(center,z['world_cg'])
    parent_triangle_sha=sha256(triangles.tobytes()).hexdigest()
    if parent_triangle_sha!=reference['panel_sha256']:
        raise ValueError('original physical parent triangle SHA differs before numerical refinement')
    refinement_geometry=None
    if args.refinement:
        triangles,refinement_geometry=refine_once(triangles)
    orders=(4,8,16) if args.refinement==0 else (4,8)
    mesh, odd = mirrored_mesh(triangles)
    panels = TrianglePanels.from_mesh(mesh, normal_sign=-1)
    nf,nv,nb = len(panels),len(mesh.vertices),len(triangles)
    lower_load=body_load(mesh,panels,np.arange(nb),center)
    B=odd.T@lower_load
    pair_bound=nf*nf*sum(q*q for q in orders)
    memory=8*(8*nf*nf+8*nf*nv+12*nv*nv+nf*max(orders)**2*4)+32*1024**2
    if pair_bound>args.maximum_pairs or memory>args.maximum_bytes:
        raise RuntimeError('mirrored energy source/workspace preflight exceeds budget')
    factor_dry=cholesky(cached['dry_mass'],lower=True)

    def whiten(matrix):
        left=solve_triangular(factor_dry,matrix,lower=True)
        return solve_triangular(factor_dry,left.T,lower=True).T

    qpoints,qweights=triangle_samples(triangles,4)
    bary,_=triangle_samples(np.eye(3)[None],4)
    bary=bary[0]
    normal=panels.normals[:nb,None,:]
    T=np.concatenate((np.broadcast_to(normal,qpoints.shape),np.cross(qpoints-center,normal)),axis=2).reshape(-1,6)
    weights=(panels.areas[:nb,None]*qweights).ravel()
    arrays=dict(vertices=mesh.vertices,faces=mesh.faces,odd_lift=odd,lower_load=lower_load,
        odd_lower_load=B,triangles=triangles,world_cg=center,**cached)
    report.update(complete=False,qualified=False,body_panels=nb,mirrored_panels=nf,
        refinement_level=args.refinement,parent_triangle_sha256=parent_triangle_sha,
        current_triangle_sha256=sha256(triangles.tobytes()).hexdigest(),
        refinement_geometry=refinement_geometry,predeclared_orders=list(orders),
        closed_nodes=nv,odd_nodes=odd.shape[1],exact_seam_nodes=int((mesh.vertices[:,2]==0).sum()),
        source_pair_bound=pair_bound,array_workspace_bound_bytes=memory,density_kg_m3=WATER_DENSITY,
        memory_scope='Array/workspace preflight excludes Python/Numba/library baseline; not a process-RSS cap',
        geometry='Exact reflection and reversed upper winding; exact coordinate equality only; HullMesh closed/winding/volume validation',
        energy='Half of existing CLOSED exterior energy restricted to odd nodal trace; seam trace exactly zero',
        cached_reference_raw_asymmetry_max={name:float(abs(value-value.T).max()) for name,value in cached.items()
            if name.endswith('added_mass')},
        body_orders=[],archived_refinement=reference['archived_body_refinement'])
    checkpoint(args,report,arrays)
    previous=None
    for order in orders:
        started=perf_counter();row=dict(order=order,accepted=False,source_pair_reservation=nf*nf*order*order)
        report['body_orders'].append(row)
        try:
            closed,diagnostics=assemble_energy(mesh.faces,panels,order=order)
            half=.5*(odd.T@closed@odd)
            arrays[f'q{order}_closed_energy']=closed
            arrays[f'q{order}_half_energy']=half
            row['assembly_diagnostics']=diagnostics
            error=float(abs(half-half.T).max())
            allowance=128*np.finfo(float).eps*(nf+nv)*float(abs(half).max())
            row['half_energy_transpose_gate']=dict(maximum_absolute=error,arithmetic_allowance=allowance)
            if error>allowance:
                raise ValueError('half energy exceeds transpose arithmetic allowance; no repair')
            spectrum=eigvalsh(half)
            row['half_energy_eigenvalues']=spectrum.tolist()
            if spectrum[0]<=0:
                raise ValueError('odd half energy not positive definite; no jitter/clipping')
            factor=cholesky(half,lower=True)
            response=cho_solve((factor,True),B)
            A=WATER_DENSITY*B.T@response
            full=odd@response
            nodal=full[mesh.faces[:nb]]
            potential=np.einsum('qj,fjc->fqc',bary,nodal).reshape(-1,6)
            physical_port=WATER_DENSITY*T.T@(weights[:,None]*potential)
            port_error=float(abs(A-physical_port).max())
            scale=float((WATER_DENSITY*abs(T).T@(weights[:,None]*abs(potential))).max())
            port_allowance=128*np.finfo(float).eps*(len(weights)+nv)*scale
            if port_error>port_allowance:
                raise ValueError('independent polynomial P1 load-port quadrature parity failed')
            row.update(accepted=True,diagonal=A.diagonal().tolist(),
                raw_added_mass_asymmetry_max=float(abs(A-A.T).max()),
                half_energy_condition=float(spectrum[-1]/spectrum[0]),
                independent_q4_load_port_parity_max=port_error,
                independent_q4_load_port_arithmetic_allowance=port_allowance)
            for name in ('ritz_added_mass','direct_added_mass','capytaine_added_mass'):
                difference=A-cached[name];white=whiten(difference)
                row[name+'_comparison']=dict(diagonal_relative_difference=[ratio(a,b) for a,b in
                    zip(difference.diagonal(),cached[name].diagonal(),strict=True)],
                    dry_mass_whitened_difference=white.tolist(),
                    dry_mass_whitened_relative_operator_norm=ratio(np.linalg.norm(white,2),np.linalg.norm(whiten(cached[name]),2)))
            row['archived_refinement_diagonal_comparisons']=[]
            for item in reference['archived_body_refinement']:
                row['archived_refinement_diagonal_comparisons'].append(dict(
                    level=item['level'],panels=item['panels'],
                    direct_relative_difference=[ratio(a-b,b) for a,b in zip(A.diagonal(),item['direct_diagonal'],strict=True)],
                    source_relative_difference=[ratio(a-b,b) for a,b in zip(A.diagonal(),item['source_diagonal'],strict=True)]))
            if previous is not None:
                row['previous_order_added_mass_difference_max']=float(abs(A-previous).max())
                row['previous_order_diagonal_relative_change']=[ratio(a-b,b) for a,b in zip(A.diagonal(),previous.diagonal(),strict=True)]
                row['previous_order_dry_mass_whitened_change_max']=float(abs(whiten(A-previous)).max())
            previous=A
            arrays.update({f'q{order}_added_mass':A,f'q{order}_odd_response':response,
                f'q{order}_physical_load_port':physical_port})
        except (ValueError,np.linalg.LinAlgError) as error:
            row['rejected']=dict(type=type(error).__name__,message=str(error))
        row['seconds']=perf_counter()-started
        checkpoint(args,report,arrays)
        print(json.dumps(row),flush=True)
    report.update(complete=True,process_peak_rss=peak_rss_bytes(),limits=[
        'Existing assembler averages raw V quadrature asymmetry as part of its symmetric Galerkin definition; diagnostics retained',
        'No new kernel, fitted-force correction, energy eigenvalue clipping, or jitter',
        'Physical q4 port checks polynomial lower-P1 load reconstruction, not the half-energy factor or continuum added-mass accuracy',
        'The half-energy factor follows odd-reflection energy symmetry and is tested physically through external added-mass comparisons',
        'Independent references retain their own discretization and raw asymmetry',
        'Refinement changes both P1 trace and P0 auxiliary flux spaces; strict monotonic added mass is not asserted',
        'Whole-plane zero-potential instantaneous/infinite-frequency limit only; no finite-frequency, moving, nonlinear6DOF or trajectory qualification'])
    checkpoint(args,report,arrays)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('family','reference','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--maximum-pairs',type=int,default=10_000_000)
    parser.add_argument('--maximum-bytes',type=int,default=128*1024**2)
    parser.add_argument('--refinement',type=int,choices=(0,1),default=0)
    args=parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.source').exists() or min(args.maximum_pairs,args.maximum_bytes)<=0:
        parser.error('fresh source/output and positive bounds required')
    archive=args.output.with_suffix('.source');archive.mkdir(parents=True)
    hashes={}
    for source in [*Path('src/topwater').glob('*.py'),*Path('tests').glob('*.py'),*Path(__file__).parent.glob('*.py')]:
        relative=source.resolve().relative_to(Path.cwd());target=archive/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        hashes[str(relative)]=sha256(source.read_bytes()).hexdigest()
    report=dict(complete=False,producer_sha256=hashes,parents={name:sha256(getattr(args,name).read_bytes()).hexdigest()
        for name in ('family','reference')})
    started=perf_counter()
    try:
        run(args,report)
    except BaseException as error:
        report['failure']=dict(type=type(error).__name__,message=str(error));raise
    finally:
        report['total_seconds']=perf_counter()-started
        args.output.write_text(json.dumps(report,indent=2,allow_nan=False))
