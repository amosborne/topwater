"""One actual six-mode Ainf source producer, with bounded dense lifetimes."""
from time import perf_counter
from hashlib import sha256
import numpy as np
from .paired_body import scaled_check,capytaine_radius_check
from .frozen_source_quadrature import pressure_matrix


def canonical_matrix_hash(value):
    h=sha256()
    for row in value:h.update(np.ascontiguousarray(row).view(np.uint8))
    return h.hexdigest()


def validate_infinite_basis(phi,sigma,boundary,panels,completed,tolerance=1e-10):
    if completed!=6:raise ValueError('all six freshly solved infinite modes required')
    for array in (phi,sigma,boundary):
        if array.shape!=(panels,6) or not np.isfinite(array).all():raise ValueError('missing/nonfinite six-mode infinite basis')
        if np.any(abs(array.imag)>tolerance*(1+abs(array.real))):raise ValueError('unexpected complex infinite-frequency basis')


def export_basis(cpt,green,triangles,center,mass,original_A,scale,budget,progress,on_solved):
    from tests.capytaine_reference import body_from_triangles,MODES
    begin=perf_counter();n=len(triangles)
    body=body_from_triangles(triangles,center,mass)
    if body.mesh_including_lid is not body.mesh:raise ValueError('plain no-lid dense allocation path required')
    record=dict(radius_error_m=capytaine_radius_check(body,triangles),modes=list(MODES),cache_reuse=[])
    solver=cpt.BEMSolver(green_function=green)
    if solver.engine.exportable_settings['linear_solver']!='lu_decomposition':raise ValueError('declared non-overwriting LU required')
    phi=np.zeros((n,6),complex);sigma=phi.copy();q=phi.copy();A=np.zeros((6,6))
    arrays=dict(source_triangles=triangles,native_potential=phi,normalized_sources=sigma,native_boundary=q,native_Ainf=A,native_points=body.mesh.faces_centers)
    for j,mode in enumerate(MODES):
        budget();pb=cpt.RadiationProblem(body=body,radiating_dof=mode,omega=np.inf,rho=1025.,g=9.81,water_depth=np.inf,forward_speed=0.)
        request=dict(free_surface=pb.free_surface,water_depth=pb.water_depth,wavenumber=pb.wavenumber,
            adjoint_double_layer=True,diagonal_term_in_double_layer=True)
        expected=(body.mesh_including_lid,body.mesh_including_lid,request)
        if j==0:
            original_S,original_K=solver.engine.build_matrices(body.mesh_including_lid,body.mesh_including_lid,**request);budget()
            if (not isinstance(original_S,np.ndarray) or not isinstance(original_K,np.ndarray)
                    or original_S.shape!=(n,n) or original_K.shape!=(n,n)
                    or original_S.dtype!=np.dtype(complex) or original_K.dtype!=np.dtype(complex)):raise ValueError('unmodeled matrix layout')
            S=original_S.copy();K=original_K.copy();budget()
            record['pristine_operators']={name:dict(shape=list(value.shape),dtype=str(value.dtype),sha256=canonical_matrix_hash(value)) for name,value in [('S',S),('K',K)]}
            if solver.engine.last_computed_inputs!=expected:raise ValueError('warmup cache key differs from solve request')
        budget();result=solver.solve(pb,method='indirect',keep_details=True);on_solved(j+1);budget()
        if solver.engine.last_computed_inputs!=expected or solver.engine.last_computed_matrices[0] is not original_S:raise ValueError('unexpected matrix rebuild invalidates allocation model')
        record['cache_reuse'].append(True)
        factor=-1j*pb.encounter_omega
        phi[:,j]=np.asarray(result.potential/factor);sigma[:,j]=np.asarray(result.sources/factor);q[:,j]=np.asarray(pb.boundary_condition/factor)
        A[:,j]=[result.added_mass[m] for m in MODES]
        if any(result.radiation_damping[m]!=0 for m in MODES):raise ValueError('nonzero infinite-frequency damping')
        progress(arrays,j+1);budget()
    validate_infinite_basis(phi,sigma,q,n,6)
    points=body.mesh.faces_centers;normal=-body.mesh.faces_normals;N=np.c_[normal,np.cross(points-center,normal)]
    produced_phi=S@sigma;produced_q=K@sigma;native_load=pressure_matrix(points,normal,body.mesh.faces_areas,phi,center,1025.)
    record['source_checks']=dict(rigid_normal=scaled_check(q,-N,1.,1e-10),potential=scaled_check(produced_phi,phi,1.,1e-10),
        normal=scaled_check(produced_q,q,1.,1e-10),native_load=scaled_check(native_load,A,scale,1e-10),
        original_native_Ainf=scaled_check(A,original_A,scale,1e-10))
    arrays.update(operator_potential=produced_phi,operator_normal=produced_q,potential_residual=produced_phi-phi,normal_residual=produced_q-q,native_integrated_A_complex=native_load)
    record['source_export_seconds']=perf_counter()-begin
    del S,K,original_S,original_K,result,pb,solver,body
    return arrays,record
