"""Exact-current mixed field derivative on one guarded, fixed nodal chart.

The only time stencil supplied here is the prescribed incident boundary load.
No offset fluid operator or field is used. Body acceleration is deliberately
zero; the current acceleration-potential columns supply added mass separately.
"""
from time import perf_counter
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from tests.reference_boat import WATER_DENSITY
from .cuda_directional import energy_direction_action
from .boat_world_field import boundary_data


def incident_load_rate(base,plus,minus,values,twist,time,step):
    """Cheap h-stencil of the SAME prescribed incident Neumann quadrature."""
    if not np.isfinite(step) or step<=0:
        raise ValueError('positive finite incident derivative step required')
    basis=('faces','body_faces','free_faces','free_nodes','body_source_group',
           'free_source_parent','free_parent_indices','node_chart_key')
    if 'node_chart_key' not in base:
        raise ValueError('explicit current chart ownership required')
    for moved in (plus,minus):
        if any(not np.array_equal(base[name],moved[name]) for name in basis):
            raise ValueError('incident derivative must preserve current chart/basis')
    loads=[]
    for moved,offset in ((plus,step),(minus,-step)):
        boundary=boundary_data(moved,values,twist,time=time+offset)
        loads.append(boundary['right']-boundary['body_load']@twist)
    return (loads[0]-loads[1])/(2*step)


def field_momentum_rate(arrays,assembly,field,values,twist,time,eta_rate,psi_rate,
                        point_jacobian,map_jacobian,body_load_jacobian,incident_rate,
                        *,body_pose_load,density=WATER_DENSITY):
    """Return ALE phi_dot and dF_V/dt for the current zero-acceleration path.

    This contained prototype factors E_UU locally, because current_field does
    not expose its existing factor. It never mistakes the assembly's V factor
    for E_UU. Production retention can eliminate that duplicate factor later.
    """
    started=perf_counter()
    points,faces=arrays['points'],arrays['faces']
    nv,np_=len(points),len(arrays['parent_nodes']);nd=6+np_
    twist=np.asarray(twist,dtype=float)
    inputs=((twist,(6,)),(eta_rate,(np_,)),(psi_rate,(np_,)),
            (point_jacobian,(nv,3,nd)),(map_jacobian,(len(arrays['free_nodes']),np_,nd)),
            (body_load_jacobian,(nv,6,nd)),(incident_rate,(nv,)),(body_pose_load,(6,)))
    if any(np.shape(value)!=shape or not np.isfinite(value).all() for value,shape in inputs):
        raise ValueError('finite current chart tangents, rates and loads required')
    if not np.isfinite(density) or density<=0:
        raise ValueError('positive finite density required')
    if (not assembly.uniform_quadrature or not np.array_equal(faces,assembly.faces)
            or not np.array_equal(points[faces],assembly.triangles)
            or np.any(assembly.normal_sign!=1)):
        raise ValueError('field derivative requires exact current uniform-quadrature assembly')
    boundary=boundary_data(arrays,values,twist,time=time)
    for key,owned in (('parent_map',boundary['parent_map']),('body_load',boundary['body_load']),('right',boundary['right'])):
        if not np.array_equal(field[key],owned):
            raise ValueError('field does not own the current physical boundary data')
    phi=np.asarray(field['potential'],dtype=float)
    fixed=arrays['free_nodes'];unknown=np.setdiff1d(np.arange(nv),fixed)
    if (phi.shape!=(nv,) or not np.isfinite(phi).all()
            or not np.array_equal(phi[fixed],field['parent_map']@arrays['parent_psi'])):
        raise ValueError('field does not own the current prescribed Dirichlet trace')
    residual=assembly.energy@phi-field['right']
    relative=np.linalg.norm(residual[unknown])/max(np.linalg.norm(field['right'][unknown]),1e-30)
    if relative>1e-9:
        raise ValueError('current mixed field residual exceeds ownership check')
    direction=np.r_[twist,eta_rate]
    velocity=np.einsum('nid,d->ni',point_jacobian,direction)
    map_rate=np.einsum('npd,d->np',map_jacobian,direction)
    b_rate=np.einsum('nkd,d->nk',body_load_jacobian,direction)
    prescribed_rate=map_rate@arrays['parent_psi']+field['parent_map']@psi_rate
    phi_rate=np.zeros(nv);phi_rate[fixed]=prescribed_rate
    before=perf_counter()
    e_rate_phi,diagnostics=energy_direction_action(points,faces,velocity,phi,assembly=assembly)
    action_seconds=perf_counter()-before
    # Known-load path: twist_dot=0, not a total accelerated-body derivative.
    right_rate=b_rate@twist+incident_rate
    rhs=right_rate-e_rate_phi-assembly.energy@phi_rate
    before=perf_counter()
    factor=cho_factor(assembly.energy[np.ix_(unknown,unknown)])
    phi_rate[unknown]=cho_solve(factor,rhs[unknown])
    solve_seconds=perf_counter()-before
    differentiated_residual=assembly.energy@phi_rate+e_rate_phi-right_rate
    momentum=-density*(field['body_load'].T@phi)
    transport=-density*(b_rate.T@phi+field['body_load'].T@phi_rate)
    angular=np.r_[np.zeros(3),-np.cross(twist[3:],momentum[3:])]
    known_load=np.asarray(body_pose_load)+transport+angular
    if not np.isfinite(np.r_[phi_rate,known_load]).all():
        raise FloatingPointError('nonfinite exact-current field/momentum derivative')
    return dict(potential_rate_ale=phi_rate,mesh_velocity=velocity,body_load_rate=b_rate,
        momentum_rate=transport,angular_transport=angular,known_load=known_load,
        energy_rate_action=e_rate_phi,right_rate=right_rate,
        diagnostics=dict(diagnostics,current_field_relative_residual=float(relative),
            differentiated_field_relative_residual=float(np.linalg.norm(differentiated_residual[unknown])/max(np.linalg.norm(right_rate[unknown]),1e-30)),
            directional_action_seconds=action_seconds,mixed_factor_and_solve_seconds=solve_seconds,
            total_seconds=perf_counter()-started))
