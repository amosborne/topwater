"""Current-domain rigid-transport acceleration potential, prescribed acceleration.

With phi=H+R and V=U+Omega cross r, Z=phi_t+V dot grad(R) is harmonic.
Its body data require only Fourier derivatives of H and current source velocity.
The full moving-operator finite-difference oracle remains a separate harness.
"""
from time import perf_counter
import numpy as np

from topwater.constants import GRAVITY, WATER_DENSITY
from topwater.frames import rotation_matrix
from .boat_graph_fields import source_velocity_body, source_velocity_bulk
from .cuda_graph import harmonic_velocity, cp


def normal_transport_lift(context, spectrum, state):
    """Velocity of H and dn(V dot grad(H)), using existing continuation maps."""
    panels = context.panels
    omega = rotation_matrix(state.orientation)@state.angular_velocity
    rigid = state.linear_velocity+np.cross(omega, panels.centroids-state.position)
    size = context.graph.padded_shape[0]
    lift_velocity = cp.asnumpy(harmonic_velocity(context.potential, spectrum, size, context.length))
    kx, ky = cp.meshgrid(2*np.pi*cp.fft.fftfreq(size, d=context.length/size),
                         2*np.pi*cp.fft.fftfreq(size, d=context.length/size))
    derivatives = (1j*kx.ravel(), 1j*ky.ravel(), cp.hypot(kx, ky).ravel())
    normal_hessian_velocity = sum(rigid[:, axis]*cp.asnumpy(
        (context.normal@(spectrum*d)).real) for axis, d in enumerate(derivatives))
    normal_transport = np.sum(np.cross(omega, panels.normals)*lift_velocity, axis=1)+normal_hessian_velocity
    return lift_velocity, normal_transport


def prescribed_pressure(context, density, spectrum, state, acceleration, extension, chi_exterior):
    """Recover phi_t and pressure with the same CURRENT body/graph operator.

    acceleration is [world CG acceleration, world angular acceleration].
    chi_exterior is Eulerian phi_t data at physical exterior graph nodes:
    psi_t-u_z*eta_t for a forced graph, or -g*eta-.5*|u|^2 for physical p=0.
    No sponge/state forcing is silently interpreted as physical pressure.
    """
    started = perf_counter()
    acceleration = np.asarray(acceleration, dtype=float)
    chi_exterior = np.asarray(chi_exterior, dtype=float)
    if (acceleration.shape != (6,) or not np.isfinite(acceleration).all()
            or chi_exterior.shape != (len(extension.exterior),) or not np.isfinite(chi_exterior).all()):
        raise ValueError('finite world six-acceleration and physical exterior phi_t data required')
    panels = context.panels
    omega = rotation_matrix(state.orientation)@state.angular_velocity
    relative = panels.centroids-state.position
    rigid = state.linear_velocity+np.cross(omega, relative)
    lift_velocity, normal_transport = normal_transport_lift(context, spectrum, state)
    source_velocity = source_velocity_body(context, density)
    total_velocity = lift_velocity+source_velocity
    rigid_acceleration = (acceleration[:3]+np.cross(acceleration[3:], relative)
                          -np.cross(omega, state.linear_velocity))
    normal_data = np.sum(rigid_acceleration*panels.normals, axis=1)-normal_transport
    source_at_surface = source_velocity_bulk(context, density, context.surface_points)
    rigid_at_surface = state.linear_velocity+np.cross(omega, context.surface_points-state.position)
    known_transport = np.sum(rigid_at_surface*source_at_surface, axis=1)
    source_at_samples = source_velocity_bulk(context, density, context.source_points)
    rigid_at_samples = state.linear_velocity+np.cross(omega, context.source_points-state.position)
    transport_samples = np.sum(rigid_at_samples*source_at_samples, axis=1)
    trace = extension.extend(chi_exterior)+context.project_source_trace(transport_samples)
    grid = context.shape[0]
    kx0, ky0 = cp.meshgrid(2*np.pi*cp.fft.fftfreq(grid, d=context.length/grid),
                           2*np.pi*cp.fft.fftfreq(grid, d=context.length/grid))
    r_trace = (context.trace@cp.asarray(density)).reshape(context.shape)
    r_spectrum = cp.fft.fft2(r_trace)
    eta_spectrum = cp.fft.fft2(cp.asarray(context.surface_points[:, 2].reshape(context.shape)))
    spectral_r_gradient = np.column_stack([cp.asnumpy(cp.fft.ifft2(1j*k*r_spectrum).real).ravel()
                                           for k in [kx0, ky0]])
    eta_gradient = np.column_stack([cp.asnumpy(cp.fft.ifft2(1j*k*eta_spectrum).real).ravel()
                                    for k in [kx0, ky0]])
    analytic_r_gradient = source_at_surface[:, :2]+source_at_surface[:, 2, None]*eta_gradient
    derivative_alias = analytic_r_gradient-spectral_r_gradient
    transport_alias = np.sum(rigid_at_surface[:, :2]*derivative_alias, axis=1)
    field_seconds = perf_counter()-started
    solution = context.solve(trace, normal_data)
    phi_t = solution.potential-np.sum(rigid*source_velocity, axis=1)
    pressure = -WATER_DENSITY*(phi_t+.5*np.sum(total_velocity**2, axis=1)
                              +GRAVITY*panels.centroids[:, 2])
    weighted = panels.areas[:, None]*np.column_stack((panels.normals,
        np.cross(relative, panels.normals)))
    return dict(phi_t=phi_t, pressure=pressure, pressure_functionals=weighted.T@pressure,
        transport_solution=solution, known_source_transport=known_transport,
        diagnostics=dict(current_field_preparation_seconds=field_seconds,
            pressure_solve=solution.diagnostics, total_seconds=perf_counter()-started,
            additional_operator_builds=0, operator='same current body/graph operator',
            surface_trace_ownership='E(chi_exterior)+full known Vbody dot grad(R); known source transport is not refilled',
            source_trace_sampling=dict(
                exterior_gradient_max_difference_m_per_s=float(np.max(np.linalg.norm(derivative_alias[extension.exterior], axis=1))),
                exterior_gradient_rms_difference_m_per_s=float(np.sqrt(np.mean(derivative_alias[extension.exterior]**2))),
                exterior_transport_max_difference_m2_per_s2=float(np.max(abs(transport_alias[extension.exterior]))),
                scope='Partial tangential diagnostic: analytic graph derivative of R minus derivative of projected source trace. Global exterior maximum includes periodic seams; not a bound on full transport forcing error.'),
            variable='Z=phi_t+Vbody dot grad(R); known harmonic-lift transport removed'))
