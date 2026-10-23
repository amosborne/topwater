"""Independent h/h2 parent-action verification; not a moving solver branch."""
import numpy as np
from scipy.spatial.transform import Rotation
from topwater import InitialState
from topwater.constants import GRAVITY, WATER_DENSITY
from .boat_geometry_template import InvalidTemplate
from .boat_world_field import boundary_data
from .boat_parent_action import parent_mass, factor_parent_mass, port_transport_covector, port_seam_load


def action_rates(arrays, field, gradient, values, twist, time, geometry, *,
                 combine_incident_ports=False, port_covector=None):
    """Candidate all-parent rates from one CURRENT fixed-nodal energy gradient.

    geometry(pose, eta) must replay the same owned/refined chart at fixed psi.
    Both h/h2 results are returned; this does not choose a physical tolerance.
    """
    points, fixed = arrays['points'], arrays['free_nodes']
    if gradient.shape != points.shape or not np.isfinite(gradient).all():
        raise ValueError('finite current vertex gradient required')
    eta, psi = arrays['parent_nodes'][:, 2], arrays['parent_psi']
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    mass = parent_mass(arrays)
    solve, condition = factor_parent_mass(mass)
    eta_rate = solve(field['parent_map'].T@field['flux'])
    geometry_evaluations = 0
    def checked(current_pose, height):
        nonlocal geometry_evaluations
        geometry_evaluations += 1
        current = geometry(current_pose, height)
        for key in ('faces', 'free_nodes', 'free_faces', 'body_faces', 'wall_faces', 'bottom_faces', 'node_chart_key'):
            if not np.array_equal(current[key], arrays[key]):
                raise ValueError('parent derivative leaves the owned chart/basis')
        return current
    for halving in range(13):
        dt = 1e-5/2**halving
        try:
            timed = []
            for offset in (dt, -dt, dt/2, -dt/2):
                orientation = (Rotation.from_rotvec(offset*twist[3:])
                    *Rotation.from_quat(pose.orientation, scalar_first=True)).as_quat(scalar_first=True)
                timed.append(checked(InitialState(pose.position+offset*twist[:3], orientation), eta+offset*eta_rate))
        except InvalidTemplate:
            continue
        break
    else:
        raise ValueError('no valid owned mass-transport stencil')
    time_mass = [parent_mass(current) for current in timed]
    cdot = [(time_mass[0]-time_mass[1])/(2*dt), (time_mass[2]-time_mass[3])/dt]
    covector = (port_transport_covector(arrays, field['potential'], values, time)
        if port_covector is None else np.asarray(port_covector))
    if covector.shape != points.shape or not np.isfinite(covector).all():
        raise ValueError('finite current port covector required')
    count = len(eta)
    names = ('kinetic', 'body_load', 'moving_trace', 'mass_transport', 'gravity')
    names += ('incident_port',) if combine_incident_ports else ('incident_load', 'port_work')
    terms = {name: np.zeros((2, count)) for name in names}
    steps = []
    for index in range(count):
        direction = np.eye(count)[index]
        for halving in range(13):
            h = 1e-5/2**halving
            try:
                moved = [checked(pose, eta+offset*direction) for offset in (h, -h, h/2, -h/2)]
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no valid owned parent-height stencil')
        steps.append(h)
        loads = [boundary_data(current, values, twist, time=time,
            incident_load=not combine_incident_ports) for current in moved]
        masses = [parent_mass(current) for current in moved]
        for level in range(2):
            step = h/2**level
            plus, minus = moved[2*level:2*level+2]
            bp, bm = loads[2*level:2*level+2]
            moving = (plus['points']-minus['points'])/(2*step)
            dbody = ((bp['body_load']-bm['body_load'])@twist)/(2*step)
            dincident = (bp['right']-bm['right'])/(2*step)-dbody
            dtrace = ((bp['parent_map']-bm['parent_map'])@psi)/(2*step)
            dc = (masses[2*level]-masses[2*level+1])/(2*step)
            terms['kinetic'][level, index] = -np.sum(gradient*moving)
            terms['body_load'][level, index] = WATER_DENSITY*(field['potential']@dbody)
            terms['moving_trace'][level, index] = -WATER_DENSITY*(field['flux']@dtrace)
            terms['gravity'][level, index] = -WATER_DENSITY*GRAVITY*(mass@eta)[index]
            terms['mass_transport'][level, index] = WATER_DENSITY*(psi@dc@eta_rate-(cdot[level].T@psi)[index])
            if combine_incident_ports:
                terms['incident_port'][level, index] = -WATER_DENSITY*np.sum(covector*moving)
            else:
                terms['incident_load'][level, index] = WATER_DENSITY*(field['potential']@dincident)
                terms['port_work'][level, index] = -WATER_DENSITY*(field['potential']@dincident+np.sum(covector*moving))
    corrected = sum(terms.values())
    raw = None if combine_incident_ports else corrected-terms['port_work']
    return dict(eta_rate=eta_rate, psi_rate=np.stack([solve(load/WATER_DENSITY) for load in corrected]),
        parent_mass=mass, mass_time_derivative=np.stack(cdot), raw_action_weak=raw,
        corrected_action_weak=corrected, terms=terms, height_steps_m=np.array(steps),
        transport_step_seconds=dt, geometry_evaluations=geometry_evaluations,
        parent_mass_condition=float(condition), port_seam_weak=port_seam_load(arrays, field, values, time),
        port_seam_weak_order6=port_seam_load(arrays, field, values, time, order=6))
