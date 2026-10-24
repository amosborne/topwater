"""Body action derivatives for frozen audits and the opt-in moving driver.

World translation is the CG position; rotations act on the left about that CG.
The caller owns the CURRENT energy gradient and the matching offset fields.
"""
import numpy as np
from scipy.spatial.transform import Rotation

from topwater import InitialState
from tests.reference_boat import GRAVITY, WATER_DENSITY
from topwater.surface import surface_integrals
from .boat_geometry_template import InvalidTemplate
from .boat_parent_action import parent_mass, port_transport_covector
from .boat_world_field import boundary_data


def fluid_gravity(arrays):
    volume, center = surface_integrals(arrays['points'], arrays['faces'])
    return WATER_DENSITY*GRAVITY*volume*center[2]


def body_derivatives(arrays, field, gradient, values, twist, time, eta_rate, geometry, *,
                     combine_incident_ports=False, port_covector=None):
    """Complete fixed-world pose derivatives at h and h/2, no fluid solves."""
    pose = InitialState(arrays['pose_position'], arrays['pose_orientation'])
    rotation = Rotation.from_quat(pose.orientation, scalar_first=True)
    eta, psi = arrays['parent_nodes'][:, 2], arrays['parent_psi']
    covector = (port_transport_covector(arrays, field['potential'], values, time)
        if port_covector is None else np.asarray(port_covector))
    if covector.shape != arrays['points'].shape or not np.isfinite(covector).all():
        raise ValueError('finite current port covector required')
    names = ('kinetic', 'boundary_load', 'moving_trace', 'surface_one_form', 'gravity')
    names += ('incident_port',) if combine_incident_ports else ('port_work',)
    terms = {name: np.zeros((2, 6)) for name in names}
    steps, evaluations = [], 0
    for index in range(6):
        for halving in range(13):
            h = 1e-5/2**halving
            moved = []
            try:
                for offset in (h, -h, h/2, -h/2):
                    direction = offset*np.eye(3)[index % 3]
                    current_pose = (InitialState(pose.position+direction, pose.orientation) if index < 3
                        else InitialState(pose.position, (Rotation.from_rotvec(direction)*rotation).as_quat(scalar_first=True)))
                    evaluations += 1
                    current = geometry(current_pose, eta)
                    for key in ('faces', 'free_nodes', 'free_faces', 'body_faces', 'wall_faces', 'bottom_faces', 'node_chart_key'):
                        if not np.array_equal(current[key], arrays[key]):
                            raise ValueError('body derivative leaves the owned chart/basis')
                    moved.append(current)
            except InvalidTemplate:
                continue
            break
        else:
            raise ValueError('no valid owned body derivative stencil')
        steps.append(h)
        loads = [boundary_data(current, values, twist, time=time,
            incident_load=not combine_incident_ports) for current in moved]
        masses = [parent_mass(current) for current in moved]
        gravity = [fluid_gravity(current) for current in moved]
        for level in range(2):
            step = h/2**level
            plus, minus = moved[2*level:2*level+2]
            bp, bm = loads[2*level:2*level+2]
            moving = (plus['points']-minus['points'])/(2*step)
            dright = (bp['right']-bm['right'])/(2*step)
            dincident = dright-(bp['body_load']-bm['body_load'])@twist/(2*step)
            dtrace = (bp['parent_map']-bm['parent_map'])@psi/(2*step)
            dc = (masses[2*level]-masses[2*level+1])/(2*step)
            terms['kinetic'][level, index] = -np.sum(gradient*moving)
            terms['boundary_load'][level, index] = WATER_DENSITY*(field['potential']@dright)
            terms['moving_trace'][level, index] = -WATER_DENSITY*(field['flux']@dtrace)
            terms['surface_one_form'][level, index] = WATER_DENSITY*(psi@dc@eta_rate)
            terms['gravity'][level, index] = -(gravity[2*level]-gravity[2*level+1])/(2*step)
            if combine_incident_ports:
                terms['incident_port'][level, index] = -WATER_DENSITY*np.sum(covector*moving)
            else:
                terms['port_work'][level, index] = -WATER_DENSITY*(field['potential']@dincident+np.sum(covector*moving))
    return dict(terms=terms, pose_load=np.sum(list(terms.values()), axis=0),
                steps=np.array(steps), geometry_evaluations=evaluations)


def momentum(field):
    """F_V, including the full nodal Neumann map and its contact rows."""
    return -WATER_DENSITY*(field['body_load'].T@field['potential'])


def known_load(derivatives, field, twist, plus, minus, step, level):
    """Zero-body-acceleration action load; includes time-varying B and phi."""
    transport = (momentum(plus)-momentum(minus))/(2*step)
    angular = np.r_[np.zeros(3), -np.cross(twist[3:], momentum(field)[3:])]
    return derivatives['pose_load'][level]+transport+angular, transport, angular
