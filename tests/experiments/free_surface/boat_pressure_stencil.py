"""One owner for prescribed-state operator comparisons on the canonical boat.

Offsets are supplied explicitly: this compares operators on a common physical
direction, rather than advancing each model's own free-surface trajectory.
"""
import numpy as np

from topwater import InitialState
from topwater.frames import rotation_matrix
from tests.reference_boat import WATER_DENSITY, GRAVITY
from .boat_world_boundary import body_material_samples, interpolate_body_trace
from .boat_world_field import current_field, surface_rates
from .body_pressure import evaluate_body_pressure, quadrature
from .boat_ale_pressure import differentiate


def evaluate(body, arrays, state, values, time, energy, offsets, offset_energies, *, charted=False):
    pose = InitialState(state[:3], state[3:7]/np.linalg.norm(state[3:7]))
    twist = state[7:13]
    if (not np.isfinite(energy).all()
            or np.linalg.norm(energy-energy.T) > 1e-10*np.linalg.norm(energy)):
        raise ValueError('nonfinite or asymmetric predicted energy')
    field = current_field(arrays, energy, values, twist, time=time)
    eta, psi, _ = surface_rates(arrays['points'], arrays['faces'][arrays['free_faces']],
        arrays['free_nodes'], field['parent_map'], arrays['parent_psi'], field['flux'])
    samples, groups = body_material_samples(arrays, quadrature()[0])
    hints = np.broadcast_to(np.arange(len(arrays['body_faces']))[:, None], groups.shape)
    rotation = rotation_matrix(pose.orientation)
    mesh_velocity = twist[:3]+np.cross(twist[3:], arrays['points']-pose.position)
    pressures = []
    for h in sorted((dt for dt in offsets if dt > 0), reverse=True):
        traces = []
        for dt in (h, -h):
            moved = offsets[dt]
            moved_field = current_field(moved, offset_energies[dt], values, twist,
                                        acceleration_columns=False, time=time+dt)
            traces.append(moved_field['potential'] if charted else
                interpolate_body_trace(moved, moved_field['potential'], samples, groups, face_hints=hints))
        diagnostics = {}
        if charted:
            rate, mesh_velocity, diagnostics = differentiate(arrays, offsets[h], offsets[-h],
                traces[0], traces[1], h, pose.position, twist)
        else:
            rate = (traces[0]-traces[1])/(2*h)
        pressure = evaluate_body_pressure(arrays['points'], arrays['faces'], arrays['body_faces'],
            pose.position, body.mass, rotation@body.inertia_cg@rotation.T, twist,
            field['potential'], rate, mesh_velocity,
            field['acceleration_potentials'], density=WATER_DENSITY, gravity=GRAVITY)
        pressures.append(dict(step=h, load_world=pressure.hydrodynamic_load.tolist(),
            acceleration_world=pressure.acceleration.tolist(), added_mass=pressure.added_mass.tolist(),
            derivative_representation='charted nodal' if charted else 'material quadrature',
            derivative_diagnostics=diagnostics))
    return field, eta, psi, pressures
