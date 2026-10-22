"""Small 2500 kg boat contact-chart preflight; no volume mesh or fluid solve.

Uses the existing planar-facet chart unchanged. Failures are retained, including
real chine/vertex transitions: this diagnostic never remeshes or relaxes guards.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.spatial.transform import Rotation

from topwater import sea_states
from tests.reference_boat import BEAM, GRAVITY, equilibrium, reference_body
from .canoe_contact import CanoeContact
from .domain import tank
from .p4_canoe_contact import P4CanoeContact
from .p4_fluid_ordering import fresh


def make_chart(*, enriched=False):
    """Explicit 12 m x 6 m x 2 m surface box around the actual loaded boat.

    This box is solely a contact-geometry diagnostic, not an accepted wave tank.
    No tetrahedra, volume ALE factor or fluid matrix are built.
    """
    body = reference_body()
    if body.mass != 2500:
        raise ValueError('the updated 2500 kg reference boat is required')
    rest, _ = equilibrium()
    box, top, top_nodes = tank(12, 6, 2, depth=2.)
    points = box.vertices*np.array([6/(2*np.pi), 3/np.pi, 1.])
    chart = CanoeContact.make(body, rest, points[box.faces[~top]],
                              points[top_nodes.ravel()], subdivisions=1)
    return body, rest, P4CanoeContact.make(chart) if enriched else chart


def target_heights(points, waves):
    """Actual catalog realization's deep-water eta at t=0, without recentering."""
    k = waves.angular_frequencies**2/GRAVITY
    phase = (points[:, :1]*np.cos(waves.directions)
             +points[:, 1:2]*np.sin(waves.directions))*k+waves.phases
    return np.cos(phase)@waves.amplitudes


def pose(rest, displacement, rotation):
    orientation = (Rotation.from_rotvec(rotation)
                   *Rotation.from_quat(np.roll(rest.orientation, -1))).as_quat()
    return rest.position+displacement, np.roll(orientation, 1)


def run(output):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    waves = sea_states.realize(beam_m=BEAM, gravity=GRAVITY, seed=42)
    report = dict(qualified=False, complete=False, cases=[], derivative_checks=[], construction=[],
        scope='Existing contact geometry only; no fluid PDE, dynamics or operating-envelope qualification',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('boat_contact_preflight.py', 'canoe_contact.py', 'p4_canoe_contact.py')},
        fixture_sha256=sha256(Path(__file__).resolve().parents[2].joinpath('reference_boat.py').read_bytes()).hexdigest(),
        sea_realization= json.loads(sea_states.dumps_realization(waves)))
    arrays = {}
    try:
        for enriched in (False, True):
            label = 'P4' if enriched else 'P1'
            before = perf_counter()
            try:
                body, rest, chart = make_chart(enriched=enriched)
            except Exception as error:
                report['construction'].append(dict(chart=label, complete=False,
                    seconds=perf_counter()-before, error_type=type(error).__name__, message=str(error)))
                continue
            physical = chart.linear if enriched else chart
            height = target_heights(chart.reference_points[chart.free_nodes], waves)
            pattern = height/max(abs(height))
            report['construction'].append(dict(chart=label, complete=True, seconds=perf_counter()-before,
                nodes=len(chart.reference_points), free_nodes=len(chart.free_nodes),
                contact_nodes=len(chart.contact_nodes), body_triangles=len(physical.body_faces),
                free_triangles=len(physical.free_faces), reference_target_height_min=float(height.min()),
                reference_target_height_max=float(height.max())))
            arrays[label+'_reference_points'] = chart.reference_points
            arrays[label+'_free_nodes'] = chart.free_nodes
            arrays[label+'_target_heights'] = height
            # Authored small geometry perturbations, not rescaled catalog cases.
            cases = [('rest', np.zeros(3), np.zeros(3), np.zeros_like(height))]
            for value in (.001, .01):
                for sign in (-1, 1):
                    cases.extend([(f'heave_{sign*value}', np.array([0., 0., sign*value]), np.zeros(3), np.zeros_like(height)),
                        (f'roll_{sign*value}', np.zeros(3), np.array([sign*value, 0., 0.]), np.zeros_like(height)),
                        (f'pitch_{sign*value}', np.zeros(3), np.array([0., sign*value, 0.]), np.zeros_like(height)),
                        (f'target_shape_max_{sign*value}_m', np.zeros(3), np.zeros(3), sign*value*pattern)])
            cases.extend([('six_axis_small', np.full(3, .001), np.array([.001, -.001, .001]), .001*pattern),
                          ('actual_catalog_seed42_t0', np.zeros(3), np.zeros(3), height)])
            for name, translation, rotation, eta in cases:
                before = perf_counter()
                row = dict(chart=label, name=name)
                try:
                    position, orientation = pose(rest, translation, rotation)
                    state = chart.evaluate(position, orientation, eta)
                    row.update(complete=True, maximum_node_displacement=float(
                        np.max(np.linalg.norm(state.points-chart.reference_points, axis=1))))
                except Exception as error:
                    row.update(complete=False, error_type=type(error).__name__, message=str(error))
                row['seconds'] = perf_counter()-before
                report['cases'].append(row)
            base = chart.evaluate(rest.position, rest.orientation, np.zeros_like(height))
            for direction in range(7):
                twist = np.eye(6)[direction] if direction < 6 else np.zeros(6)
                eta_direction = np.zeros_like(height) if direction < 6 else pattern
                exact = base.velocity(eta_direction, twist)
                for h in (1e-6, 5e-7):
                    row = dict(chart=label, direction=direction, step=h)
                    try:
                        values = []
                        for sign in (-1, 1):
                            p, q = pose(rest, sign*h*twist[:3], sign*h*twist[3:])
                            values.append(chart.evaluate(p, q, sign*h*eta_direction).points)
                        difference = (values[1]-values[0])/(2*h)-exact
                        row.update(complete=True, maximum_absolute_error=float(np.max(abs(difference))),
                            relative_l2_error=float(np.linalg.norm(difference)/np.linalg.norm(exact)))
                    except Exception as error:
                        row.update(complete=False, error_type=type(error).__name__, message=str(error))
                    report['derivative_checks'].append(row)
        report['complete'] = True
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.output), indent=2))
