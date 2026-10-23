"""Fresh planar-boat contact charts across flat-water chine transitions.

Prescribed WORLD trace functions are independent of current triangle IDs. They
are transfer controls, not a claimed transfer scheme for a solved wave history.
No volume mesher, fluid operator, LU, dynamics or nonflat-water clipping is used.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import InitialState
from topwater.frames import rotation_matrix
from topwater.hydrostatics import hydrostatics
from tests.laplace_reference import triangle_samples, affine, normal_work
from tests.reference_boat import GRAVITY, WATER_DENSITY
from .boat_contact_preflight import make_chart, pose
from .p4_fluid_ordering import fresh




def smooth(points):
    return 1.7+.2*np.cos(.3*np.asarray(points)[..., 0]+.2*np.asarray(points)[..., 1])




def diagnostics(body, state, chart):
    points = chart.reference_points
    triangles = points[chart.faces]
    independent = hydrostatics(body, state, density=WATER_DENSITY, gravity=GRAVITY)
    area_normal = np.cross(triangles[:, 1]-triangles[:, 0],
                           triangles[:, 2]-triangles[:, 0])/2
    fluid_volume = float(np.einsum('fi,fi->', area_normal, triangles.mean(axis=1))/3)
    wet = triangles[chart.body_faces]
    # Verify each original physical triangle's wet area, not just the total.
    area = np.linalg.norm(area_normal[chart.body_faces], axis=1)
    native_area = np.linalg.norm(np.cross(independent.wetted_triangles[:, 1]-independent.wetted_triangles[:, 0],
                                          independent.wetted_triangles[:, 2]-independent.wetted_triangles[:, 0]), axis=1)/2
    actual_parent = np.bincount(chart.body_source_faces, weights=area, minlength=len(body.mesh.faces))
    expected_parent = np.bincount(independent.source_faces, weights=native_area, minlength=len(body.mesh.faces))
    twist = np.array([.03, -.02, .01, .01, -.015, .02])
    current = chart.evaluate(state.position, state.orientation, np.zeros(len(chart.free_nodes)))
    velocity = current.velocity(np.zeros(len(chart.free_nodes)), twist)
    # Integrate the chart velocity and nodal affine trace independently by the
    # exact P1 triangle polynomial formula. Normal orientation is fluid-outward.
    nodes = chart.faces[chart.body_faces]
    nodal_phi = affine(points[nodes])
    nodal_flux = np.einsum('fai,fi->fa', velocity[nodes], area_normal[chart.body_faces])
    chart_work = np.sum((nodal_phi.sum(axis=1)*nodal_flux.sum(axis=1)
                         +np.sum(nodal_phi*nodal_flux, axis=1))/12)
    exact_work = -normal_work(independent.wetted_triangles, state.position, twist, affine)
    # Constant gauge shift changes the potential work by C times net rigid flux;
    # do not subtract means or silently impose a different gauge on a fresh mesh.
    shift = .6
    shifted_phi = nodal_phi+shift
    shifted_work = np.sum((shifted_phi.sum(axis=1)*nodal_flux.sum(axis=1)
                           +np.sum(shifted_phi*nodal_flux, axis=1))/12)
    net_flux = -normal_work(independent.wetted_triangles, state.position, twist,
                           lambda x: np.ones(x.shape[:-1]))
    # P1 trace interpolation's actual error against the same smooth WORLD
    # function, rather than calling direct resampling exact state transfer.
    free_triangles = triangles[chart.free_faces]
    sample, weights = triangle_samples(free_triangles, 8)
    barycentric, _ = triangle_samples(np.eye(3)[None], 8)
    nodal = smooth(free_triangles)
    interpolated = nodal@barycentric[0].T
    area = np.linalg.norm(area_normal[chart.free_faces], axis=1)
    trace_rms = np.sqrt(np.einsum('f,q,fq->', area, weights, (interpolated-smooth(sample))**2)/area.sum())
    return dict(fluid_volume=fluid_volume, expected_fluid_volume=12*6*2-independent.displaced_volume,
        fluid_volume_error=fluid_volume-(12*6*2-independent.displaced_volume),
        wet_parent_area_max_error=float(np.max(abs(actual_parent-expected_parent))),
        wet_area_error=float(actual_parent.sum()-independent.wetted_area),
        affine_body_normal_work=float(chart_work), independent_body_normal_work=exact_work,
        normal_work_error=float(chart_work-exact_work),
        constant_gauge_work_error=float(shifted_work-chart_work-shift*net_flux),
        prescribed_world_smooth_trace_interpolation_area_rms=float(trace_rms),
        rim_clearance_flat_water=independent.rim_clearance,
        physical_parent_ids=np.flatnonzero(expected_parent).tolist())


def run(output):
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(qualified=False, complete=False, cases=[],
        scope='Flat-water fresh contact geometry and prescribed-world-trace controls only; no evolved-state transfer or fluid solve',
        surface_triangulation='Constrained PSLG; exact input-point bijection and all prescribed boundary edges; no Steiner points',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
            ('boat_chart_rebuild.py', 'boat_contact_preflight.py', 'canoe_contact.py', 'canoe.py')},
        fixture_sha256=sha256(Path(__file__).resolve().parents[2].joinpath('reference_boat.py').read_bytes()).hexdigest())
    arrays = {}
    try:
        before = perf_counter()
        body, rest, frozen = make_chart()
        report['cold_reference_construction_seconds'] = perf_counter()-before
        report['reference_body'] = dict(mass=body.mass, cg=body.cg.tolist(),
            rest_position=rest.position.tolist(), rest_orientation=rest.orientation.tolist())
        arrays.update(reference_points=frozen.reference_points, reference_faces=frozen.faces)
        for name, translation, rotation in (
            ('rest', np.zeros(3), np.zeros(3)),
            ('pitch_0.001', np.zeros(3), np.array([0., .001, 0.])),
            ('heave_0.01', np.array([0., 0., .01]), np.zeros(3)),
            ('roll_0.01', np.zeros(3), np.array([.01, 0., 0.])),
            ('pitch_0.01', np.zeros(3), np.array([0., .01, 0.]))):
            position, orientation = pose(rest, translation, rotation)
            state = InitialState(position, orientation)
            row = dict(name=name)
            try:
                frozen.evaluate(position, orientation, np.zeros(len(frozen.free_nodes)))
                row['frozen_chart_valid'] = True
            except Exception as error:
                row.update(frozen_chart_valid=False, frozen_error_type=type(error).__name__, frozen_error=str(error))
            try:
                times = []
                for _ in range(2):
                    before = perf_counter()
                    _, _, chart = make_chart(state=state)
                    times.append(perf_counter()-before)
                row.update(complete=True, construction_seconds=times, nodes=len(chart.reference_points),
                    free_nodes=len(chart.free_nodes), body_triangles=len(chart.body_faces),
                    contact_nodes=len(chart.contact_nodes),
                    physical_vertex_displacement_max=float(np.max(np.linalg.norm(
                        (body.mesh.vertices-body.cg)@(rotation_matrix(orientation)-rotation_matrix(rest.orientation)).T
                        +position-rest.position, axis=1))), **diagnostics(body, state, chart))
                arrays[name+'_points'] = chart.reference_points
                arrays[name+'_faces'] = chart.faces
                arrays[name+'_free_nodes'] = chart.free_nodes
                arrays[name+'_world_trace'] = smooth(chart.reference_points[chart.free_nodes])
            except Exception as error:
                row.update(complete=False, error_type=type(error).__name__, message=str(error))
            report['cases'].append(row)
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
