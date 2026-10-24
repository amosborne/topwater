"""A-priori static geometry screen for the unchanged canonical boat.

Uses native clipped-mesh GL2 hydrostatics only: no BEM solve, motion integrator,
or altered hull/mass/stiffness. Execute only under an allocated numerical lease.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
import time
import traceback
from importlib.metadata import version
from pathlib import Path

for thread_variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[thread_variable] = '1'

import numpy as np
import capytaine as cpt
from scipy.io import loadmat
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from threadpoolctl import threadpool_limits


class GeometryConsistencyError(RuntimeError):
    pass


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clip_polygon(points, normal, offset):
    """Clip a planar polygon to normal.dot(point)+offset <= 0."""
    if not len(points):
        return np.empty((0, 3))
    result = []
    previous = points[-1]
    d_previous = float(previous@normal+offset)
    for current in points:
        d_current = float(current@normal+offset)
        if (d_current <= 0) != (d_previous <= 0):
            result.append(previous+(current-previous)*d_previous/(d_previous-d_current))
        if d_current <= 0:
            result.append(current)
        previous, d_previous = current, d_current
    return np.asarray(result).reshape((-1, 3))


def polygon_moments(points):
    """Independent exact area, first and second moments by triangle fans."""
    area, first, second = 0., np.zeros(3), np.zeros((3, 3))
    for index in range(1, len(points)-1):
        triangle = points[[0, index, index+1]]
        a, b, c = triangle
        triangle_area = np.linalg.norm(np.cross(b-a, c-a))/2
        total = triangle.sum(axis=0)
        area += triangle_area
        first += triangle_area*total/3
        second += triangle_area*(np.outer(total, total)+triangle.T@triangle)/12
    return area, first, second


def feature_edges(vertices, faces, threshold):
    edges = {}
    normals = []
    for index, face in enumerate(faces):
        triangle = vertices[face]
        normal = np.cross(triangle[1]-triangle[0], triangle[2]-triangle[0])
        normals.append(normal/np.linalg.norm(normal))
        for a, b in zip(face, np.roll(face, -1)):
            edges.setdefault(tuple(sorted((int(a), int(b)))), []).append(index)
    rim, knuckles = [], []
    for edge, adjacent in sorted(edges.items()):
        if len(adjacent) == 1:
            rim.append(edge)
        elif len(adjacent) == 2:
            angle = np.rad2deg(np.arccos(np.clip(np.dot(normals[adjacent[0]], normals[adjacent[1]]), -1, 1)))
            if angle >= threshold:
                knuckles.append(edge)
        else:
            raise ValueError('Physical hull contains a non-manifold edge')
    return rim, knuckles, np.asarray(normals)


def box_checks():
    """Independent wall-sided box expectations check clipping and wet area."""
    vertices = np.array([[-1., -.5, -.4], [1., -.5, -.4], [1., .5, -.4], [-1., .5, -.4],
                         [-1., -.5, .6], [1., -.5, .6], [1., .5, .6], [-1., .5, .6]])
    faces = [[0, 2, 1], [0, 3, 2], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
             [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]]
    rows = []
    for heave in (-.1, 0., .1):
        posed = vertices+np.array([0., 0., heave])
        mesh = cpt.Mesh(vertices=posed, faces=faces, quadrature_method='Gauss-Legendre 2',
                        auto_clean=False, auto_check=False).immersed_part()
        depth = .4-heave
        checks = dict(volume=bool(abs(mesh.volume-2*depth) < 1e-12),
                      cb=bool(np.max(abs(mesh.center_of_buoyancy-[0., 0., -depth/2])) < 1e-12),
                      waterplane_area=bool(abs(mesh.waterplane_area-2) < 1e-12),
                      wetted_area=bool(abs(mesh.surface_integral(1.)-(2+6*depth)) < 1e-12))
        rows.append(dict(heave_m=heave, checks=checks))
    return rows


def connected_zero_component(rows):
    zero = next(index for index, row in enumerate(rows) if row['heel_deg'] == 0)
    if not rows[zero].get('equilibrium_converged'):
        return []
    left = right = zero
    while left > 0 and rows[left-1].get('equilibrium_converged'):
        left -= 1
    while right+1 < len(rows) and rows[right+1].get('equilibrium_converged'):
        right += 1
    return rows[left:right+1]


def component_logic_checks():
    grid = [index*.5 for index in range(-24, 25)]
    outer_failures = [dict(heel_deg=angle, equilibrium_converged=abs(angle) < 12) for angle in grid]
    inner_gaps = [dict(heel_deg=angle, equilibrium_converged=abs(angle) != .5) for angle in grid]
    zero_failure = [dict(heel_deg=angle, equilibrium_converged=angle != 0) for angle in grid]
    outer = connected_zero_component(outer_failures)
    # Independent toy righting R=Keff*phi has zero secant/tangent error on the
    # valid inner branch; failed outer points cannot reduce its extent to zero.
    angles = np.deg2rad([row['heel_deg'] for row in outer])
    tangent = np.gradient(7*angles, angles, edge_order=2)
    return dict(outer_failures_retain_11_5_degree_inner_branch=outer[0]['heel_deg'] == -11.5 and outer[-1]['heel_deg'] == 11.5,
                toy_linear_inner_tangent_passes=bool(np.max(abs(tangent/7-1)) < 1e-12),
                inner_gaps_do_not_bridge_islands=[row['heel_deg'] for row in connected_zero_component(inner_gaps)] == [0.],
                unresolved_zero_has_no_connected_branch=connected_zero_component(zero_failure) == [])


def make_plots(output, baseline, free_roll, pose_rows, effective_roll, mass, gravity):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plot_keys = ['GZ_m', 'secant_linear_deviation_fraction', 'tangent_deviation_fraction',
                 'wetted_area_change_fraction', 'wetted_set_symmetric_difference_fraction',
                 'waterplane_area_change_fraction', 'heave_m', 'trim_deg']
    valid = [row if row.get('equilibrium_converged') else
             dict(heel_deg=row['heel_deg'], **{key: np.nan for key in plot_keys}) for row in free_roll]
    angles = np.asarray([row['heel_deg'] for row in valid])
    radians = np.deg2rad(angles)
    fig, axes = plt.subplots(3, 2, figsize=(12, 11), layout='constrained')
    axes[0, 0].plot(angles, [row['GZ_m'] for row in valid], label='Exact planar hydrostatics')
    axes[0, 0].plot(angles, effective_roll*radians/(mass*gravity), '--', label='Mean K times heel')
    axes[0, 0].plot(angles, effective_roll*np.sin(radians)/(mass*gravity), ':', label='GM0 sin(heel)')
    axes[0, 0].set(ylabel='World-x righting arm GZ (m)', xlabel='Heel (deg)')
    axes[0, 0].legend()
    axes[0, 1].plot(angles, [100*row['secant_linear_deviation_fraction'] for row in valid], label='Secant error vs mean K')
    axes[0, 1].plot(angles, [100*row['tangent_deviation_fraction'] if row['tangent_deviation_fraction'] is not None else np.nan
                           for row in valid], label='Tangent error vs mean K')
    axes[0, 1].axhline(5, color='grey', linestyle=':')
    axes[0, 1].axhline(10, color='grey', linestyle='--')
    axes[0, 1].set(ylabel='Deviation (%)', xlabel='Heel (deg)')
    axes[0, 1].legend()
    axes[1, 0].plot(angles, [100*row['wetted_area_change_fraction'] for row in valid], label='Net wetted area change')
    axes[1, 0].plot(angles, [100*row['wetted_set_symmetric_difference_fraction'] for row in valid], label='Material wet-set turnover')
    axes[1, 0].plot(angles, [100*row['waterplane_area_change_fraction'] for row in valid], label='Waterplane area change')
    axes[1, 0].set(ylabel='Change relative to mean (%)', xlabel='Heel (deg)')
    axes[1, 0].legend()
    axes[1, 1].plot(angles, [1000*row['heave_m'] for row in valid], label='Free heave (mm)')
    axes[1, 1].plot(angles, [row['trim_deg'] for row in valid], label='Free trim (deg)')
    axes[1, 1].set(ylabel='Free equilibrium (labelled units)', xlabel='Heel (deg)')
    axes[1, 1].legend()
    heave_rows = sorted([row for row in pose_rows if row['valid'] and row['heel_deg'] == 0 and row['trim_deg'] == 0],
                        key=lambda row: row['heave_m'])
    axes[2, 0].plot([row['heave_m'] for row in heave_rows], [row['exact_generalized_load'][2] for row in heave_rows], label='Exact net buoyancy')
    axes[2, 0].plot([row['heave_m'] for row in heave_rows], [row['linear_generalized_load'][2] for row in heave_rows], '--', label='Mean linear K')
    axes[2, 0].set(ylabel='Vertical restoring load (N)', xlabel='Heave (m), heel/trim zero')
    axes[2, 0].legend()
    trim_rows = sorted([row for row in pose_rows if row['valid'] and row['heel_deg'] == 0 and row['heave_m'] == 0],
                       key=lambda row: row['trim_deg'])
    axes[2, 1].plot([row['trim_deg'] for row in trim_rows], [row['exact_generalized_load'][4] for row in trim_rows], label='Exact pitch torque')
    axes[2, 1].plot([row['trim_deg'] for row in trim_rows], [row['linear_generalized_load'][4] for row in trim_rows], '--', label='Mean linear K')
    axes[2, 1].set(ylabel='Pitch restoring load (Nm)', xlabel='Trim (deg), heel/heave zero')
    axes[2, 1].legend()
    for axis in axes.flat:
        axis.grid(True, alpha=.3)
    fig.suptitle('Same canonical hull: a-priori static geometry/restoring screen\n5%/10% lines are engineering labels, not dynamic accuracy proof')
    fig.savefig(output/'geometry-curves.png', dpi=160)
    plt.close(fig)


class GeometryScreen:
    def __init__(self, case, native, controls, started):
        self.case, self.controls, self.started = case, controls, started
        self.rho = case['water']['density_kg_m3']
        self.gravity = case['water']['gravity_m_s2']
        self.mass = case['body_mass_kg']
        self.cg = np.asarray(case['world_cg_m'])
        hull = case['minimum_hull']
        quaternion = hull['equilibrium_quaternion_wxyz']
        rotation = Rotation.from_quat([*quaternion[1:], quaternion[0]]).as_matrix()
        self.vertices = ((np.asarray(hull['vertices_material_m'])-hull['cg_material_m'])@rotation.T+self.cg)
        self.faces = np.asarray([face for index, face in enumerate(hull['triangular_faces'])
                                 if index not in hull['nonphysical_closure_faces']])
        self.rim, self.knuckles, self.normals = feature_edges(
            self.vertices, self.faces, controls['knuckle_dihedral_degrees'])
        self.rim_vertices = sorted(set(vertex for edge in self.rim for vertex in edge))
        self.baseline_polygons = [clip_polygon(self.vertices[face], np.array([0., 0., 1.]), 0.)
                                  for face in self.faces]
        self.baseline_areas = np.asarray([polygon_moments(polygon)[0] for polygon in self.baseline_polygons])
        self.K = np.asarray(native['K'])
        self.length = float(np.ptp(self.vertices[:, 0]))
        self.beam = float(np.ptp(self.vertices[:, 1]))
        self.evaluations = 0

    def evaluate(self, heel_degrees, heave_m, trim_degrees):
        if time.perf_counter()-self.started > self.controls['maximum_elapsed_seconds']:
            raise TimeoutError('Declared geometry screen time budget exceeded')
        self.evaluations += 1
        rotation = Rotation.from_euler('xyz', [heel_degrees, trim_degrees, 0.], degrees=True).as_matrix()
        cg = self.cg+np.array([0., 0., heave_m])
        vertices = (self.vertices-self.cg)@rotation.T+cg
        # The physical opening is a boundary of the hull after artificial deck
        # closure triangles are removed. Never count immersed closure buoyancy.
        rim_clearance = float(vertices[self.rim_vertices, 2].min())
        plane_normal = rotation[2, :]
        plane_offset = float(cg[2]-plane_normal@self.cg)
        wetted_areas, intersections = [], []
        integral_normal_first = 0.
        integral_normal_second = np.zeros(3)
        waterplane_area = 0.
        waterplane_first = np.zeros(3)
        waterplane_second = np.zeros((3, 3))
        for index, face in enumerate(self.faces):
            polygon = clip_polygon(self.vertices[face], plane_normal, plane_offset)
            wetted_areas.append(polygon_moments(polygon)[0])
            common = clip_polygon(polygon, np.array([0., 0., 1.]), 0.)
            intersections.append(polygon_moments(common)[0])
            posed = (polygon-self.cg)@rotation.T+cg
            area, first, second = polygon_moments(posed)
            normal = rotation@self.normals[index]
            integral_normal_first += np.dot(normal, first)/3
            integral_normal_second += normal*np.diag(second)
            waterplane_area -= normal[2]*area
            waterplane_first -= normal[2]*first
            waterplane_second -= normal[2]*second
        wetted_areas, intersections = np.asarray(wetted_areas), np.asarray(intersections)
        turnover = wetted_areas+self.baseline_areas-2*intersections
        edge_classes, contacts = [], []
        for index, edge in enumerate(self.knuckles):
            heights = vertices[list(edge), 2]
            edge_classes.append(-1 if max(heights) < 0 else 1 if min(heights) > 0 else 0)
            if min(heights) <= 0 <= max(heights):
                contacts.append(index)
        q = np.array([0., 0., heave_m, np.deg2rad(heel_degrees), np.deg2rad(trim_degrees), 0.])
        linear = -self.K@q
        result = dict(heel_deg=float(heel_degrees), heave_m=float(heave_m), trim_deg=float(trim_degrees),
                      heave_relative_to_test_water_plane_m=float(heave_m), test_water_plane_z_m=0.,
                      rim_clearance_m=rim_clearance, downflooding=bool(rim_clearance <= 0),
                      vertex_waterline_contacts=np.where(abs(vertices[:, 2]) <= self.controls['contact_tolerance_m'])[0].tolist(),
                      knuckle_waterline_edges=contacts, knuckle_edge_wet_class=edge_classes,
                      wetted_physical_area_m2=float(sum(wetted_areas)),
                      wetted_material_face_area_m2=wetted_areas.tolist(),
                      wetted_set_symmetric_difference_m2=float(sum(turnover)),
                      wetted_set_symmetric_difference_fraction=float(sum(turnover)/sum(self.baseline_areas)),
                      linear_generalized_load=linear.tolist())
        if rim_clearance <= 0:
            result.update(valid=False, invalid_reason='Opening rim at/below water plane; no closed-cap buoyancy accepted')
            return result
        mesh = cpt.Mesh(vertices=vertices, faces=self.faces, quadrature_method='Gauss-Legendre 2',
                        auto_clean=False, auto_check=False).immersed_part()
        volume = float(mesh.volume)
        if mesh.nb_faces == 0 or volume <= 0:
            result.update(valid=False, invalid_reason='No valid immersed body volume')
            return result
        cb = np.asarray(mesh.center_of_buoyancy)
        analytic_cb = integral_normal_second/(2*integral_normal_first)
        native_area = float(mesh.surface_integral(1.))
        # Independent analytic polygon moments check native clipping/GL2, not
        # a second call to the same quadrature implementation.
        checks = dict(volume=bool(abs(volume-integral_normal_first) < 1e-10),
                      center_of_buoyancy=bool(np.max(abs(cb-analytic_cb)) < 1e-10),
                      wetted_area=bool(abs(native_area-sum(wetted_areas)) < 1e-10),
                      waterplane_area=bool(abs(float(mesh.waterplane_area)-waterplane_area) < 1e-10))
        if not all(checks.values()):
            raise GeometryConsistencyError('Native geometry does not match independent analytic moments')
        buoyancy = self.rho*self.gravity*volume
        force = np.array([0., 0., buoyancy-self.mass*self.gravity])
        torque = np.cross(cb-cg, np.array([0., 0., buoyancy]))
        exact = np.r_[force, torque]
        result.update(valid=True, volume_m3=volume, cb_world_m=cb.tolist(), world_cg_m=cg.tolist(),
                      buoyancy_N=buoyancy, exact_generalized_load=exact.tolist(),
                      load_minus_linear=(exact-linear).tolist(),
                      load_error_over_weight_and_weight_length=((exact-linear)/
                        np.array([self.mass*self.gravity]*3+[self.mass*self.gravity*self.beam,
                                  self.mass*self.gravity*self.length, self.mass*self.gravity*self.length])).tolist(),
                      waterplane_area_m2=float(waterplane_area),
                      waterplane_first_xy_m3=waterplane_first[:2].tolist(),
                      waterplane_second_xy_m4=waterplane_second[:2, :2].tolist(),
                      native_vs_analytic_checks=checks, righting_world_x_Nm=float(-torque[0]),
                      righting_about_heel_axis_Nm=float(-np.dot(rotation[:, 0], torque)),
                      equivalent_arm_m=float(-torque[0]/(self.mass*self.gravity)),
                      equivalent_arm_axis='World x; true GZ label only after mass-balanced equilibrium',
                      heel_axis_equivalent_arm_m=float(-np.dot(rotation[:, 0], torque)/(self.mass*self.gravity)))
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', type=Path, default=Path('reference/capytaine/canonical-boat.json'))
    parser.add_argument('--native', type=Path, default=Path('artifacts/reference/minimum-20261024/native-gl2'))
    parser.add_argument('--controls', type=Path, default=Path('reference/capytaine/geometry-envelope-controls.json'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    controls = json.loads(args.controls.read_text())
    if controls['review_status'] != 'ready_for_coordinator_execution_allocation':
        raise ValueError('Geometry controls/source require independent review before allocation')
    case = json.loads(args.case.read_text())
    if sha(args.case) != controls['physical_input_sha256']:
        raise ValueError('Canonical physical identity changed')
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.output/'source'
    source.mkdir()
    for path in (Path(__file__), args.case, args.controls):
        shutil.copyfile(path, source/path.name)
    started = time.perf_counter()
    report = dict(complete=False, scope='Static geometry/buoyancy operating-envelope screen only; no dynamics or BEM',
                  command=[sys.executable, *sys.argv], input_sha256={str(path): sha(path) for path in
                  (Path(__file__), args.case, args.controls, args.native/'native.mat')},
                  versions={name: version(name) for name in ('capytaine', 'numpy', 'scipy', 'threadpoolctl', 'matplotlib')},
                  new_BEM_solves=0, body_simulations=0, native_threads=1)
    rows, fixed_roll, free_roll = [], [], []

    def save():
        (args.output/'manifest.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        (args.output/'geometry-results.json').write_text(json.dumps(dict(fixed_roll=fixed_roll,
            free_roll=free_roll, pose_grid=rows), indent=2, allow_nan=False)+'\n')

    save()
    try:
        native = loadmat(args.native/'native.mat', simplify_cells=True)
        screen = GeometryScreen(case, native, controls, started)
        if not np.array_equal(native['M'], case['dry_mass_matrix_world']) or not np.array_equal(native['cg'], case['world_cg_m']):
            raise ValueError('Frozen full dry mass/CG differs from canonical geometry')
        if version('capytaine') != controls['capytaine_version']:
            raise ValueError('Capytaine version differs from declared geometry controls')
        with threadpool_limits(1):
            component_checks = component_logic_checks()
            if not all(component_checks.values()):
                raise ValueError('Independent connected-coverage interpretation checks failed')
            boxes = box_checks()
            if not all(all(row['checks'].values()) for row in boxes):
                raise ValueError('Independent box clipping/hydrostatic checks failed')
            baseline = screen.evaluate(0., 0., 0.)
            selected = [2, 3, 4]
            K = screen.K[np.ix_(selected, selected)]
            free = [0, 2]  # heave and pitch in this 3x3 block
            effective_roll = float(K[1, 1]-K[1, free]@np.linalg.solve(K[np.ix_(free, free)], K[free, 1]))
            if not baseline['valid'] or effective_roll <= 0:
                raise ValueError('Invalid/non-restoring mean geometry')
            perturbations = [controls['mean_fd_heave_m'], controls['mean_fd_angle_rad'], controls['mean_fd_angle_rad']]
            finite_K = np.empty((3, 3))
            for column, delta in enumerate(perturbations):
                pose = np.zeros(3)
                pose[column] = delta if column == 0 else np.rad2deg(delta)
                # evaluate signature is heel,heave,trim; pose is heave,heel,trim.
                positive = screen.evaluate(pose[1], pose[0], pose[2])
                negative = screen.evaluate(-pose[1], -pose[0], -pose[2])
                finite_K[:, column] = -(np.array(positive['exact_generalized_load'])[selected]-
                                       np.array(negative['exact_generalized_load'])[selected])/(2*delta)
            mean_checks = dict(volume_matches_frozen=bool(abs(baseline['volume_m3']-native['volume']) < 1e-10),
                               cb_matches_frozen=bool(np.max(abs(np.array(baseline['cb_world_m'])-native['cb'])) < 1e-10),
                               finite_difference_K_matches=bool(np.allclose(finite_K, K, rtol=1e-5, atol=1e-3)))
            if not all(mean_checks.values()):
                raise ValueError('Mean full-hull hydrostatics differ from accepted frozen reference')
            report.update(connected_coverage_interpretation_checks=component_checks,
                          independent_box_checks=boxes, mean_checks=mean_checks,
                          geometry_evaluations=screen.evaluations)
            save()
            guesses = {-1: np.zeros(2), 1: np.zeros(2)}
            for heel in sorted(controls['heel_degrees'], key=lambda angle: (abs(angle), angle)):
                fixed_roll.append(screen.evaluate(heel, 0., 0.))
                sign = -1 if heel < 0 else 1

                def residual(values):
                    trial = screen.evaluate(heel, values[0], np.rad2deg(values[1]))
                    if not trial['valid']:
                        raise ValueError(trial['invalid_reason'])
                    return [(trial['volume_m3']-screen.mass/screen.rho)/(screen.mass/screen.rho),
                            (trial['cb_world_m'][0]-screen.cg[0])/screen.length]

                try:
                    solution = least_squares(residual, guesses[sign],
                        bounds=([-controls['free_heave_bound_m'], -np.deg2rad(controls['free_trim_bound_degrees'])],
                                [controls['free_heave_bound_m'], np.deg2rad(controls['free_trim_bound_degrees'])]),
                        max_nfev=controls['maximum_root_evaluations'], xtol=1e-12, ftol=1e-12, gtol=1e-12)
                    result = screen.evaluate(heel, solution.x[0], np.rad2deg(solution.x[1]))
                    result['equilibrium_converged'] = bool(solution.success and np.max(np.abs(residual(solution.x))) < 1e-9)
                    result['equilibrium_residuals'] = residual(solution.x)
                    result['GZ_m'] = result['equivalent_arm_m'] if result['equilibrium_converged'] else None
                    result['GZ_definition'] = 'World-x arm on mass-balanced equilibrium; heel-axis GZ separately reported'
                    result['heel_axis_GZ_m'] = result['heel_axis_equivalent_arm_m'] if result['equilibrium_converged'] else None
                    result['equilibrium_status'] = 'Converged on attempted continuation branch' if result['equilibrium_converged'] else 'Solver/tolerance failure; equilibrium existence unknown'
                    if result['equilibrium_converged']:
                        guesses[sign] = solution.x
                except ValueError as error:
                    result = dict(heel_deg=float(heel), valid=False, equilibrium_converged=False, error=str(error),
                                  equilibrium_status='Invalid solver trial or algorithm failure; equilibrium existence unknown')
                free_roll.append(result)
                if len(free_roll) % 10 == 0:
                    report['geometry_evaluations'] = screen.evaluations
                    save()
            for heel in controls['pose_grid_heel_degrees']:
                for heave in controls['pose_grid_heave_m']:
                    for trim in controls['pose_grid_trim_degrees']:
                        rows.append(screen.evaluate(heel, heave, trim))
            fixed_roll.sort(key=lambda row: row['heel_deg'])
            free_roll.sort(key=lambda row: row['heel_deg'])
            valid = [row for row in free_roll if row.get('equilibrium_converged')]
            complete_free_curve = len(valid) == len(free_roll)
            connected = connected_zero_component(free_roll)
            angles = np.deg2rad([row['heel_deg'] for row in connected])
            restoring = np.asarray([row['righting_world_x_Nm'] for row in connected])
            tangent_values = np.gradient(restoring, angles, edge_order=2) if len(connected) >= 3 else [None]*len(connected)
            tangent_by_angle = {row['heel_deg']: tangent for row, tangent in zip(connected, tangent_values)}
            for row in valid:
                angle, moment = np.deg2rad(row['heel_deg']), row['righting_world_x_Nm']
                on_connected = row['heel_deg'] in tangent_by_angle
                tangent = tangent_by_angle.get(row['heel_deg'])
                if angle == 0 and on_connected:
                    tangent = effective_roll  # independently checked mean derivative limit
                row['effective_mean_roll_stiffness_Nm_rad'] = effective_roll
                row['secant_linear_deviation_fraction'] = float(abs(moment-effective_roll*angle)/(effective_roll*abs(angle))) if angle else 0.
                row['tangent_deviation_fraction'] = float(abs(tangent-effective_roll)/effective_roll) if tangent is not None else None
                row['GZ_over_GM_sin_heel'] = float(moment/(effective_roll*np.sin(angle))) if angle else 1.
                row['tangent_method'] = 'Mean derivative limit' if angle == 0 else 'Finite differences on declared heel grid'
                row['connected_to_mean_equilibrium'] = on_connected
                row['tangent_boundary_rule'] = 'Second-order one-sided at connected-component boundary; centered inside'
                row['restoring_branch'] = bool(moment*angle > 0 or angle == 0)
                row['prospective_restoring_shape_label'] = bool(on_connected and tangent is not None and row['restoring_branch'] and
                    row['secant_linear_deviation_fraction'] <= controls['engineering_assumptions']['secant_fraction'] and
                    row['tangent_deviation_fraction'] <= controls['engineering_assumptions']['tangent_fraction'])
                row['prospective_label_status'] = ('Unknown outside connected mean branch or insufficient tangent coverage'
                    if not on_connected or tangent is None else 'Passes exploratory restoring-shape label'
                    if row['prospective_restoring_shape_label'] else 'Exceeds exploratory label or non-restoring branch')
                row['wetted_area_change_fraction'] = row['wetted_physical_area_m2']/baseline['wetted_physical_area_m2']-1
                row['waterplane_area_change_fraction'] = row['waterplane_area_m2']/baseline['waterplane_area_m2']-1
            by_angle = {row['heel_deg']: row for row in connected}
            parity_rows = []
            for angle in sorted(angle for angle in by_angle if angle > 0 and -angle in by_angle):
                positive, negative = by_angle[angle], by_angle[-angle]
                first, second = np.array(positive['waterplane_second_xy_m4']), np.array(negative['waterplane_second_xy_m4'])
                checks = dict(free_heave_even=abs(positive['heave_m']-negative['heave_m']) < 1e-8,
                              free_trim_even=abs(positive['trim_deg']-negative['trim_deg']) < 1e-6,
                              righting_odd=abs(positive['righting_world_x_Nm']+negative['righting_world_x_Nm']) < 1e-5,
                              cb_y_odd=abs(positive['cb_world_m'][1]+negative['cb_world_m'][1]) < 1e-8,
                              wetted_area_even=abs(positive['wetted_physical_area_m2']-negative['wetted_physical_area_m2']) < 1e-8,
                              waterplane_area_even=abs(positive['waterplane_area_m2']-negative['waterplane_area_m2']) < 1e-8,
                              waterplane_Ixx_Iyy_even=bool(np.max(abs(np.diag(first)[:2]-np.diag(second)[:2])) < 1e-7),
                              waterplane_Ixy_odd=bool(abs(first[0, 1]+second[0, 1]) < 1e-7),
                              wetset_turnover_even=abs(positive['wetted_set_symmetric_difference_m2']-negative['wetted_set_symmetric_difference_m2']) < 1e-8)
                parity_rows.append(dict(heel_magnitude_deg=angle, checks=checks))
            if not all(all(row['checks'].values()) for row in parity_rows):
                raise ValueError('Independent physical reflection parity checks failed')
            half_range = 0.
            if connected:
                for angle in sorted(angle for angle in by_angle if angle > 0 and -angle in by_angle):
                    if by_angle[angle]['prospective_restoring_shape_label'] and by_angle[-angle]['prospective_restoring_shape_label']:
                        half_range = angle
                    else:
                        break
            make_plots(args.output, baseline, free_roll, rows, effective_roll, screen.mass, screen.gravity)
            report.update(complete=True, elapsed_seconds=time.perf_counter()-started,
                          geometry_evaluations=screen.evaluations, mean_checks=mean_checks,
                          connected_coverage_interpretation_checks=component_checks,
                          independent_box_checks=boxes, independent_parity_checks=parity_rows,
                          mean_finite_difference_K=finite_K.tolist(), effective_mean_roll_stiffness_Nm_rad=effective_roll,
                          baseline=baseline, free_equilibria_converged=len(valid), free_equilibria_requested=len(free_roll),
                          free_curve_coverage_complete=complete_free_curve,
                          connected_mean_equilibrium_heel_degrees=[row['heel_deg'] for row in connected],
                          engineering_assumptions=controls['engineering_assumptions'],
                          connected_prospective_restoring_shape_half_range_deg=half_range,
                          target_5deg_passes_prospective_restoring_shape_label=bool(half_range >= 5),
                          target_10deg_passes_prospective_restoring_shape_label=bool(half_range >= 10),
                          feature_edges=dict(rim=screen.rim, knuckles=screen.knuckles),
                          physical_limit='Actual opening immersion or non-restoring branch rejects that physical state; failed root leaves equilibrium existence unknown. Chine contact alone is diagnostic.',
                          interpretation='A-priori still-water geometry/restoring screen only. Does not establish finite-wave diffraction/radiation/excitation linearity or dynamic trajectory accuracy. For future buoy sea interpret immersion relative to local water elevation and body-versus-wave slopes; large global heave alone is not flooding or failure when hull follows eta. Calm-water heel curves are conservative fixed-mean-assumption diagnostics, not a nonlinear wave-surface model. No phase/shape fitting.')
    except BaseException as error:
        report.update(error=dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()),
                      elapsed_seconds=time.perf_counter()-started)
        if 'screen' in locals():
            report['geometry_evaluations'] = screen.evaluations
    save()
    print(json.dumps({key: report.get(key) for key in ('complete', 'geometry_evaluations', 'elapsed_seconds', 'error')}))
    return 0 if report['complete'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
