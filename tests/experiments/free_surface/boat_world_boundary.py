"""One coarse full-boat nonflat boundary, using the exact parent reference.

Actual seed42 incident elevations are sampled at fixed WORLD P1 parent nodes.
The body stays at its independently solved mean-water rest pose. This is a
geometry snapshot, not an initialized incident/scattered field or a simulation.
Fraction arithmetic is a bounded reference, not the production meshing path.

Run: python -m tests.experiments.free_surface.boat_world_boundary --output artifacts/boat-world-boundary-seed42.json
Frozen successful reference: artifacts/boat-world-boundary-seed42-symmetry.{json,npz}.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from topwater import sea_states
from topwater.frames import rotation_matrix
from topwater.hydrostatics import FloodingError, LossOfContactError
from tests.reference_boat import BEAM, GRAVITY, equilibrium, reference_body
from .boat_chart_rebuild import affine, normal_work
from .boat_parent_clip_proof import PlanePoints, GROUPS, columns, hull_planes, cell_polygons, triangulate, sub, cross, dot
from .p4_fluid_ordering import fresh


def rim_clearance(registry, world, rim, cells, roofs):
    """Exact line/parent intersection; includes the virtual graph under the boat."""
    minimum = None
    for a, b in zip(rim, np.roll(rim, -1)):
        start, end = tuple(world[a]), tuple(world[b])
        delta = sub(end, start)
        intervals = []
        for constraints, (roof, _) in zip(cells, roofs):
            lo, hi = F(0), F(1)
            for plane, sign in constraints[2:]:
                value = sign*registry.value(plane, start)
                change = sign*(registry.value(plane, end)-registry.value(plane, start))
                if change > 0:
                    hi = min(hi, -value/change)
                elif change < 0:
                    lo = max(lo, -value/change)
                elif value > 0:
                    hi = -1
            if lo > hi:
                continue
            intervals.append((lo, hi))
            for t in (lo, hi):
                point = tuple(x+t*d for x, d in zip(start, delta))
                clearance = registry.value(roof, point)/registry.planes[roof][2]
                minimum = clearance if minimum is None else min(minimum, clearance)
        covered = F(0)
        for lo, hi in sorted(intervals):
            if lo > covered:
                raise ValueError('rim leaves the fixed parent graph')
            covered = max(covered, hi)
        if covered != 1:
            raise ValueError('rim leaves the fixed parent graph')
    if minimum is None or minimum <= 0:
        raise FloodingError('nonflat graph reaches the actual opening rim')
    return float(minimum)


def parent_weights(point, triangle):
    """Exact affine XY weights; state coefficients remain on fixed parents."""
    a, b = sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0])
    d = sub(point, triangle[0])
    determinant = a[0]*b[1]-a[1]*b[0]
    u, v = (d[0]*b[1]-d[1]*b[0])/determinant, (a[0]*d[1]-a[1]*d[0])/determinant
    return 1-u-v, u, v


def source_weights(point, triangle):
    """Exact original-facet barycentrics; no off-plane projection permitted."""
    a, b, d = sub(triangle[1], triangle[0]), sub(triangle[2], triangle[0]), sub(point, triangle[0])
    normal = cross(a, b)
    if dot(normal, d):
        raise ValueError('material point is not on its original physical facet')
    drop = max(range(3), key=lambda j: abs(normal[j]))
    axes = [j for j in range(3) if j != drop]
    i, j = axes
    determinant = a[i]*b[j]-a[j]*b[i]
    u, v = (d[i]*b[j]-d[j]*b[i])/determinant, (a[i]*d[j]-a[j]*d[i])/determinant
    return 1-u-v, u, v


def material_coordinates(registry_points, rr, cg, position):
    """Invert the SAME rational affine pose, not rounded R.T as an inverse."""
    vectors = cross(rr[1], rr[2]), cross(rr[2], rr[0]), cross(rr[0], rr[1])
    determinant = dot(rr[0], vectors[0])
    return [tuple(cg[j]+sum(vectors[i][j]*(point[i]-position[i]) for i in range(3))/determinant
                  for j in range(3)) for point in registry_points]


def construct(output, *, body=None, pose=None, parent_xy=None, parents=None,
              parent_eta=None, parent_psi=None, export_template=False):
    """Construct one checkpointed snapshot; coefficients stay on WORLD parents.

    Returns ``(mesh, arrays, report)``. Pose translations locate the body's CG;
    orientations use the existing world quaternion convention. Supplied bodies
    must retain the canonical boat's authored 18-triangle physical groups.
    """
    output = fresh(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    report = dict(qualified=False, complete=False, parents_completed=[],
        scope='Full nonflat physical boundary reference only; no PDE, force model, evolved-state transfer or runtime qualification',
        limits=dict(maximum_planes=128, maximum_point_owners=12000, maximum_intersection_cache=50000,
                    soft_parent_budget_seconds=25),
        budget_scope='Elapsed checks between parents; an in-flight parent/edge pass is not preempted',
        source_sha256={name: sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in
                       ('boat_world_boundary.py', 'boat_parent_clip_proof.py')},
        fixture_sha256=sha256(Path(__file__).resolve().parents[2].joinpath('reference_boat.py').read_bytes()).hexdigest())
    arrays, polygons = {}, []
    registry = None
    try:
        supplied_body = body is not None
        body = reference_body() if body is None else body
        rest = equilibrium()[0] if pose is None else pose
        if body.mesh.faces.shape != (18, 3):
            raise ValueError('this reference requires the authored base boat facets')
        rotation = rotation_matrix(rest.orientation)
        source = [[F(float(x)) for x in p] for p in body.mesh.vertices]
        cg = [F(float(x)) for x in body.cg]
        position = [F(float(x)) for x in rest.position]
        rr = [[F(float(x)) for x in row] for row in rotation]
        world = np.array([tuple(sum(rr[j][k]*(p[k]-cg[k]) for k in range(3))+position[j]
                                 for j in range(3)) for p in source], dtype=object)
        if parent_xy is None:
            x, y = np.meshgrid(np.linspace(-6, 6, 5), np.linspace(-3, 3, 3), indexing='ij')
            xy = np.column_stack((x.ravel(), y.ravel()))
        else:
            xy = np.asarray(parent_xy, dtype=float)
        if parents is None:
            if parent_xy is not None:
                raise ValueError('supplied parent nodes require triangle connectivity')
            parents = []
            for i in range(4):
                for j in range(2):
                    a = i*3+j
                    parents.extend(((a, a+3, a+4), (a, a+4, a+1)))
        parents = np.asarray(parents, dtype=np.int64)
        if (xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all()
                or parents.ndim != 2 or parents.shape[1] != 3
                or np.any(parents < 0) or np.any(parents >= len(xy))):
            raise ValueError('invalid fixed parent nodes/connectivity')
        parents = [tuple(row) for row in parents]
        waves = sea_states.realize(beam_m=BEAM, gravity=GRAVITY, seed=42)
        k = waves.angular_frequencies**2/GRAVITY
        phase = (xy[:, :1]*np.cos(waves.directions)+xy[:, 1:]*np.sin(waves.directions))*k+waves.phases
        eta = np.cos(phase)@waves.amplitudes if parent_eta is None else np.asarray(parent_eta, dtype=float)
        if eta.shape != (len(xy),) or not np.isfinite(eta).all():
            raise ValueError('parent_eta must contain one finite height per parent node')
        # A prescribed trace control, not a claim of nonlinear field preparation.
        psi_mean = np.sin(phase)@(GRAVITY/waves.angular_frequencies*waves.amplitudes)
        psi = ((np.sin(phase)*np.exp(eta[:, None]*k))@(GRAVITY/waves.angular_frequencies*waves.amplitudes)
               if parent_psi is None else np.asarray(parent_psi, dtype=float))
        if psi.shape != (len(xy),) or not np.isfinite(psi).all():
            raise ValueError('parent_psi must contain one finite potential per parent node')
        roof = np.array([[F(float(v)) for v in row] for row in np.column_stack((xy, eta))], dtype=object)
        arrays.update(parent_nodes=np.asarray(roof, dtype=float), parents=np.array(parents), parent_psi=psi,
            parent_psi_mean_water=psi_mean,
            body_vertices=np.asarray(world, dtype=float), source_body_vertices=body.mesh.vertices,
            source_body_faces=body.mesh.faces, source_body_cg=body.cg,
            source_body_rim=body.mesh.opening_rim,
            pose_position=rest.position, pose_orientation=rest.orientation)
        report.update(body_mass=body.mass, actual_sea_realization=json.loads(sea_states.dumps_realization(waves)),
            pose=('actual mean-water equilibrium, unchanged by the selected wave snapshot'
                  if pose is None else 'supplied physical CG position and world quaternion'),
            supplied_inputs=dict(body=supplied_body, pose=pose is not None,
                                 parent_xy=parent_xy is not None, parent_eta=parent_eta is not None,
                                 parent_psi=parent_psi is not None),
            potential_trace_control='P1 nodal trace of incident Airy potential at actual graph height; mean-water values saved separately; no nonlinear field initialization claim',
            affine_pose_roundoff_max=float(np.max(abs(np.asarray(world, dtype=float)
                -((body.mesh.vertices-body.cg)@rotation.T+rest.position)))))
        registry = PlanePoints(max_points=12000, max_intersections=50000, capture_constraints=export_template)
        hull = hull_planes(registry, world, body.mesh.faces)
        cells, roofs, outer = columns(registry, roof, parents)
        if len(registry.planes) > 128:
            raise RuntimeError('reference plane cap reached before clipping')
        report['planes'] = len(registry.planes)
        report['rim_clearance'] = rim_clearance(registry, world, body.mesh.opening_rim, cells, roofs)
        bbox_lo = [min(world[:, j]) for j in (0, 1)]
        bbox_hi = [max(world[:, j]) for j in (0, 1)]
        if (min(world[:, 2]) <= -2 or any(bbox_lo[j] <= (-6, -3)[j]
                or bbox_hi[j] >= (6, 3)[j] for j in (0, 1))):
            raise ValueError('actual hull reaches the bed or outer walls')
        body_seen = set()
        solid_volume, column_volume = 0., 0.
        solid_first, column_first = np.zeros(3), np.zeros(3)
        for parent_id, constraints in enumerate(cells):
            if perf_counter()-start > 25:
                raise RuntimeError('reference parent budget exhausted')
            before = perf_counter()
            parent = parents[parent_id]
            points = roof[list(parent)]
            # Strict disjointness only: equality can carry a FINITE body facet.
            far = any(max(points[:, j]) < bbox_lo[j] or min(points[:, j]) > bbox_hi[j] for j in (0, 1))
            pieces, volume, first, duplicates = cell_polygons(registry, hull, constraints,
                outer[parent_id], parent_id, body_seen, outside_bbox=far)
            polygons.extend(pieces)
            solid_volume += volume; solid_first += first
            float_points = np.asarray(points, dtype=float)
            area = abs(np.linalg.det(float_points[1:, :2]-float_points[:1, :2]))/2
            bary = np.full((3, 3), 1/6); np.fill_diagonal(bary, 2/3)
            sample = bary@float_points
            column_volume += area*(sample[:, 2]+2).mean()
            column_first += area*np.r_[(sample[:, :2]*(sample[:, 2:]+2)).mean(axis=0),
                                       ((sample[:, 2]**2-4)/2).mean()]
            report['parents_completed'].append(dict(parent=parent_id, seconds=perf_counter()-before,
                bbox_disjoint=far, polygons=len(pieces), point_owners=len(registry.points), duplicate_body=duplicates))
        report['parent_clipping_seconds'] = perf_counter()-start
        critical = set()
        if export_template:
            # Dry vertices are essential: a new wet region can be born before
            # any currently active wet-boundary vertex changes arrangement.
            for parent, constraints in enumerate(cells):
                if perf_counter()-start > 25:
                    raise RuntimeError('reference critical-witness budget exhausted')
                vertices = registry.polytope([*hull, *constraints[2:]])
                critical.update(vertices)
                for vertex in vertices:
                    registry.relevant[vertex].add(roofs[parent][0])
        mesh, meta = triangulate(registry, polygons)
        triangles = mesh.vertices[mesh.faces]
        tetra = np.einsum('fi,fi->f', triangles[:, 0], np.cross(triangles[:, 1], triangles[:, 2]))/6
        volume = tetra.sum()
        first = (tetra[:, None]*triangles.sum(axis=1)/4).sum(axis=0)
        free = np.flatnonzero(np.isin(meta['planes'], [i for i, _ in roofs])&(meta['kinds'] == 'column'))
        body_faces = np.flatnonzero(meta['kinds'] == 'body')
        bottom_ids = [i for i, labels in enumerate(registry.labels) if any(label[0] == 'bottom' for label in labels)]
        bottom = np.flatnonzero(np.isin(meta['planes'], bottom_ids)&(meta['kinds'] == 'column'))
        if not len(body_faces) or solid_volume <= 0:
            raise LossOfContactError('the reference snapshot has no submerged body')
        walls = np.setdiff1d(np.arange(len(mesh.faces)), np.r_[free, body_faces, bottom])
        support = {}
        for face in free:
            parent_id = meta['source'][face]
            if meta['planes'][face] != roofs[parent_id][0]:
                raise ValueError('free face lost its original parent')
            for node in mesh.faces[face]:
                point = meta['points'][meta['active'][node]]
                weights = parent_weights(point, roof[list(parents[parent_id])])
                assert min(weights) >= 0 and sum(weights) == 1
                assert sum(w*p[2] for w, p in zip(weights, roof[list(parents[parent_id])])) == point[2]
                owner = tuple(sorted((i, w) for i, w in zip(parents[parent_id], weights) if w))
                if node in support:
                    assert support[node] == owner, 'shared node has conflicting parent state weights'
                support[node] = owner
        free_nodes = np.array(sorted(support))
        indices, weights = np.full((len(free_nodes), 3), -1), np.zeros((len(free_nodes), 3))
        for row, node in enumerate(free_nodes):
            for column, (i, value) in enumerate(support[node]):
                indices[row, column], weights[row, column] = i, float(value)
        twist = np.array([.03, -.02, .01, .01, -.015, .02])
        expected_work = np.array([.2, -.1, .3])@(twist[:3]*(column_volume-solid_volume)
            +np.cross(twist[3:], column_first-solid_first-rest.position*(column_volume-solid_volume)))
        work = normal_work(triangles, rest.position, twist, affine)
        volume_error = float(volume-(column_volume-solid_volume))
        first_error = float(np.max(abs(first-(column_first-solid_first))))
        work_error = float(work-expected_work)
        # Geometry roundoff controls, not hydrodynamic accuracy/qualification.
        checks = np.array([volume_error, first_error, work_error])
        if not np.isfinite(checks).all() or np.any(abs(checks) > [1e-10, 1e-9, 1e-10]):
            raise ValueError('independent closed-boundary moment/work check failed')
        arrays.update(points=mesh.vertices, faces=mesh.faces, body_faces=body_faces, free_faces=free,
            wall_faces=walls, bottom_faces=bottom, body_source_group=meta['source'][body_faces],
            free_source_parent=meta['source'][free], free_nodes=free_nodes,
            free_parent_indices=indices, free_parent_weights=weights)
        # Material labels survive independent offset triangulations. Each child
        # corner additionally retains one containing ORIGINAL source triangle.
        material = material_coordinates([meta['points'][i] for i in meta['active']], rr, cg, position)
        source_ids, source_bary = [], []
        for face in body_faces:
            ids, bary = [], []
            for node in mesh.faces[face]:
                for source_id in GROUPS[meta['source'][face]]:
                    weights = source_weights(material[node], [source[i] for i in body.mesh.faces[source_id]])
                    if min(weights) >= 0:
                        ids.append(source_id); bary.append(weights)
                        break
                else:
                    raise ValueError('wet child vertex leaves its original physical facet group')
            source_ids.append(ids); source_bary.append(bary)
        arrays.update(body_material_vertices=np.array([[material[node] for node in mesh.faces[f]]
                                                      for f in body_faces], dtype=float),
                      body_vertex_source_triangle=np.asarray(source_ids, dtype=np.int64),
                      body_vertex_source_barycentric=np.asarray(source_bary, dtype=float))
        if export_template:
            from .boat_geometry_template import export_template as export
            arrays.update(export(registry, meta, arrays, critical))
            report['template_export'] = dict(dry_and_wet_critical_vertices=len(critical),
                source_sha256=sha256(Path(__file__).with_name('boat_geometry_template.py').read_bytes()).hexdigest())
        report.update(complete=True, seconds=perf_counter()-start, nodes=len(mesh.vertices), triangles=len(mesh.faces),
            body_triangles=len(body_faces), free_triangles=len(free), wall_triangles=len(walls), bottom_triangles=len(bottom),
            source_physical_triangle_groups=[list(group) for group in GROUPS[:-1]],
            volume_error=volume_error, first_moment_max_error=first_error,
            affine_rigid_work_error=work_error, parent_height_and_state_inheritance_exact=True,
            geometry_roundoff_checks_passed=True,
            body_polygon_merge=meta['body_merge'],
            geometry_check_absolute_tolerances=dict(volume_m3=1e-10, first_moment_m4=1e-9, affine_rigid_work=1e-10))
    except BaseException as error:
        report['failure'] = dict(type=type(error).__name__, message=str(error))
        raise
    finally:
        if registry is not None:
            arrays['partial_intersections'] = np.array(registry.points, dtype=float).reshape(-1, 3)
            arrays['partial_polygon_vertices'] = np.array([i for polygon, *_ in polygons for i in polygon])
            arrays['partial_polygon_offsets'] = np.r_[0, np.cumsum([len(polygon) for polygon, *_ in polygons])]
            report['point_owners'] = len(registry.points)
            report['intersection_cache'] = len(registry.cache)
        np.savez(output.with_suffix('.npz'), **arrays)
        output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return mesh, arrays, report


def run(output):
    return construct(output)[2]


def body_material_samples(arrays, barycentric):
    """Material XYZ and original group for the chosen hull quadrature rule.

    The returned shapes are (body faces, samples, 3) and (body faces, samples).
    Use these SAME material samples on every offset snapshot, even if its wet
    triangulation differs. No identity between offset nodes is assumed.
    """
    barycentric = np.asarray(barycentric, dtype=float)
    if (barycentric.ndim != 2 or barycentric.shape[1] != 3
            or not np.isfinite(barycentric).all() or np.any(barycentric < 0)
            or not np.allclose(barycentric.sum(axis=1), 1., rtol=0, atol=4*np.finfo(float).eps)):
        raise ValueError('invalid hull quadrature barycentrics')
    points = np.einsum('qa,fai->fqi', barycentric, arrays['body_material_vertices'])
    groups = np.broadcast_to(arrays['body_source_group'][:, None], points.shape[:-1])
    return points, groups


def interpolate_body_trace(arrays, potential, material_points, groups):
    """P1 trace at matching material samples, rejecting missing wet coverage.

    Only same-original-facet-group wet triangles may contain a sample. Boundary
    classification allows coordinate summation roundoff, never extrapolation
    over a physically dry patch. Gauge alignment is the caller's responsibility.
    """
    potential = np.asarray(potential, dtype=float)
    points = np.asarray(material_points, dtype=float)
    groups = np.asarray(groups)
    if (potential.shape != (len(arrays['points']),) or not np.isfinite(potential).all()
            or points.shape[-1:] != (3,) or groups.shape != points.shape[:-1]
            or not np.isfinite(points).all()):
        raise ValueError('invalid material trace samples or nodal potential')
    result = np.full(points.shape[:-1], np.nan)
    flat, labels, values = points.reshape(-1, 3), groups.ravel(), result.ravel()
    triangles = arrays['body_material_vertices']
    face_values = potential[arrays['faces'][arrays['body_faces']]]
    scale = max(1., np.max(abs(arrays['source_body_vertices'])))
    coordinate_roundoff = 64*np.finfo(float).eps*scale
    for group in np.unique(labels):
        if group not in range(len(GROUPS)-1):
            raise ValueError('unknown original physical facet group')
        rows = np.flatnonzero(labels == group)
        source_face = arrays['source_body_faces'][GROUPS[int(group)][0]]
        original = arrays['source_body_vertices'][source_face]
        normal = np.cross(original[1]-original[0], original[2]-original[0])
        normal /= np.linalg.norm(normal)
        if np.any(abs((flat[rows]-original[0])@normal) > coordinate_roundoff):
            raise ValueError('material trace query leaves its original facet plane')
        axes = np.flatnonzero(np.arange(3) != np.argmax(abs(normal)))
        best = np.full(len(rows), -np.inf)
        for face in np.flatnonzero(arrays['body_source_group'] == group):
            projected = triangles[face][:, axes]
            transform = (projected[1:]-projected[0]).T
            inverse = np.linalg.inv(transform)
            uv = (flat[rows][:, axes]-projected[0])@inverse.T
            bary = np.column_stack((1-uv.sum(axis=1), uv))
            # Conditioning scales a COORDINATE roundoff bound into weights;
            # it does not introduce a fixed geometric containment tolerance.
            tolerance = coordinate_roundoff*np.linalg.norm(inverse, ord=np.inf)*2
            minimum = bary.min(axis=1)
            use = (minimum >= -tolerance)&(minimum > best)
            values[rows[use]] = bary[use]@face_values[face]
            best[use] = minimum[use]
    if not np.isfinite(result).all():
        missing = np.flatnonzero(~np.isfinite(values))
        raise ValueError(f'material samples lack wet original-facet coverage: {missing[:8].tolist()}')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    print(json.dumps(run(parser.parse_args().output), indent=2))
