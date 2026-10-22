"""Static canoe coupling on a complete finite fluid boundary; not a wave solver."""
import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import solve
from scipy.spatial import ConvexHull, Delaunay
from threadpoolctl import threadpool_limits

from topwater import HullMesh, SimulationInputs
from topwater.hydrostatics import equilibrate_heave
from topwater.rankine import TrianglePanels
from tests.reference_hull import input_data, surface_integrals
from .domain import tank
from .operator import assemble_energy
from .coupling import check_flow


def fluid_domain(wet, resolution=(8, 4, 4), domain_scale=1.):
    """Join the canonical convex flat waterline to a rectangular tank.

    This unconstrained Delaunay construction is deliberately restricted to the
    canonical canoe. Closed-mesh validation rejects lost waterline edges; no
    tolerance welding or general moving-waterline remeshing is provided.
    """
    box, top, node_map = tank(*resolution)
    points = box.vertices*domain_scale
    walls = points[box.faces[~top]]
    return join_wetted_hull(wet.wetted_triangles, wet.waterline_segments,
                            walls, points[node_map.ravel()])


def join_wetted_hull(wetted_triangles, segments, walls, grid):
    """Join the canonical convex flat waterline to supplied local box surfaces.

    Shared by the all-BEM reference and the spatial-pocket experiment. This
    remains a restricted mesh construction, not a general remeshing interface.
    """
    line = np.unique(segments.reshape(-1, 3), axis=0)
    if not np.all(line[:, 2] == 0):
        raise ValueError('this experiment requires a horizontal waterline at z=0')
    if np.any(line[:, :2].min(axis=0) <= grid[:, :2].min(axis=0)) or np.any(
            line[:, :2].max(axis=0) >= grid[:, :2].max(axis=0)):
        raise ValueError('waterline must lie strictly inside the tank')
    equations = ConvexHull(line[:, :2]).equations
    midpoints = segments[:, :, :2].mean(axis=1)
    distance = np.max(midpoints@equations[:, :2].T+equations[:, 2], axis=1)
    if np.any(np.abs(distance) > 64*np.finfo(float).eps*np.ptp(line, axis=0).max()):
        raise ValueError('this experiment requires a single convex waterline')
    outside = np.max(grid[:, :2]@equations[:, :2].T+equations[:, 2], axis=1) > 0
    surface = np.unique(np.vstack((grid[outside], line)), axis=0)
    triangles = surface[Delaunay(surface[:, :2]).simplices]
    centers = triangles.mean(axis=1)
    keep = np.max(centers[:, :2]@equations[:, :2].T+equations[:, 2], axis=1) > 0
    triangles = triangles[keep]
    reverse = np.cross(triangles[:, 1]-triangles[:, 0],
                       triangles[:, 2]-triangles[:, 0])[:, 2] < 0
    triangles[reverse] = triangles[reverse, ::-1]
    # Physical wet hull is body-outward; the complete domain is fluid-outward.
    all_triangles = np.concatenate((walls, wetted_triangles[:, ::-1], triangles))
    vertices, indices = np.unique(all_triangles.reshape(-1, 3), axis=0,
                                  return_inverse=True)
    mesh = HullMesh(vertices, indices.reshape(-1, 3))
    body = np.arange(len(walls), len(walls)+len(wetted_triangles))
    free = np.arange(body[-1]+1, len(mesh.faces))
    return mesh, body, free


def body_load(mesh, panels, body, center):
    """Exact P1 weak load for six rigid motions, including waterline nodes."""
    points = mesh.vertices[mesh.faces[body]]
    normal = panels.normals[body, None]
    modes = np.concatenate((np.broadcast_to(normal, points.shape),
                            np.cross(points-center, normal)), axis=2)
    load = np.zeros((len(mesh.vertices), 6))
    total = modes.sum(axis=1)
    for corner in range(3):
        values = panels.areas[body, None]*(total+modes[:, corner])/12
        np.add.at(load, mesh.faces[body, corner], values)
    return load


def solve_mixed(energy, fixed, trace, normal_load):
    """Prescribed nodal potential plus weak non-free-surface flux; one factor."""
    unknown = np.setdiff1d(np.arange(len(energy)), fixed)
    rhs = normal_load[unknown]
    if np.any(trace):
        rhs -= energy[np.ix_(unknown, fixed)]@trace
    potential = np.empty_like(normal_load)
    potential[fixed] = trace
    potential[unknown] = solve(energy[np.ix_(unknown, unknown)], rhs, assume_a='pos')
    return potential


def solve_inertia(mesh, body, free, center, order=3):
    """One mixed factorization: six body modes plus four independent fields."""
    panels = TrianglePanels.from_mesh(mesh)
    energy, diagnostics = assemble_energy(mesh.faces, panels, order)
    fixed = np.unique(mesh.faces[free])
    unknown = np.setdiff1d(np.arange(len(mesh.vertices)), fixed)
    load = body_load(mesh, panels, body, center)
    # phi=x,y,z,1: matching body translations and outer-wall normal flux.
    # Keep all body/wall load rows, including shared free-surface nodes.
    walls = np.setdiff1d(np.arange(len(mesh.faces)), np.r_[body, free])
    fields = np.column_stack((mesh.vertices, np.ones(len(mesh.vertices))))
    velocity = np.zeros((6, 4))
    velocity[:3, :3] = np.eye(3)
    field_load = load@velocity
    for corner in range(3):
        np.add.at(field_load[:, :3], mesh.faces[walls, corner],
                  panels.areas[walls, None]*panels.normals[walls]/3)
    rhs = np.column_stack((load, field_load))
    trace = np.zeros((len(fixed), 10))
    trace[:, 6:] = fields[fixed]
    potential = solve_mixed(energy, fixed, trace, rhs)
    added_mass = 1000*load.T@potential[:, :6]
    volume, _ = surface_integrals(mesh.vertices, mesh.faces)
    depth = -mesh.vertices[:, 2].min()
    affine = dict(
        potential_max_depth_normalized=float(np.max(np.abs(
            potential[unknown, 8]-mesh.vertices[unknown, 2]))/depth),
        energy_relative_error=float((field_load[:, 2]@potential[:, 8]-volume)/volume),
        exact_trace_energy_relative_error=float(
            (mesh.vertices[:, 2]@energy@mesh.vertices[:, 2]-volume)/volume),
        pressure_load=(-load.T@potential[:, 8]).tolist(),
    )
    diagnostics['coupling'] = check_flow(mesh, panels, body, free, center, energy,
                                         load, potential[:, 6:], field_load)
    return added_mass, affine, diagnostics


def compare(refinement=4, resolution=(8, 4, 4), domain_scale=1., order=3):
    from tests.capytaine_reference import instantaneous_inertia

    data = SimulationInputs.from_json(json.dumps(input_data(refinement)))
    start = perf_counter()
    rest, wet = equilibrate_heave(data.body, data.initial_state)
    mesh, body, free = fluid_domain(wet, resolution, domain_scale)
    actual, affine, diagnostics = solve_inertia(mesh, body, free, rest.position, order)
    seconds = perf_counter()-start
    # Fresh reference, outside candidate timing, with identical physical input.
    expected, damping = instantaneous_inertia(wet.wetted_triangles, rest.position, data.body.mass)
    assert np.all(damping == 0)
    scale = np.array([1., 1., 1., 3.5, 3.5, 3.5])
    normalized = (actual-expected)/scale[:, None]/scale[None, :]/(1000*3.5**3)
    return dict(qualified=False, refinement=refinement, resolution=resolution,
        domain_scale=domain_scale, quadrature_order=order, faces=len(mesh.faces),
        vertices=len(mesh.vertices), wetted_faces=len(body), free_faces=len(free),
        density=1000., center_of_mass=rest.position.tolist(), capytaine='3.0.0',
        candidate_seconds=seconds, actual=actual.tolist(), reference=expected.tolist(),
        diagonal_relative_error=(np.diag(actual-expected)/np.diag(expected)).tolist(),
        max_normalized_matrix_error=float(np.max(np.abs(normalized))),
        affine=affine, affine_pressure_reference=(np.r_[wet.force, wet.moment]/9810).tolist(),
        diagnostics=diagnostics,
        scope='finite rigid tank versus open infinite-depth reference; not a correlation gate')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', action='store_true',
                        help='also vary outer spacing, domain size and quadrature')
    args = parser.parse_args()
    cases = [(3, (8, 4, 4), 1., 3), (4, (8, 4, 4), 1., 3)]
    if args.study:
        cases.extend([(4, (16, 8, 4), 1., 3), (4, (16, 8, 8), 1., 3),
                      (4, (32, 16, 8), 2., 3), (4, (8, 4, 4), 1., 6)])
    rows = []
    destination = Path('artifacts')
    destination.mkdir(exist_ok=True)
    with threadpool_limits(1):
        for case in cases:
            row = compare(*case)
            rows.append(row)
            (destination/'free-surface-canoe.json').write_text(
                json.dumps(rows, indent=2)+'\n', encoding='utf-8')
            print(json.dumps(row), flush=True)


if __name__ == '__main__':
    main()
