"""Old/new source-matrix parity and component ABBA on actual convex boats.

The literal old builder is a verification oracle only. Current velocity and
prescribed-pressure parity use existing physics owners, not a second solver.
"""
import argparse
from contextlib import contextmanager
from hashlib import sha256
import json
from pathlib import Path
import shutil
from time import perf_counter
import numpy as np
from topwater import InitialState
from topwater.rankine import TrianglePanels, single_layer_trace
from topwater.numerics import LaplaceControls
from topwater.frames import rotation_matrix
from .boat_spatial_work import source_matrix
from .panel_gradient_reference import checked_influence


def old_source_matrix(panels):
    """Unchanged pre-8e31da0 per-row implementation, independent of trace API."""
    count = len(panels)
    matrix = np.empty((count, count))
    for row in range(count):
        other = np.arange(count) != row
        subset = TrianglePanels(panels.triangles[other], panels.normal_sign[other])
        real = checked_influence(panels.centroids[row:row+1], subset)
        matrix[row, other] = real.gradient_single_layer[0]@panels.normals[row]
        matrix[row, row] = .5
    mirror = panels.centroids*np.array([1., 1., -1.])
    for first in range(0, count, 128):
        last = min(first+128, count)
        image = checked_influence(mirror[first:last], panels)
        matrix[first:last] -= np.einsum('tnj,tj->tn', image.gradient_single_layer,
                            panels.normals[first:last]*np.array([1., 1., -1.]))
    return matrix


def difference(actual, expected):
    error = float(np.max(abs(actual-expected)))
    relative = float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected), 1e-30))
    if not np.isfinite(actual).all():
        raise ValueError('nonfinite computational parity result')
    return dict(max_absolute=error, relative_l2=relative)


def matrix_accuracy_request(panels):
    """Two qualified gS routes, component requests contracted with |n_target|.

    This uses the existing requested integration accuracy, not a certified
    physical force bound. The identical image calculation cancels in comparison.
    """
    controls = LaplaceControls()
    budget = np.empty((len(panels), len(panels)))
    for first in range(0, len(panels), controls.target_block_size):
        last = min(first+controls.target_block_size, len(panels))
        gradient = single_layer_trace(panels, np.arange(first, last)).gradient_single_layer
        component = 2*(controls.quadrature_atol+controls.quadrature_rtol*abs(gradient))/(1-controls.quadrature_rtol)
        budget[first:last] = np.einsum('tnj,tj->tn', component, abs(panels.normals[first:last]))
    np.fill_diagonal(budget, 0.)
    return budget


def geometry(inputs, time, grid=65, length=64.):
    from .boat_graph_geometry import constant_twist_pose, wet_snapshot
    from .boat_hos_continuation import tilted_geometry, tangent_graph
    from .boat_graph_images import image_height
    initial = inputs.initial_state
    pose = (InitialState(initial.position, initial.orientation) if time == 0.
            else constant_twist_pose(initial, time))
    _, level, slope = tilted_geometry(inputs)
    normal = np.r_[-slope, 1.]; normal /= np.linalg.norm(normal)
    tangent = np.array([1., 0., slope[0]]); tangent /= np.linalg.norm(tangent)
    frame = np.vstack((tangent, np.cross(normal, tangent), normal))
    origin = initial.position[:2]-length/2
    x, y = np.meshgrid(origin[0]+np.arange(grid)*length/grid, origin[1]+np.arange(grid)*length/grid)
    xy = np.column_stack((x.ravel(), y.ravel()))
    panels, inside, _, _ = wet_snapshot(inputs.body, pose, np.r_[initial.position[:2], level],
                                       frame, xy, length, 1)
    eta = tangent_graph(inputs, x, y, level, slope)[0]
    image = image_height(panels.triangles, eta, [length]*2, origin)[0]
    image = max(image, float(np.max(panels.triangles[..., 2]+level
                    +(panels.triangles[..., :2]-initial.position[:2])@slope))/2)
    return dict(pose=pose, panels=panels, inside=inside, eta=eta, origin=origin,
                xy=xy, length=length, image=float(np.nextafter(image, np.inf)),
                height=float(panels.triangles[..., 2].max()))


def matrix_audit(data):
    triangles = data['panels'].triangles-[0., 0., data['image']]
    signs = data['panels'].normal_sign
    def timed(builder):
        before = perf_counter()
        matrix = builder(TrianglePanels(triangles, signs))
        return matrix, perf_counter()-before
    old, first_old = timed(old_source_matrix)
    new, first_new = timed(source_matrix)
    parity = difference(new, old)
    requested = matrix_accuracy_request(TrianglePanels(triangles, signs))
    rounding = 64*np.finfo(float).eps*np.maximum(1., abs(old))
    if np.any(abs(new-old) > requested+rounding):
        raise ValueError('matrix difference exceeds existing contracted gS accuracy request')
    parity['maximum_fraction_of_existing_request'] = float(np.max(abs(new-old)/(requested+rounding)))
    rows = []
    for name, builder in [('old', old_source_matrix), ('new', source_matrix),
                          ('new', source_matrix), ('old', old_source_matrix)]:
        matrix, elapsed = timed(builder)
        difference(matrix, old)
        rows.append(dict(builder=name, seconds=elapsed))
    return dict(panel_count=len(triangles), matrix_shape=list(old.shape),
        panel_sha256=sha256(triangles.tobytes()+signs.tobytes()).hexdigest(),
        matrix_parity=parity, diagonal_max_difference=float(np.max(abs(np.diag(new-old)))),
        first_use_seconds=dict(old=first_old, new=first_new), abba=rows,
        timing_scope='Fresh TrianglePanels plus COMPLETE real and unchanged image assembly; geometry preparation separately counted. First-use includes JIT if needed. Component only, not full-RHS throughput.')


@contextmanager
def builder_owner(module, builder):
    previous = module.source_matrix
    module.source_matrix = builder
    try:
        yield
    finally:
        module.source_matrix = previous


def pressure_stage(inputs, data, builder):
    from . import boat_graph_bvp as module
    from .cuda_graph import ExteriorTrace, cp
    from .boat_transport_pressure import prescribed_pressure
    before = perf_counter()
    with builder_owner(module, builder):
        context = module.GraphBody(data['panels'], data['eta'], data['length'], data['origin'],
                                   data['height'], data['image'], 4)
    extension = ExteriorTrace(data['inside'].reshape(data['eta'].shape))
    ext = .1*np.cos(data['xy'][extension.exterior]@(2*np.pi*np.array([3., 1.])/data['length']))
    trace = extension.extend(ext)
    pose, panels = data['pose'], data['panels']
    omega = rotation_matrix(pose.orientation)@pose.angular_velocity
    velocity = pose.linear_velocity+np.cross(omega, panels.centroids-pose.position)
    solution = context.solve(trace, np.sum(velocity*panels.normals, axis=1))
    # One nonzero forced boundary derivative, same data for both builders.
    chi = .02*np.sin(data['xy'][extension.exterior]@(2*np.pi*np.array([3., 1.])/data['length']))
    pressure = prescribed_pressure(context, solution.density, solution.spectrum, pose,
                                   np.array([.01, -.02, .03, .005, -.004, .003]), extension, chi)
    cp.cuda.get_current_stream().synchronize()
    values = dict(density=solution.density, body_potential=solution.potential,
                  phi_t=pressure['phi_t'], pressure=pressure['pressure'],
                  all_six_loads=pressure['pressure_functionals'])
    return values, dict(total_seconds=perf_counter()-before, velocity=solution.diagnostics,
                       pressure=pressure['diagnostics'],
                       scope='Complete current velocity BVP and prescribed six-acceleration pressure stage; forced graph/trace, no evolved trajectory or free-body RHS qualification.')


def run(output, pressure=True):
    from tests.reference_boat import simulation_inputs
    from tests.reference_convex_boat import box_inputs
    result = dict(complete=False, qualified=False, scope=__doc__, cases=[])
    output.write_text(json.dumps(result))
    for name, inputs in [('canonical', simulation_inputs()), ('box', box_inputs())]:
        output.with_suffix('.source').joinpath(name+'-inputs.json').write_text(inputs.to_json())
        for time in [0., .15]:
            before = perf_counter()
            data = geometry(inputs, time)
            row = dict(boat=name, prescribed_time=time, geometry_seconds=perf_counter()-before,
                       input_sha256=sha256(inputs.to_json().encode()).hexdigest(), **matrix_audit(data))
            if pressure and time:
                old, old_diagnostics = pressure_stage(inputs, data, old_source_matrix)
                new, new_diagnostics = pressure_stage(inputs, data, source_matrix)
                row['pressure_stage_parity'] = {key: difference(new[key], old[key]) for key in old}
                row['pressure_stages'] = dict(old=old_diagnostics, new=new_diagnostics)
            result['cases'].append(row)
            output.write_text(json.dumps(result, indent=2, allow_nan=False))
            print(json.dumps(row), flush=True)
    result['complete'] = True
    output.write_text(json.dumps(result, indent=2, allow_nan=False))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--matrix-only', action='store_true', help='CPU source-matrix controls only')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('fresh output required')
    archive = args.output.with_suffix('.source')
    archive.mkdir(parents=True, exist_ok=False)
    # Freeze all local experiment owners rather than omit a transitive physical helper.
    sources = [*Path('src/topwater').glob('*.py'),
               *[p for p in Path(__file__).parent.iterdir() if p.suffix in ('.py', '.cu', '.cuh')],
               *[Path('tests')/name for name in ('reference_boat.py', 'reference_convex_boat.py',
                                                'reference_hull.py', 'laplace_reference.py', 'benchmark_rankine.py')]]
    hashes = {}
    for source in sources:
        relative = source.resolve().relative_to(Path.cwd())
        target = archive/relative; target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[str(relative)] = sha256(source.read_bytes()).hexdigest()
    result = run(args.output, pressure=not args.matrix_only)
    result['producer_sha256'] = hashes
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False))
