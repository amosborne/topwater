"""Batched version of reviewed interval geometry; original scalar oracle retained."""
import time
import numpy as np
from actual_wave_geometry_core import AnalyticSea, BudgetExceeded, area_bounds


class BatchSea(AnalyticSea):
    def __init__(self, bank, gravity, budget):
        super().__init__(bank, gravity, budget)
        self.constant_roundoff_scale = 4096*np.finfo(float).eps*(1+float(np.sum(abs(self.a))))
        self.maximum_phase = float(np.max(abs(self.phase)))
        self.maximum_omega = float(np.max(self.w))
        self.maximum_k = float(np.max(np.linalg.norm(self.k, axis=1)))
        self.timings = dict(field_seconds=0., clipping_seconds=0., adaptation_seconds=0.)

    def field(self, xy, t):
        started = time.perf_counter()
        result = super().field(xy, t)
        self.timings['field_seconds'] += time.perf_counter()-started
        return result

    def roundoff_margin(self, xy, t):
        return self.constant_roundoff_scale*(1+self.maximum_phase+self.maximum_omega*abs(t)
            +self.maximum_k*float(np.max(np.linalg.norm(xy, axis=-1))))

    def triangle_margins(self, xy, t):
        return self.constant_roundoff_scale*(1+self.maximum_phase+self.maximum_omega*abs(t)
            +self.maximum_k*np.max(np.linalg.norm(xy, axis=2), axis=1))


def batch_clip(rows, counts, column, thresholds):
    """Same Sutherland-Hodgman order/<= rule as scalar clip; padded capacity6."""
    size, capacity, width = rows.shape
    index = np.arange(capacity)[None, :]
    valid = index < counts[:, None]
    following = np.take_along_axis(rows, np.broadcast_to(
        ((index+1) % np.maximum(counts[:, None], 1))[:, :, None], rows.shape), axis=1)
    threshold = np.broadcast_to(np.asarray(thresholds), (size,))[:, None]
    first_inside = rows[:, :, column] <= threshold
    second_inside = following[:, :, column] <= threshold
    crossing = valid & (first_inside != second_inside)
    fraction = np.divide(threshold-rows[:, :, column], following[:, :, column]-rows[:, :, column],
                         out=np.zeros((size, capacity)), where=crossing)
    intersection = rows+fraction[:, :, None]*(following-rows)
    candidates = np.stack((rows, intersection), axis=2).reshape(size, 2*capacity, width)
    mask = np.stack((valid & first_inside, crossing), axis=2).reshape(size, 2*capacity)
    ranks = np.cumsum(mask, axis=1)-1
    output_counts = mask.sum(axis=1)
    if np.any(output_counts > capacity):
        raise ValueError('Convex clipping exceeded reviewed padded vertex capacity')
    result = np.zeros_like(rows)
    row_ids = np.broadcast_to(np.arange(size)[:, None], mask.shape)
    result[row_ids[mask], ranks[mask]] = candidates[mask]
    return result, output_counts


def batch_polygon_area(rows, counts):
    size, capacity, _ = rows.shape
    index = np.arange(capacity)[None, :]
    xyz = rows[:, :, :3]-rows[:, :1, :3]
    following = np.take_along_axis(xyz, np.broadcast_to(
        ((index+1) % np.maximum(counts[:, None], 1))[:, :, None], xyz.shape), axis=1)
    area_vector = np.where((index < counts[:, None])[:, :, None], np.cross(xyz, following), 0).sum(axis=1)/2
    return np.where(counts >= 3, np.linalg.norm(area_vector, axis=1), 0.)


def batch_area_bounds(mean, f, epsilon):
    size = len(mean)
    rows = np.zeros((size, 6, 5))
    rows[:, :3, :3], rows[:, :3, 3], rows[:, :3, 4] = mean, f, mean[:, :, 2]
    counts = np.full(size, 3)
    lower, lower_count = batch_clip(rows, counts, 3, -epsilon)
    upper, upper_count = batch_clip(rows, counts, 3, epsilon)
    mean_rows, mean_count = batch_clip(rows, counts, 4, 0.)
    low_intersection, low_intersection_count = batch_clip(lower, lower_count, 4, 0.)
    up_intersection, up_intersection_count = batch_clip(upper, upper_count, 4, 0.)
    L, U = batch_polygon_area(lower, lower_count), batch_polygon_area(upper, upper_count)
    M = batch_polygon_area(mean_rows, mean_count)
    IL, IU = batch_polygon_area(low_intersection, low_intersection_count), batch_polygon_area(up_intersection, up_intersection_count)
    total = np.linalg.norm(np.cross(mean[:, 1]-mean[:, 0], mean[:, 2]-mean[:, 0]), axis=1)/2
    turn_lower = np.maximum(0., M-IU+L-IL)
    turn_upper = np.minimum(total, np.maximum(turn_lower, M-IL+U-IU))
    return np.column_stack((L, U, turn_lower, turn_upper))


def evaluate(mean, posed, face, depth, sea, t, plane):
    sea.budget.check(len(mean))
    field = sea.field(posed[:, :, :2].reshape(-1, 2), t).reshape(-1, 3, 3)
    started = time.perf_counter()
    edges = posed[:, [1, 2, 0], :2]-posed[:, [0, 1, 2], :2]
    diameter2 = np.max(np.sum(edges*edges, axis=2), axis=1)
    margin = sea.triangle_margins(posed[:, :, :2], t)
    epsilon = sea.C2*diameter2/2+margin
    bounds = batch_area_bounds(mean, posed[:, :, 2]-field[:, :, 0], epsilon)
    departure = field[:, :, 0]-(plane[0]+(posed[:, :, :2]-plane[1:3])@plane[3:5])
    sampled = np.max(abs(departure), axis=1)
    curvature = np.column_stack((np.maximum(0., sampled-margin), sampled+epsilon))
    sea.timings['clipping_seconds'] += time.perf_counter()-started
    return dict(mean=mean, posed=posed, face=face, depth=depth, bounds=bounds, curvature=curvature)


def split_batch(triangles):
    a, b, c = triangles[:, 0], triangles[:, 1], triangles[:, 2]
    ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
    # Parent-major ordering, matched by np.repeat for face/depth arrays.
    return np.stack((np.stack((a, ab, ca), axis=1), np.stack((ab, b, bc), axis=1),
                     np.stack((ca, bc, c), axis=1), np.stack((ab, bc, ca), axis=1)), axis=1).reshape(-1, 3, 3)


def summarize(leaves, faces):
    bounds = np.zeros((faces, 4))
    np.add.at(bounds, leaves['face'], leaves['bounds'])
    return dict(wet_area_bounds_m2=bounds[:, :2].sum(axis=0).tolist(),
        wet_set_turnover_bounds_m2=bounds[:, 2:].sum(axis=0).tolist(),
        per_material_face_bounds_m2=bounds.tolist(),
        footprint_curvature_departure_bounds_m=np.max(leaves['curvature'], axis=0).tolist(),
        footprint_curvature_status='ANALYTIC_SPATIAL_BOUND_WITH_REPORTED_WIDTH_NO_TARGET_ADOPTED',
        leaf_count=len(leaves['face']), maximum_depth=int(np.max(leaves['depth'])))


def wetted_bounds_batch(mean, posed, sea, t, plane, tolerance, maximum_depth):
    faces = len(mean)
    leaves = evaluate(mean, posed, np.arange(faces), np.zeros(faces, dtype=int), sea, t, plane)
    previous = None
    while True:
        started = time.perf_counter()
        result = summarize(leaves, faces)
        area, turn = result['wet_area_bounds_m2'], result['wet_set_turnover_bounds_m2']
        result['area_target_met'] = max(area[1]-area[0], turn[1]-turn[0]) <= tolerance
        result['previous_iteration_bounds'] = previous
        if result['area_target_met']:
            result['status'] = 'BOUNDED_TO_DECLARED_SPATIAL_TOLERANCE'
            sea.timings['adaptation_seconds'] += time.perf_counter()-started
            return result
        widths = leaves['bounds'][:, [1, 3]]-leaves['bounds'][:, [0, 2]]
        eligible = np.flatnonzero((leaves['depth'] < maximum_depth) & (widths.max(axis=1) > 1e-15))
        if len(eligible) == 0:
            result['status'] = 'UNKNOWN_TARGET_NOT_MET_AT_DEPTH_CAP'
            sea.timings['adaptation_seconds'] += time.perf_counter()-started
            return result
        selected = eligible[np.argsort(widths[eligible].sum(axis=1))[::-1][:256]]
        if len(leaves['face'])+3*len(selected) > 65536:
            result['status'] = 'UNKNOWN_TARGET_NOT_MET_AT_ACTIVE_LEAF_MEMORY_CAP'
            sea.timings['adaptation_seconds'] += time.perf_counter()-started
            return result
        child_mean, child_posed = split_batch(leaves['mean'][selected]), split_batch(leaves['posed'][selected])
        child_faces, child_depths = np.repeat(leaves['face'][selected], 4), np.repeat(leaves['depth'][selected]+1, 4)
        sea.timings['adaptation_seconds'] += time.perf_counter()-started
        try:
            children = evaluate(child_mean, child_posed, child_faces, child_depths, sea, t, plane)
        except BudgetExceeded:
            result['status'] = 'UNKNOWN_BUDGET_EXHAUSTED_WITH_VALID_COARSE_BOUNDS'
            return result
        started = time.perf_counter()
        keep = np.ones(len(leaves['face']), dtype=bool)
        keep[selected] = False
        leaves = {key: np.concatenate((value[keep], children[key])) for key, value in leaves.items()}
        previous = dict(wet_area_bounds_m2=area, wet_set_turnover_bounds_m2=turn)
        sea.timings['adaptation_seconds'] += time.perf_counter()-started


def batch_oracle_fixtures(mean_triangles):
    # Static inputs only, no synthetic trajectory: parallel/contact/all-wet/dry,
    # ordinary cuts, equal-valued vertices and translated/degenerate triangles.
    base = np.array([[0., 0., -.5], [1., 0., .5], [0., 1., .5]])
    examples = [base, base+np.array([3., -2., 1.]), np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]])]
    examples.extend(mean_triangles)
    mean, f, epsilon = [], [], []
    for triangle in examples:
        for values in ([0., 0., 0.], [-1., -1., -1.], [1., 1., 1.], [-1., 0., 1.],
                       [0., 0., 1.], [0., 0., -1.], [1e-13, -1e-13, 0.], [-.3, .7, -.2]):
            for margin in (0., 1e-12, .1):
                mean.append(triangle)
                f.append(values)
                epsilon.append(margin)
    mean, f, epsilon = np.asarray(mean), np.asarray(f), np.asarray(epsilon)
    actual = batch_area_bounds(mean, f, epsilon)
    expected = np.array([area_bounds(a, a, b, e) for a, b, e in zip(mean, f, epsilon)])
    np.testing.assert_allclose(actual, expected, rtol=3e-12, atol=3e-12)
    assert np.isfinite(actual).all() and np.all(actual[:, 0] <= actual[:, 1]+3e-12)
    assert np.all(actual[:, 2] <= actual[:, 3]+3e-12)
    return dict(cases=len(mean), scalar_oracle_exact_algorithm=True,
                maximum_abs_error=float(np.max(abs(actual-expected))), fixtures='degenerate/parallel/contact/flat/material triangles')
