"""Geometry-only analytic-wave bounds; no forces, integration of motion or BEM."""
import time
import numpy as np


class BudgetExceeded(RuntimeError):
    pass


class Budget:
    def __init__(self, seconds, triangles):
        self.started = time.perf_counter()
        self.seconds, self.maximum = seconds, triangles
        self.triangles, self.points = 0, 0

    def check(self, triangles=0):
        if time.perf_counter()-self.started > self.seconds or self.triangles+triangles > self.maximum:
            raise BudgetExceeded('Declared time/triangle budget exhausted; no tolerance relaxation')
        self.triangles += triangles


class AnalyticSea:
    def __init__(self, bank, gravity, budget):
        self.a = np.asarray(bank['amplitudes'])
        self.w = np.asarray(bank['angular_frequencies'])
        self.phase = np.asarray(bank['phases'])
        theta = np.asarray(bank['directions'])
        assert all(v.shape == (384,) and np.isfinite(v).all() for v in (self.a, self.w, self.phase, theta))
        self.k = self.w[:, None]**2/gravity*np.column_stack((np.cos(theta), np.sin(theta)))
        self.C2 = float(np.sum(abs(self.a)*np.sum(self.k*self.k, axis=1)))
        self.budget = budget

    def roundoff_margin(self, xy, t):
        phase_scale = 1+float(np.max(abs(self.phase)))+float(np.max(self.w))*abs(t)
        phase_scale += float(np.max(np.linalg.norm(self.k, axis=1)))*float(np.max(np.linalg.norm(xy, axis=-1)))
        return 4096*np.finfo(float).eps*(1+float(np.sum(abs(self.a))))*phase_scale

    def field(self, xy, t):
        xy = np.asarray(xy).reshape(-1, 2)
        self.budget.check()
        output = []
        for start in range(0, len(xy), 2048):
            self.budget.check()
            psi = self.k@xy[start:start+2048].T+self.phase[:, None]-self.w[:, None]*t
            output.append(np.column_stack((self.a@np.cos(psi),
                (-(self.a[:, None]*self.k).T@np.sin(psi)).T)))
        self.budget.points += len(xy)
        return np.concatenate(output)


def polygon_area(points):
    points = np.asarray(points)
    if len(points) < 3:
        return 0.
    return float(sum(np.linalg.norm(np.cross(points[i]-points[0], points[i+1]-points[0]))/2
                     for i in range(1, len(points)-1)))


def clip(rows, column, threshold):
    """Rows contain materialxyz, affine current f, and affine meanplane z."""
    if len(rows) == 0:
        return np.empty((0, 5))
    result = []
    for first, second in zip(rows, np.roll(rows, -1, axis=0)):
        inside_first, inside_second = first[column] <= threshold, second[column] <= threshold
        if inside_first:
            result.append(first)
        if inside_first != inside_second:
            fraction = (threshold-first[column])/(second[column]-first[column])
            result.append(first+fraction*(second-first))
    return np.asarray(result).reshape(-1, 5)


def area_bounds(mean, posed, f, epsilon):
    rows = np.column_stack((mean, f, mean[:, 2]))
    total = polygon_area(mean)
    lower, upper = clip(rows, 3, -epsilon), clip(rows, 3, epsilon)
    L, U = polygon_area(lower[:, :3]), polygon_area(upper[:, :3])
    mean_area = polygon_area(clip(rows, 4, 0)[:, :3])
    IL = polygon_area(clip(lower, 4, 0)[:, :3])
    IU = polygon_area(clip(upper, 4, 0)[:, :3])
    # Disjoint sets: meanwet/currentdry and meandry/currentwet.
    turn_lower = max(0., mean_area-IU+L-IL)
    turn_upper = min(total, max(turn_lower, mean_area-IL+U-IU))
    return np.array([L, U, turn_lower, turn_upper])


def split(triangle):
    a, b, c = triangle
    ab, bc, ca = (a+b)/2, (b+c)/2, (c+a)/2
    return (np.array([a, ab, ca]), np.array([ab, b, bc]),
            np.array([ca, bc, c]), np.array([ab, bc, ca]))


def evaluate_leaves(entries, sea, t, plane):
    sea.budget.check(len(entries))
    posed = np.stack([entry['posed'] for entry in entries])
    data = sea.field(posed[:, :, :2].reshape(-1, 2), t).reshape(-1, 3, 3)
    for entry, values in zip(entries, data):
        tri = entry['posed']
        diameter2 = max(np.sum((tri[i, :2]-tri[j, :2])**2) for i in range(3) for j in range(i))
        margin = sea.roundoff_margin(tri[:, :2], t)
        epsilon = sea.C2*diameter2/2+margin
        f = tri[:, 2]-values[:, 0]
        entry['bounds'] = area_bounds(entry['mean'], tri, f, epsilon)
        departure = values[:, 0]-(plane[0]+(tri[:, :2]-plane[1:3])@plane[3:5])
        sampled_departure = float(np.max(abs(departure)))
        entry['curvature_lower'] = max(0., sampled_departure-margin)
        entry['curvature_upper'] = sampled_departure+epsilon
    return entries


def summarize(leaves, faces):
    bounds = np.zeros((faces, 4))
    for entry in leaves:
        bounds[entry['face']] += entry['bounds']
    return dict(wet_area_bounds_m2=bounds[:, :2].sum(axis=0).tolist(),
        wet_set_turnover_bounds_m2=bounds[:, 2:].sum(axis=0).tolist(),
        per_material_face_bounds_m2=bounds.tolist(),
        footprint_curvature_departure_bounds_m=[max(v['curvature_lower'] for v in leaves),
                                                max(v['curvature_upper'] for v in leaves)],
        footprint_curvature_status='ANALYTIC_SPATIAL_BOUND_WITH_REPORTED_WIDTH_NO_TARGET_ADOPTED',
        leaf_count=len(leaves), maximum_depth=max(v['depth'] for v in leaves))


def wetted_bounds(mean, posed, sea, t, plane, tolerance, maximum_depth):
    leaves = evaluate_leaves([dict(mean=a, posed=b, face=i, depth=0)
                             for i, (a, b) in enumerate(zip(mean, posed))], sea, t, plane)
    previous = None
    while True:
        result = summarize(leaves, len(mean))
        area = result['wet_area_bounds_m2']
        turn = result['wet_set_turnover_bounds_m2']
        result['area_target_met'] = max(area[1]-area[0], turn[1]-turn[0]) <= tolerance
        result['previous_iteration_bounds'] = previous
        if result['area_target_met']:
            result['status'] = 'BOUNDED_TO_DECLARED_SPATIAL_TOLERANCE'
            return result
        candidates = [v for v in leaves if v['depth'] < maximum_depth and
                      max(v['bounds'][1]-v['bounds'][0], v['bounds'][3]-v['bounds'][2]) > 1e-15]
        if not candidates:
            result['status'] = 'UNKNOWN_TARGET_NOT_MET_AT_DEPTH_CAP'
            return result
        selected = sorted(candidates, key=lambda v: v['bounds'][1]-v['bounds'][0]+v['bounds'][3]-v['bounds'][2], reverse=True)[:32]
        if len(leaves)+3*len(selected) > 8192:
            result['status'] = 'UNKNOWN_TARGET_NOT_MET_AT_ACTIVE_LEAF_MEMORY_CAP'
            return result
        children = []
        for parent in selected:
            children.extend(dict(mean=a, posed=b, face=parent['face'], depth=parent['depth']+1)
                            for a, b in zip(split(parent['mean']), split(parent['posed'])))
        try:
            evaluate_leaves(children, sea, t, plane)
        except BudgetExceeded:
            result['status'] = 'UNKNOWN_BUDGET_EXHAUSTED_WITH_VALID_COARSE_BOUNDS'
            return result
        selected_ids = {id(v) for v in selected}
        leaves = [v for v in leaves if id(v) not in selected_ids]+children
        previous = dict(wet_area_bounds_m2=area, wet_set_turnover_bounds_m2=turn)


def rim_bounds(edges, sea, t, tolerance, maximum_depth):
    leaves = []

    def evaluate(edge, depth):
        field = sea.field(edge[:, :2], t)[:, 0]
        f = edge[:, 2]-field
        margin = sea.roundoff_margin(edge[:, :2], t)
        epsilon = sea.C2*np.sum((edge[1, :2]-edge[0, :2])**2)/8+margin
        return dict(edge=edge, depth=depth, lower=float(min(f)-epsilon), upper=float(min(f)+margin))

    for edge in edges:
        leaves.append(evaluate(edge, 0))
    while True:
        lower, upper = min(v['lower'] for v in leaves), min(v['upper'] for v in leaves)
        result = dict(clearance_min_bounds_m=[lower, upper],
            opening_status='PROVEN_DRY_AT_SAVED_TIME' if lower > 0 else
                           'PROVEN_IMMERSED_AT_SAVED_TIME' if upper < 0 else 'UNKNOWN_CONTACT_OR_IMMERSION',
            numerical_target_met=upper-lower <= tolerance, segments=len(leaves))
        if result['numerical_target_met']:
            return result
        candidates = [v for v in leaves if v['lower'] < upper-tolerance and v['depth'] < maximum_depth]
        if not candidates:
            result['numerical_status'] = 'UNKNOWN_TARGET_NOT_MET_AT_DEPTH_CAP'
            return result
        parent = min(candidates, key=lambda v: v['lower'])
        a, b = parent['edge']
        middle = (a+b)/2
        try:
            children = [evaluate(np.array([a, middle]), parent['depth']+1),
                        evaluate(np.array([middle, b]), parent['depth']+1)]
        except BudgetExceeded:
            result['numerical_status'] = 'UNKNOWN_BUDGET_EXHAUSTED_WITH_VALID_COARSE_BOUNDS'
            return result
        leaves = [v for v in leaves if v is not parent]+children
