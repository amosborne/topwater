"""Independent positive-depth Poisson-work check, not a Galerkin V assembly.

Full-plane compact P1 traces only: no hull-hole correction or time evolution.
Quadpack estimates and layer accuracy requests are NOT certified error bounds.
"""
from time import perf_counter
import numpy as np
from scipy.integrate import quad
from .compact_p1_trace import CompactP1Trace


def h1_seminorm_squared(trace):
    """Exact triangle integration of |grad f|² for the declared P1 trace."""
    nodal = trace.values.ravel()[trace.faces]
    edges = trace.panels.triangles[:, 1:, :2]-trace.panels.triangles[:, :1, :2]
    gradient = np.linalg.solve(edges, (nodal[:, 1:]-nodal[:, :1])[..., None])[..., 0]
    return float(np.sum(trace.panels.areas*np.sum(gradient**2, axis=1)))


def poisson_work(trace, depth, *, absolute_tolerance, max_evaluations):
    """Integrate E_h=int f(x)*u_z[f](x,-h) dx with adaptive target quadrature.

    Parseval gives 0 <= E_0-E_h <= h*int|grad f|² in the CONTINUUM. The
    returned interval inserts numerical E_h, whose integration/layer errors
    remain separate. No extrapolation, positivity repair, or mean projection.
    """
    if (not isinstance(trace, CompactP1Trace) or not np.isfinite(depth) or depth <= 0
            or not np.isfinite(absolute_tolerance) or absolute_tolerance <= 0
            or isinstance(max_evaluations, bool) or not isinstance(max_evaluations, int)
            or max_evaluations <= 0):
        raise ValueError('compact trace, positive depth/tolerance and integer evaluation budget required')
    started = perf_counter()
    nodal = trace.values.ravel()[trace.faces]
    active = np.flatnonzero(np.any(nodal != 0, axis=1))
    per_face = absolute_tolerance/max(1, len(active))
    evaluations = 0
    work, estimate = 0., 0.
    for face in active:
        triangle, values = trace.panels.triangles[face], nodal[face]
        edge1, edge2 = triangle[1]-triangle[0], triangle[2]-triangle[0]
        inner_max = 0.

        def outer(radial):
            nonlocal inner_max

            def integrand(transverse):
                nonlocal evaluations
                if evaluations >= max_evaluations:
                    raise RuntimeError(f'Poisson work target budget exhausted: {evaluations}/{max_evaluations}')
                point = triangle[0]+radial*((1-transverse)*edge1+transverse*edge2)
                value = (1-radial)*values[0]+radial*((1-transverse)*values[1]+transverse*values[2])
                evaluations += 1
                velocity = trace.bulk(np.array([[point[0], point[1], -depth]]))[1][0, 2]
                return 2*trace.panels.areas[face]*radial*value*velocity

            answer = quad(integrand, 0., 1., epsabs=per_face/4, epsrel=0., limit=32, full_output=1)
            if len(answer) != 3:
                raise RuntimeError(f'Poisson work inner quadrature did not converge: {answer[3]}')
            inner_max = max(inner_max, float(answer[1]))
            return answer[0]

        answer = quad(outer, 0., 1., epsabs=per_face/2, epsrel=0., limit=32, full_output=1)
        if len(answer) != 3:
            raise RuntimeError(f'Poisson work outer quadrature did not converge: {answer[3]}')
        work += answer[0]
        estimate += answer[1]+inner_max
    seminorm = h1_seminorm_squared(trace)
    width = depth*seminorm
    if not np.isfinite([work, estimate, seminorm]).all():
        raise ValueError('nonfinite Poisson-work result')
    return dict(depth_m=float(depth), work=float(work), h1_seminorm_squared=seminorm,
        continuum_error_width=width, inserted_continuum_interval=[float(work), float(work+width)],
        quadrature_error_estimate=float(estimate),
        estimate_expanded_interval=[float(work-estimate), float(work+width+estimate)],
        target_evaluations=evaluations, active_triangles=len(active),
        seconds=perf_counter()-started, certified_total_error=False,
        error_scope='Continuum Parseval width is exact; inserted work uses adaptive target quadrature estimates and existing Rankine accuracy requests, not interval arithmetic.')


def hat_pair(hx, hy, offset):
    """Two uniform-mesh interior hats; return their sum and difference traces."""
    offset = np.asarray(offset)
    if (offset.shape != (2,) or offset.dtype.kind not in 'iu'
            or not np.isfinite([hx, hy]).all() or min(hx, hy) <= 0):
        raise ValueError('positive spacings and integer two-dimensional offset required')
    ix = np.arange(min(0, offset[0])-1, max(0, offset[0])+2)
    iy = np.arange(min(0, offset[1])-1, max(0, offset[1])+2)
    plus, minus = np.zeros((len(iy), len(ix))), np.zeros((len(iy), len(ix)))
    for point, sign in [(np.array([0, 0]), 1), (offset, -1)]:
        row, col = int(np.flatnonzero(iy == point[1])[0]), int(np.flatnonzero(ix == point[0])[0])
        plus[row, col] += 1
        minus[row, col] += sign
    return CompactP1Trace(ix*hx, iy*hy, plus), CompactP1Trace(ix*hx, iy*hy, minus)


def polarized_entry(hx, hy, offset, depth, *, absolute_tolerance, max_evaluations):
    """Independent W_ij interval from E(f_i+f_j)-E(f_i-f_j), divided by four.

    Both integrals share one finite target-evaluation budget. No translation
    kernel or Galerkin triangle-pair V quadrature is called.
    """
    plus, minus = hat_pair(hx, hy, offset)
    a = poisson_work(plus, depth, absolute_tolerance=2*absolute_tolerance,
                     max_evaluations=max_evaluations)
    remaining = max_evaluations-a['target_evaluations']
    if remaining <= 0:
        raise RuntimeError('Poisson work exhausted shared budget before negative polarization')
    b = poisson_work(minus, depth, absolute_tolerance=2*absolute_tolerance,
                     max_evaluations=remaining)
    value = (a['work']-b['work'])/4
    lower = value-b['continuum_error_width']/4
    upper = value+a['continuum_error_width']/4
    estimate = (a['quadrature_error_estimate']+b['quadrature_error_estimate'])/4
    return dict(offset=np.asarray(offset).tolist(), depth_m=float(depth), value=value,
        inserted_continuum_interval=[lower, upper], quadrature_error_estimate=estimate,
        estimate_expanded_interval=[lower-estimate, upper+estimate],
        target_evaluations=a['target_evaluations']+b['target_evaluations'],
        plus=a, minus=b, certified_total_error=False)
