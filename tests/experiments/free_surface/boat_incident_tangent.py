"""Current material derivative of the diagnostic tank's incident Neumann load.

Differentiate the existing q3/P1 boundary quadrature, including its moving
sample positions and oriented face areas. No field solve or moved chart is
needed. Coefficients and the primal velocity retain their existing owner.
"""
import math

import numpy as np
from numba import njit

from tests.laplace_reference import triangle_samples
from .incident_baseline_benchmark import fused_cpu


@njit(fastmath=False)
def _material_velocity_rate(points, moving, values, time):
    result = np.zeros((len(points), 3))
    for i in range(len(points)):
        for j in range(len(values)):
            a, omega, k, dx, dy, phase0 = values[j]
            phase = k*(dx*points[i, 0]+dy*points[i, 1])-omega*time+phase0
            phase_rate = k*(dx*moving[i, 0]+dy*moving[i, 1])-omega
            sine, cosine = math.sin(phase), math.cos(phase)
            scale = a*omega*math.exp(k*points[i, 2])
            horizontal = scale*(k*moving[i, 2]*cosine-phase_rate*sine)
            result[i, 0] += dx*horizontal
            result[i, 1] += dy*horizontal
            result[i, 2] += scale*(k*moving[i, 2]*sine+phase_rate*cosine)
    return result


def incident_rate(arrays, moving, values, time):
    """Return d(b_incident)/dt for the SAME current oriented outer faces.

    ``moving`` is the current chart's node velocity. This derivative excludes
    the body's B(t)@twist term, which the action tangent differentiates itself.
    """
    points, faces = arrays['points'], arrays['faces']
    moving, values = np.asarray(moving), np.asarray(values)
    if (moving.shape != points.shape or values.ndim != 2 or values.shape[1] != 6
            or not np.isfinite(moving).all() or not np.isfinite(values).all()
            or not np.isfinite(time)):
        raise ValueError('finite current chart velocity, wave coefficients and time required')
    outer = np.r_[arrays['wall_faces'], arrays['bottom_faces']]
    triangles, velocities = points[faces[outer]], moving[faces[outer]]
    edges = triangles[:, 1:]-triangles[:, :1]
    edge_rates = velocities[:, 1:]-velocities[:, :1]
    area = .5*np.cross(edges[:, 0], edges[:, 1])
    area_rate = .5*(np.cross(edge_rates[:, 0], edges[:, 1])
                    +np.cross(edges[:, 0], edge_rates[:, 1]))
    bary, weights = triangle_samples(np.eye(3)[None], 3)
    bary = bary[0]
    samples = np.einsum('qc,fcj->fqj', bary, triangles)
    sample_rates = np.einsum('qc,fcj->fqj', bary, velocities)
    velocity = fused_cpu(samples.reshape(-1, 3), values, time)[:, 6:].reshape(samples.shape)
    rate = _material_velocity_rate(samples.reshape(-1, 3), sample_rates.reshape(-1, 3),
                                   values, time).reshape(samples.shape)
    flux_rate = (np.einsum('fqj,fj->fq', rate, area)
                 +np.einsum('fqj,fj->fq', velocity, area_rate))
    result = np.zeros(len(points))
    for corner in range(3):
        np.add.at(result, faces[outer, corner],
                  np.einsum('q,q,fq->f', weights, bary[:, corner], flux_rate))
    if not np.isfinite(result).all():
        raise FloatingPointError('nonfinite incident boundary-load derivative')
    return result
