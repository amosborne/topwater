"""Differentiate a wet nodal field on one explicitly owned geometry chart."""
import numpy as np

from .body_pressure import quadrature


def differentiate(base, plus, minus, phi_plus, phi_minus, step, center, twist):
    """Return nodal ALE rate, chart velocity and facet normal-slip diagnostics.

    Chart nodes may slide tangentially on the hull. Their actual velocity, not
    the rigid body's velocity, must accompany the nodal potential derivative.
    The pressure owner alone applies the Eulerian chain-rule correction.
    The normal-slip RMS is an unweighted diagnostic over quadrature samples.
    """
    if not np.isfinite(step) or step <= 0:
        raise ValueError('positive finite pressure derivative step required')
    basis = ('faces', 'body_faces', 'free_faces', 'free_nodes', 'body_source_group',
             'free_source_parent', 'free_parent_indices')
    if 'node_chart_key' not in base:
        raise ValueError('nodal pressure derivative requires explicit chart ownership')
    count = len(base['points'])
    for moved in (plus, minus):
        if ('node_chart_key' not in moved
                or not np.array_equal(base['node_chart_key'], moved['node_chart_key'])
                or any(not np.array_equal(base[name], moved[name]) for name in basis)
                or np.shape(moved['points']) != (count, 3)):
            raise ValueError('nodal pressure derivative requires the same chart and node basis')
    for value, shape in ((base['points'], (count, 3)), (plus['points'], (count, 3)),
                         (minus['points'], (count, 3)), (phi_plus, (count,)),
                         (phi_minus, (count,)), (center, (3,)), (twist, (6,))):
        if np.shape(value) != shape or not np.isfinite(value).all():
            raise ValueError('finite nodal fields, geometry, center and twist required')
    rate = (np.asarray(phi_plus)-phi_minus)/(2*step)
    velocity = (plus['points']-minus['points'])/(2*step)
    return rate, velocity, normal_slip(base, velocity, center, twist)


def normal_slip(base, velocity, center, twist):
    """Compare current chart and rigid hull normal velocities at quadrature."""
    for value, shape in ((velocity, base['points'].shape), (center, (3,)), (twist, (6,))):
        if np.shape(value) != shape or not np.isfinite(value).all():
            raise ValueError('finite chart velocity, center and twist required')
    triangles = base['points'][base['faces'][base['body_faces']]]
    normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    bary, _ = quadrature()
    samples = np.einsum('qc,fcj->fqj', bary, triangles)
    chart_velocity = np.einsum('qc,fcj->fqj', bary,
                              velocity[base['faces'][base['body_faces']]])
    rigid_velocity = twist[:3]+np.cross(twist[3:], samples-center)
    slip = np.einsum('fqj,fj->fq', chart_velocity-rigid_velocity, normals)
    return dict(
        maximum_body_normal_velocity_error_m_per_s=float(np.max(np.abs(slip))),
        rms_body_normal_velocity_error_m_per_s=float(np.sqrt(np.mean(slip**2))))
