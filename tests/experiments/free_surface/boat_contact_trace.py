"""Compare two discrete normal traces on one unchanged contact geometry."""
import numpy as np

from tests.laplace_reference import triangle_samples


def compare_traces(points, faces, free_nodes, potential, recovered_density,
                   auxiliary_normal, contact_edges, gravity=9.81):
    """No trace is declared exact; compare signed values on the SAME samples."""
    local = np.full(len(points), -1, dtype=int)
    local[free_nodes] = np.arange(len(free_nodes))
    triangles = points[faces]
    delta = triangles[:, 1:, :2]-triangles[:, :1, :2]
    area = np.linalg.det(delta)/2
    slope = np.linalg.solve(delta, (triangles[:, 1:, 2]-triangles[:, :1, 2])[..., None])[..., 0]
    tangent = np.linalg.solve(delta, (potential[faces[:, 1:]]-potential[faces[:, :1]])[..., None])[..., 0]
    jacobian_square = 1+np.sum(slope*slope, axis=1)
    shape, weights = triangle_samples(np.eye(3)[None], 3)
    shape = shape[0]
    locations = np.einsum('qc,fck->fqk', shape, triangles)
    recovered = np.einsum('qc,fc->fq', shape, recovered_density[local[faces]])
    auxiliary = np.broadcast_to((auxiliary_normal*np.sqrt(jacobian_square))[:, None], recovered.shape)
    squared_distance = np.full(recovered.shape, np.inf)
    for a, b in points[contact_edges]:
        direction = b-a
        parameter = np.clip(np.sum((locations-a)*direction, axis=2)/(direction@direction), 0, 1)
        squared_distance = np.minimum(squared_distance,
            np.sum((locations-a-parameter[..., None]*direction)**2, axis=2))
    distance = np.sqrt(squared_distance)
    weighted_area = area[:, None]*weights
    tangent_dot_slope = np.sum(slope*tangent, axis=1)
    z = locations[:, :, 2]
    verticals, bernoullis = {}, {}
    for name, flux in (('reaction', recovered), ('auxiliary', auxiliary)):
        verticals[name] = (flux+tangent_dot_slope[:, None])/jacobian_square[:, None]
        bernoullis[name] = (-gravity*z-.5*np.sum(tangent*tangent, axis=1)[:, None]
                            +.5*jacobian_square[:, None]*verticals[name]**2)
    peak = np.unravel_index(np.argmax(abs(verticals['reaction'])), recovered.shape)
    report = dict(
        scope='Discrete reaction P1 recovery versus P0 auxiliary trace; neither is independently validated physical velocity.',
        sampling='Nine-point positive Gauss-Duffy rule per face, shared with surface_rates. Contact-bin integrals use sampled indicators, not geometric strip clipping.',
        auxiliary_vertical_at_reaction_peak_m_per_s=float(verticals['auxiliary'][peak]),
        reaction_vertical_at_reaction_peak_m_per_s=float(verticals['reaction'][peak]),
        bins=[])
    bounds = [0., 1e-4, 1e-3, 1e-2, .1, 1., np.inf]
    for lower, upper in zip(bounds[:-1], bounds[1:]):
        selected = (distance >= lower) & (distance < upper)
        row = dict(contact_distance_lower_m=lower,
                   contact_distance_upper_m=upper if np.isfinite(upper) else None,
                   quadrature_projected_area_m2=float(np.sum(weighted_area*selected)))
        for name, flux in (('reaction', recovered), ('auxiliary', auxiliary)):
            row[name] = dict(
                graph_flux_integral_m3_per_s=float(np.sum(weighted_area*selected*flux)),
                vertical_square_integral_m4_per_s2=float(np.sum(weighted_area*selected*verticals[name]**2)),
                bernoulli_integral_m4_per_s2=float(np.sum(weighted_area*selected*bernoullis[name])),
                maximum_absolute_vertical_m_per_s=float(np.max(abs(verticals[name][selected]))) if selected.any() else None)
        report['bins'].append(row)
    return report
