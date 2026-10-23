"""Current paired-source velocity and owned body trace; no time differentiation."""
import numpy as np

from topwater.rankine import single_layer_trace
from .cuda_graph import harmonic_velocity, cp
from .panel_gradient_reference import checked_influence
from .boat_radiation_time import harmonic_maps
from .boat_image_plane import elevated_potential
from .cuda_source_potential import paired_velocity


def _actual_harmonic_values(context, spectrum, points):
    targets = points.copy()
    targets[:, 2] -= context.reference_height
    size = context.graph.padded_shape[0]
    k = 2*np.pi*np.fft.fftfreq(size, d=context.length/size)
    maximum_exponent = float(max(0., targets[:, 2].max())*np.sqrt(2)*abs(k).max())
    if maximum_exponent >= np.log(np.finfo(float).max):
        raise ValueError('actual graph trace diagnostic would overflow its finite harmonic basis')
    potential, _, _ = harmonic_maps(targets, np.zeros_like(targets), size,
                                    context.length/size, context.origin)
    value = (potential@cp.asnumpy(spectrum)).real
    if not np.isfinite(value).all():
        raise ValueError('actual graph trace diagnostic is not finite')
    return value, maximum_exponent


def actual_graph_values(context, density, spectrum, indices):
    """True finite field at owned graph nodes, not formal surface Taylor values."""
    value, exponent = _actual_harmonic_values(context, spectrum, context.surface_points[indices])
    return value+cp.asnumpy(context.trace[indices]@cp.asarray(density)), exponent


def actual_field_values(context, density, spectrum, points):
    """Diagnostic actual finite field; positive-height exponential is bounded."""
    value, exponent = _actual_harmonic_values(context, spectrum, points)
    source = elevated_potential(points, context.panels, context.image_height)@density
    return value+source, exponent


def _image_gradient(points, panels, height):
    reflected = points.copy()
    reflected[:, 2] = 2*height-reflected[:, 2]
    return checked_influence(reflected, panels).gradient_single_layer*np.array([1., 1., -1.])


def source_velocity_bulk(context, density, points):
    return cp.asnumpy(paired_velocity(points,context.panels,context.image_height,density).values)


def source_velocity_body(context, density):
    panels = context.panels
    velocity = np.empty((len(panels), 3))
    for first in range(0, len(panels), 128):
        last = min(first+128, len(panels))
        real = single_layer_trace(panels, np.arange(first, last)).gradient_single_layer
        gradient = real-_image_gradient(panels.centroids[first:last], panels, context.image_height)
        velocity[first:last] = np.einsum('tnj,n->tj', gradient, density)
    return velocity


def body_velocity(context, density, spectrum, length):
    lift = cp.asnumpy(harmonic_velocity(context.potential, spectrum,
                                        context.graph.padded_shape[0], length))
    return lift+source_velocity_body(context, density)
