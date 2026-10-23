"""Local image-height bound for the auxiliary periodic Fourier graph.

Only the graph above the projected source support constrains its images.
This helper chooses a sufficient height; it does not compute a least height,
change the physical free surface, or initialize a moving solver.
"""
import numpy as np


def image_height(triangles, eta, lengths, origin):
    """Bound every reflected source point above the represented graph.

    A trigonometric interpolant obeys |eta_x| <= sum|kx*eta_hat|,
    and similarly for y. Nearest grid nodes are at most half a spacing away
    in each coordinate. Local node maxima plus those bounds cover the entire
    source projection, including triangle interiors. The floating allowance
    is reported explicitly; this is not an interval-arithmetic certificate.
    """
    triangles, eta = np.asarray(triangles, dtype=float), np.asarray(eta, dtype=float)
    lengths, origin = np.asarray(lengths, dtype=float), np.asarray(origin, dtype=float)
    if (triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not len(triangles)
            or eta.ndim != 2 or min(eta.shape) < 3 or any(n % 2 != 1 for n in eta.shape)
            or lengths.shape != (2,) or origin.shape != (2,) or np.any(lengths <= 0)
            or not all(np.isfinite(x).all() for x in (triangles, eta, lengths, origin))):
        raise ValueError('finite source triangles, odd Fourier graph and positive XY box lengths required')
    ny, nx = eta.shape
    spacing = lengths/np.array([nx, ny])
    xy = triangles[..., :2].reshape(-1, 2)
    low = np.floor((xy.min(axis=0)-origin)/spacing).astype(int)
    high = np.ceil((xy.max(axis=0)-origin)/spacing).astype(int)
    if np.any(low < 0) or np.any(high >= [nx, ny]):
        raise ValueError('projected source support must fit inside the grid-node rectangle')
    spectrum = np.fft.fft2(eta)/eta.size
    kx = 2*np.pi*np.fft.fftfreq(nx, d=spacing[0])[None, :]
    ky = 2*np.pi*np.fft.fftfreq(ny, d=spacing[1])[:, None]
    derivative = np.array([np.sum(np.abs(kx*spectrum)), np.sum(np.abs(ky*spectrum))])
    local_max = float(eta[low[1]:high[1]+1, low[0]:high[0]+1].max())
    source_max = float(triangles[..., 2].max())
    scale = max(1., float(np.sum(np.abs(spectrum))), float(np.max(np.abs(triangles[..., 2]))))
    allowance = 64*np.finfo(float).eps*eta.size*scale
    graph_upper = local_max+.5*float(derivative@spacing)+allowance
    height = float(np.nextafter(max(source_max, .5*(source_max+graph_upper)), np.inf))
    return height, dict(image_plane_height_m=height,
        source_maximum_height_m=source_max, local_graph_node_maximum_m=local_max,
        global_graph_node_maximum_m=float(eta.max()), graph_upper_over_source_projection_m=graph_upper,
        derivative_upper_bounds=derivative.tolist(), rounding_allowance_m=allowance,
        source_rectangle_node_low_xy=low.tolist(), source_rectangle_node_high_xy=high.tolist(),
        minimum_image_graph_gap_bound_m=2*height-source_max-graph_upper,
        scope='Analytic derivative bound for represented Fourier graph with explicit floating allowance; not interval arithmetic')
