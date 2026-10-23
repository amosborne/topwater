"""CPU tensor-query references for Fourier and structured triangular graphs.

Coordinates are (x,y), arrays (y,x). Results are (height, xy_gradient), with
gradient components on the last axis. No geometry classification is performed.
"""
import numpy as np


def _coordinates(values, name):
    result = np.asarray(values)
    if result.ndim != 1 or np.iscomplexobj(result) or not np.isfinite(result).all():
        raise ValueError(f'{name} must be a finite real coordinate vector')
    return result.astype(np.float64, copy=False)


def _fourier_data(eta, length, origin):
    eta = np.asarray(eta)
    length = np.broadcast_to(np.asarray(length, dtype=float), (2,))
    origin = np.asarray(origin, dtype=float)
    if (eta.ndim != 2 or min(eta.shape) < 3 or np.iscomplexobj(eta)
            or not np.isfinite(eta).all() or not np.isfinite(length).all()
            or np.any(length <= 0) or origin.shape != (2,) or not np.isfinite(origin).all()):
        raise ValueError('finite real graph, positive lengths and finite origin required')
    ny, nx = eta.shape
    kx = 2*np.pi*np.fft.fftfreq(nx, d=length[0]/nx)
    ky = 2*np.pi*np.fft.fftfreq(ny, d=length[1]/ny)
    spectrum = np.fft.fft2(eta.astype(np.float64))/(nx*ny)
    return kx, ky, spectrum, origin


def tensor_fourier(eta, x, y, length, origin):
    """Interpolate a periodic grid at tensor coordinates without point/mode maps.

    length is scalar or (Lx,Ly); origin locates eta[0,0]. The real Fourier
    interpolant and its derivatives are returned, including for off-grid points.
    """
    kx, ky, spectrum, origin = _fourier_data(eta, length, origin)
    x, y = _coordinates(x, 'x'), _coordinates(y, 'y')
    fx = np.exp(1j*(x-origin[0])[:, None]*kx)
    fy = np.exp(1j*(y-origin[1])[:, None]*ky)
    rows = spectrum@fx.T
    value = (fy@rows).real
    dx = (fy@(spectrum@(fx*(1j*kx)).T)).real
    dy = ((fy*(1j*ky))@rows).real
    return value, np.stack((dx, dy), axis=-1)


def fourier_points(eta, points_xy, length, origin):
    """Fourier height/xy gradient at paired points; O(Ngrid) scratch per point."""
    kx, ky, spectrum, origin = _fourier_data(eta, length, origin)
    points = np.asarray(points_xy)
    if (points.ndim != 2 or points.shape[1] != 2 or np.iscomplexobj(points)
            or not np.isfinite(points).all()):
        raise ValueError('finite real (N,2) points required')
    value, gradient = np.empty(len(points)), np.empty((len(points), 2))
    for index, point in enumerate(points):
        fx = np.exp(1j*(point[0]-origin[0])*kx)
        fy = np.exp(1j*(point[1]-origin[1])*ky)
        row = fy@spectrum
        value[index] = (row@fx).real
        gradient[index, 0] = (row@(fx*(1j*kx))).real
        gradient[index, 1] = (((fy*(1j*ky))@spectrum)@fx).real
    return value, gradient


def structured_P1(nodesx, nodesy, heights, x, y):
    """Triangular roof, each rectangular cell split bottom-left to top-right.

    Arbitrary strictly increasing axis coordinates are supported. The diagonal
    belongs to the lower-right triangle; outer endpoints belong to the last
    cell. No extrapolation or spacing-rounding is permitted.
    """
    nodesx, nodesy = _coordinates(nodesx, 'nodesx'), _coordinates(nodesy, 'nodesy')
    x, y = _coordinates(x, 'x'), _coordinates(y, 'y')
    heights = np.asarray(heights)
    if (min(len(nodesx), len(nodesy)) < 2 or np.any(np.diff(nodesx) <= 0)
            or np.any(np.diff(nodesy) <= 0) or heights.shape != (len(nodesy), len(nodesx))
            or np.iscomplexobj(heights) or not np.isfinite(heights).all()):
        raise ValueError('increasing axes and matching finite real nodal heights required')
    if (np.any(x < nodesx[0]) or np.any(x > nodesx[-1])
            or np.any(y < nodesy[0]) or np.any(y > nodesy[-1])):
        raise ValueError('query lies outside the parent roof')
    ix = np.minimum(np.searchsorted(nodesx, x, side='right')-1, len(nodesx)-2)
    iy = np.minimum(np.searchsorted(nodesy, y, side='right')-1, len(nodesy)-2)
    width, depth = np.diff(nodesx)[ix], np.diff(nodesy)[iy]
    u = ((x-nodesx[ix])/width)[None, :]
    v = ((y-nodesy[iy])/depth)[:, None]
    h00 = heights[iy[:, None], ix]
    h10 = heights[iy[:, None], ix+1]
    h01 = heights[iy[:, None]+1, ix]
    h11 = heights[iy[:, None]+1, ix+1]
    lower = v <= u
    du = np.where(lower, h10-h00, h11-h01)
    dv = np.where(lower, h11-h10, h01-h00)
    value = h00+du*u+dv*v
    return value, np.stack((du/width, dv/depth[:, None]), axis=-1)
