"""Offset field actions for a fixed density, without an offset body operator."""
from time import perf_counter
import numpy as np

from topwater.rankine import single_layer_trace
from .boat_image_plane import elevated_potential
from .panel_gradient_reference import checked_influence
from .cuda_graph import GraphDirichlet, cp


_harmonic_action = cp.RawKernel(r'''
extern "C" __global__ void harmonic_action(
    const double* points, const double* normals, const double* origin,
    const double2* spectrum, const int size, const double length, double* result) {
    const int target = blockIdx.x, lane = threadIdx.x;
    const int modes = size*size;
    double value = 0., normal = 0.;
    const double x = points[3*target]-origin[0];
    const double y = points[3*target+1]-origin[1];
    const double z = points[3*target+2];
    for (int mode=lane; mode<modes; mode+=blockDim.x) {
        int ix=mode%size, iy=mode/size;
        if (ix>size/2) ix-=size;
        if (iy>size/2) iy-=size;
        const double kx=6.2831853071795864769*ix/length;
        const double ky=6.2831853071795864769*iy/length;
        const double wave=hypot(kx,ky);
        double s,c;
        sincos(kx*x+ky*y,&s,&c);
        const double amplitude=exp(wave*z)/modes;
        const double real=amplitude*(spectrum[mode].x*c-spectrum[mode].y*s);
        const double imag=amplitude*(spectrum[mode].x*s+spectrum[mode].y*c);
        value+=real;
        normal+=normals[3*target+2]*wave*real
               -(normals[3*target]*kx+normals[3*target+1]*ky)*imag;
    }
    __shared__ double values[256], derivatives[256];
    values[lane]=value; derivatives[lane]=normal;
    __syncthreads();
    for (int stride=128; stride>0; stride/=2) {
        if (lane<stride) {
            values[lane]+=values[lane+stride];
            derivatives[lane]+=derivatives[lane+stride];
        }
        __syncthreads();
    }
    if (lane==0) { result[2*target]=values[0]; result[2*target+1]=derivatives[0]; }
}
''', 'harmonic_action')


def harmonic_action(points, normals, spectrum, size, length, origin):
    """One harmonic spectrum evaluated directly; no target-by-mode maps."""
    points, normals, origin = (np.ascontiguousarray(a, dtype=np.float64)
                               for a in (points, normals, origin))
    spectrum = cp.ascontiguousarray(spectrum, dtype=cp.complex128)
    if (points.ndim != 2 or points.shape[1:] != (3,) or normals.shape != points.shape
            or origin.shape != (2,) or size < 3 or size % 2 != 1
            or spectrum.shape != (size*size,) or length <= 0 or not np.isfinite(length)
            or not all(np.isfinite(a).all() for a in (points, normals, origin))
            or np.any(points[:, 2] > 0)):
        raise ValueError('submerged targets, matching normals and one odd-grid harmonic spectrum required')
    result = cp.empty((len(points), 2), dtype=cp.float64)
    _harmonic_action((len(points),), (256,),
        (cp.asarray(points), cp.asarray(normals), cp.asarray(origin), spectrum,
         np.int32(size), np.float64(length), result))
    return cp.asnumpy(result)


def frozen_field(panels, eta, length, origin, height, image, order, trace, density):
    """Return body phi and dn(phi) for the fixed density on current panels.

    Only scalar source values are sampled on the graph. Source gradients are
    needed at owned body centroids only. No J/S operator, factor or dense
    harmonic continuation maps are constructed.
    """
    started = perf_counter()
    grid = eta.shape[0]
    x, y = np.meshgrid(origin[0]+np.arange(grid)*length/grid,
                       origin[1]+np.arange(grid)*length/grid)
    points = np.column_stack((x.ravel(), y.ravel(), eta.ravel()))
    source_trace = np.empty(len(points))
    for first in range(0, len(points), 256):
        last = min(first+256, len(points))
        source_trace[first:last] = elevated_potential(points[first:last], panels, image)@density
    trace_seconds = perf_counter()-started
    graph = GraphDirichlet(eta-height, [length, length], order=order)
    spectrum = sum(graph.coefficients((cp.asarray(trace)-cp.asarray(source_trace)).reshape(eta.shape))).ravel()
    shifted = panels.centroids-[0., 0., height]
    harmonic = harmonic_action(shifted, panels.normals, spectrum, graph.padded_shape[0], length, origin)
    potential, normal = harmonic[:, 0], harmonic[:, 1]
    for first in range(0, len(panels), 128):
        last = min(first+128, len(panels))
        real = single_layer_trace(panels, np.arange(first, last))
        reflected = panels.centroids[first:last].copy()
        reflected[:, 2] = 2*image-reflected[:, 2]
        reflected_field = checked_influence(reflected, panels)
        potential[first:last] += (real.single_layer-reflected_field.single_layer)@density
        gradient = real.gradient_single_layer-reflected_field.gradient_single_layer*np.array([1., 1., -1.])
        normal[first:last] += np.einsum('tni,n,ti->t', gradient, density, panels.normals[first:last])
    if not np.isfinite(potential).all() or not np.isfinite(normal).all():
        raise ValueError('nonfinite frozen field action')
    return dict(potential=potential, normal=normal,
        diagnostics=dict(total_seconds=perf_counter()-started, source_trace_seconds=trace_seconds,
                         offset_operator_builds=0, dense_harmonic_maps=0))
