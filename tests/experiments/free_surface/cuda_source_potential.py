"""Guarded GPU paired source potential with bounded CPU precision recovery."""
from dataclasses import dataclass
from decimal import DecimalException
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss

from topwater.numerics import LaplaceControls
from topwater.rankine import TrianglePanels, panel_influence, QuadratureError
from topwater._rankine_extended import triangle_values
from .cuda_operator import cp
from .panel_gradient_reference import extended

_source = (Path(__file__).with_name('cuda_pair.cuh').read_text()+'\n'
           +Path(__file__).with_suffix('.cu').read_text())
_kernel = cp.RawKernel(_source, 'paired_source_potential', options=('--std=c++11', '--fmad=false'))
_gradient_kernel = cp.RawKernel(_source, 'paired_source_gradient', options=('--std=c++11', '--fmad=false'))


@dataclass(frozen=True)
class SourcePotential:
    values: object
    diagnostics: dict


@dataclass(frozen=True)
class SourceGradient:
    values: object
    diagnostics: dict


def _far_pair_gradient(point, triangle, image, controls):
    """Order-checked stable paired integrand; both points must be well separated."""
    diameter = np.linalg.norm(np.roll(triangle,-1,axis=0)-triangle,axis=1).max()
    reflected = point.copy(); reflected[2] = 2*image-point[2]
    center = triangle.mean(axis=0)
    if min(np.linalg.norm(point-center),np.linalg.norm(reflected-center)) < 2*diameter:
        return None
    area2=np.linalg.norm(np.cross(triangle[1]-triangle[0],triangle[2]-triangle[0]))
    previous=None
    for refinement in range(controls.max_refinement+1):
        nodes,weights=leggauss(4*(refinement+1))
        u,v=np.meshgrid((nodes+1)/2,(nodes+1)/2,indexing='ij')
        y=triangle[0]+u[:,:,None]*((1-v[:,:,None])*(triangle[1]-triangle[0])
                                   +v[:,:,None]*(triangle[2]-triangle[0]))
        delta=point-y; mirror=reflected-y
        radius=np.hypot(np.hypot(delta[:,:,0],delta[:,:,1]),delta[:,:,2])
        image_radius=np.hypot(np.hypot(mirror[:,:,0],mirror[:,:,1]),mirror[:,:,2])
        ir,im=1/radius,1/image_radius
        difference=4*(image-point[2])*(image-y[:,:,2])/(radius+image_radius)
        inverse_cube_difference=difference*ir*im*(ir*ir+ir*im+im*im)
        gradient=-delta*inverse_cube_difference[:,:,None]
        absolute_terms=abs(gradient)
        image_vertical=2*(image-y[:,:,2])*im**3
        gradient[:,:,2]-=image_vertical
        absolute_terms[:,:,2]+=abs(image_vertical)
        measure=np.outer(weights,weights)*u*area2/(16*np.pi)
        terms=measure[:,:,None]*gradient
        value=np.sum(terms,axis=(0,1))
        error=64*np.finfo(float).eps*np.sum(measure[:,:,None]*absolute_terms,axis=(0,1))
        if (previous is not None and np.isfinite(value).all() and np.isfinite(error).all()
                and np.all(abs(value-previous)+error <= controls.quadrature_atol+controls.quadrature_rtol*abs(value))):
            return value
        previous=value
    raise QuadratureError('far paired gradient exceeded requested accuracy/work budget')


def _recover_gradient_pair(point, panels, column, image, controls):
    far=_far_pair_gradient(point,panels.triangles[column],image,controls)
    if far is not None:
        return far
    reflected=point.copy(); reflected[2]=2*image-point[2]
    previous=None
    for refinement in range(min(controls.max_refinement+1,5)):
        try:
            real,er=extended(point,panels.triangles[column],panels.normal_sign[column],32+16*refinement)
            mirror,em=extended(reflected,panels.triangles[column],panels.normal_sign[column],32+16*refinement)
        except (DecimalException,ValueError,ZeroDivisionError):
            previous=None
            continue
        value=real[2:5]-mirror[2:5]*[1.,1.,-1.]
        error=er[2:5]+em[2:5]+4*np.finfo(float).eps*(abs(real[2:5])+abs(mirror[2:5])+abs(value))
        if np.array_equal(point,reflected):
            value[:2]=0.;error[:2]=0.
        if (previous is not None and np.isfinite(value).all() and np.isfinite(error).all()
                and np.all(abs(value-previous)+error <= controls.quadrature_atol+controls.quadrature_rtol*abs(value))):
            return value
        previous=value
    raise QuadratureError('paired gradient cancellation exceeds requested binary64 accuracy')


def _recover_pair(point, panels, column, image, controls):
    reflected = point.copy(); reflected[2] = 2*image-point[2]
    # Rejection guards have already run through public panel_influence.
    if np.array_equal(point, reflected):
        return 0.
    previous = None
    for refinement in range(min(controls.max_refinement+1, 5)):
        try:
            real, er = triangle_values(point, panels.triangles[column], panels.normals[column], 32+16*refinement)
            mirror, em = triangle_values(reflected, panels.triangles[column], panels.normals[column], 32+16*refinement)
        except (DecimalException, ValueError, ZeroDivisionError):
            previous = None
            continue
        if not all(np.isfinite(v) and v > 0 for v in [real[0], mirror[0]]):
            previous = None
            continue
        value = real[0]-mirror[0]
        error = er[0]+em[0]+4*np.finfo(float).eps*(abs(real[0])+abs(mirror[0])+abs(value))
        tolerance = controls.quadrature_atol*panels.diameters[column]+controls.quadrature_rtol*abs(value)
        if previous is not None and np.isfinite(value) and abs(value-previous)+error <= tolerance:
            return value
        previous = value
    raise QuadratureError('paired source cancellation exceeds requested binary64 accuracy')


def _inputs(points, panels, image, controls):
    controls = LaplaceControls() if controls is None else controls
    points = np.asarray(points)
    image = np.asarray(image)
    if (not isinstance(panels, TrianglePanels) or not isinstance(controls, LaplaceControls)
            or np.iscomplexobj(points) or points.ndim != 2 or points.shape[1:] != (3,)
            or not np.isfinite(points).all() or image.shape != () or image.dtype.kind not in 'fiu'
            or not np.isfinite(image)):
        raise ValueError('finite bulk points, triangle panels, controls and image height required')
    points = np.ascontiguousarray(points, dtype=np.float64)
    image = float(image)
    return points,panels,image,controls


def _evaluate(points, panels, image, controls, gradients):
    points,panels,image,controls=_inputs(points,panels,image,controls)
    count = len(points)*len(panels)
    components=3 if gradients else 1
    if count*(8*components+1) > controls.max_operator_bytes:
        raise MemoryError('paired source matrix and guard flags exceed the operator memory budget')
    shape=(len(points),len(panels),3) if gradients else (len(points),len(panels))
    values = cp.empty(shape, dtype=cp.float64)
    flags = cp.empty(count, dtype=cp.uint8)
    if count:
        geometry = [cp.asarray(np.ascontiguousarray(a)) for a in
            [panels.triangles, panels.geometric_normals, panels.edge_outward, panels.tangents,
             panels.edge_lengths, panels.areas, panels.centroids, panels.normal_sign, panels.diameters]]
        kernel=_gradient_kernel if gradients else _kernel
        kernel(((count+127)//128,), (128,), (*geometry, cp.asarray(points), np.int32(len(points)),
            np.int32(len(panels)), np.float64(image), np.float64(controls.quadrature_atol),
            np.float64(controls.quadrature_rtol), values, flags))
    indices = cp.asnumpy(cp.flatnonzero(flags))
    status = cp.asnumpy(flags[cp.asarray(indices)])
    rows = np.unique(indices//len(panels))
    # Recompute flagged rows with unchanged public geometry/accuracy guards.
    # No geometry rejection is converted into a guessed value.
    for first in range(0, len(rows), controls.target_block_size):
        selected = rows[first:first+controls.target_block_size]
        reflected = points[selected].copy(); reflected[:, 2] = 2*image-reflected[:, 2]
        for target in (points[selected],reflected):
            try:
                panel_influence(target,panels,controls=controls,gradients=gradients)
            except QuadratureError:
                if not gradients:
                    raise
                # Public geometry guards run before numerical fallback. Paired
                # gradient precision is qualified independently below; unused
                # double-layer gradients need not be returned by this API.
    for index, bits in zip(indices, status, strict=True):
        if bits & 6:
            row, column = divmod(int(index), len(panels))
            recover=_recover_gradient_pair if gradients else _recover_pair
            recovered=recover(points[row],panels,column,image,controls)
            values[row,column]=cp.asarray(recovered) if gradients else recovered
    result=SourceGradient if gradients else SourcePotential
    return result(values, dict(gpu_flagged_pairs=len(indices), cpu_guarded_rows=len(rows),
        individual_accuracy_flags=int(np.count_nonzero(status & 2)),
        paired_accuracy_flags=int(np.count_nonzero(status & 4)), geometry_guard_flags=int(np.count_nonzero(status & 1)),
        matrix_bytes=int(values.nbytes), guard_bytes=int(flags.nbytes),
        scope='Same analytic sample primitive; explicit nonself points; CPU guards/recovery for flagged pairs'))


def paired_potential(points, panels, image, controls=None):
    return _evaluate(points,panels,image,controls,False)


def paired_gradient(points, panels, image, controls=None):
    """Bulk gradient matrix; no implicit trace limits or owned-body jumps."""
    return _evaluate(points,panels,image,controls,True)


def paired_velocity(points, panels, image, density, controls=None):
    """Contract bulk gradients with density in memory-bounded target batches."""
    points,panels,image,controls=_inputs(points,panels,image,controls)
    if len(points) and 25*len(panels)>controls.max_operator_bytes:
        raise MemoryError('one paired gradient row exceeds the operator memory budget')
    density=cp.asarray(density)
    if density.shape != (len(panels),) or density.dtype.kind not in 'fiu' or not bool(cp.isfinite(density).all()):
        raise ValueError('one finite real source density per panel required')
    density=density.astype(cp.float64,copy=False)
    velocity=cp.empty((len(points),3),dtype=cp.float64)
    batch=max(1,controls.max_operator_bytes//(25*len(panels)))
    diagnostics=[]
    for first in range(0,len(points),batch):
        last=min(first+batch,len(points))
        result=paired_gradient(points[first:last],panels,image,controls)
        velocity[first:last]=cp.einsum('tnj,n->tj',result.values,density)
        diagnostics.append(result.diagnostics)
        del result
    return SourceGradient(velocity,dict(batches=len(diagnostics),
        gpu_flagged_pairs=sum(d['gpu_flagged_pairs'] for d in diagnostics),
        cpu_guarded_rows=sum(d['cpu_guarded_rows'] for d in diagnostics),
        max_matrix_bytes=max((d['matrix_bytes'] for d in diagnostics),default=0),
        scope='Contracted bulk velocity; matrix/guard storage bounded per batch, output vector separate'))
