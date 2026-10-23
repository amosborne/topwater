"""One bounded experimental raw-pair Taylor anchor; no cache/refresh policy.

All potential dependence is current. The approximate energy and its derivative
share cached raw V/H slopes, including the self jump, and exact current curls.
Caller-owned chart tokens must identify node semantics, not merely face count.
The caller must first run the ordinary physical/chart guards. This module does
not certify approximation accuracy, and never repairs a failed factorization.
"""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from .cuda_operator import cp, _curl
from cupyx import errstate
from cupyx.scipy.linalg import solve_triangular
from .cuda_shape_gradient import _geometry_jacobians

_folder=Path(__file__).parent
# Use the existing geometry view and numerical pair body, without maintaining
# a second copy of their formulas. Only the extra cache kernels live here.
_prefix=(_folder/'cuda_shape.cu').read_text().split('extern "C" __global__ void vertex_gradient(')[0]
for _name in ('cuda_dual.cuh','cuda_wide_dual.cuh','cuda_pair.cuh'):
    _prefix=_prefix.replace(f'#include "{_name}"',(_folder/_name).read_text())
_source=_prefix+(_folder/'cuda_kernel_linearization.cu').read_text()
_options=('--std=c++11','-DSHAPE_WIDTH=3')
_build=cp.RawKernel(_source,'anchor_pairs',options=_options)
_update=cp.RawKernel(_source,'linear_pairs',options=_options)
_pair_gradient=cp.RawKernel(_source,'cached_pair_gradient',options=_options)
_curl_gradient=cp.RawKernel(_source,'current_curl_gradient',options=_options)


def _validate(points,faces):
    points=np.asarray(points,dtype=float);faces=np.asarray(faces)
    if (points.ndim!=2 or points.shape[1]!=3 or not len(points)
            or not np.all(np.isfinite(points)) or faces.ndim!=2 or faces.shape[1]!=3
            or not len(faces) or faces.dtype.kind not in 'iu' or np.any(faces<0)
            or np.any(faces>=len(points)) or len(np.unique(faces))!=len(points)):
        raise ValueError('finite points and contiguous fully used triangular node numbering required')
    return points,faces


@dataclass(frozen=True,eq=False)
class LinearizedAssembly:
    """Current geometry, approximate energy and its matching factor/trace."""
    energy: np.ndarray
    diagnostics: dict
    factor: object
    mixed_trace: object
    points: np.ndarray
    single_layer: object
    curl: np.ndarray
    curl_jacobian: np.ndarray
    anchor: object

    def auxiliary(self,points,faces,potential,*,order=3):
        potential=np.asarray(potential,dtype=float)
        if (order!=self.anchor.order or not np.array_equal(points,self.points)
                or not np.array_equal(faces,self.anchor.faces)
                or potential.shape!=(len(self.points),) or not np.all(np.isfinite(potential))):
            raise ValueError('auxiliary reuse requires exact current geometry/order and finite potential')
        with errstate(linalg='raise'):
            q=solve_triangular(self.factor,self.mixed_trace@cp.asarray(potential),lower=True)
            return solve_triangular(self.factor.T,q,lower=False)

    def vertex_gradient(self,potential,*,density=1000.):
        """Derivative of this SAME approximate energy at fixed nodal phi."""
        density=float(density)
        if not np.isfinite(density) or density<=0:
            raise ValueError('finite positive density required')
        phi=np.asarray(potential,dtype=float)
        q=self.auxiliary(self.points,self.anchor.faces,phi,order=self.anchor.order)
        local_phi=phi[self.anchor.faces]
        c=cp.asarray(np.einsum('fci,fc->fi',self.curl,local_phi))
        dc=cp.asarray(np.ascontiguousarray(np.einsum('fsci,fc->fsi',self.curl_jacobian,local_phi)))
        out=cp.zeros_like(cp.asarray(self.points))
        nf=len(self.anchor.faces)
        _pair_gradient(((nf*nf+127)//128,),(128,),
            (self.anchor.jacobian,self.anchor.device_faces,c,cp.asarray(phi),q,
             np.int32(nf),np.float64(density),out))
        _curl_gradient(((nf*9+127)//128,),(128,),
            (self.anchor.device_faces,dc,self.single_layer@c,np.int32(nf),np.float64(density),out))
        result=out.get()
        if not np.all(np.isfinite(result)):
            raise FloatingPointError('nonfinite linearized energy derivative')
        return result


class RawKernelAnchor:
    """One geometry anchor, no history, fitting, eviction or silent refresh.

    ``max_pair_bytes`` bounds retained raw values/slopes (not factor/output or
    temporary workspaces). It is checked before allocating derivative tensors.
    The explicit immutable chart token is an ownership contract from the caller.
    """
    def __init__(self,points,faces,*,chart_token,order=3,max_pair_bytes=256*1024**2):
        points,faces=_validate(points,faces)
        if chart_token is None:
            raise ValueError('explicit chart token required')
        hash(chart_token)
        if not isinstance(order,(int,np.integer)) or order<1:
            raise ValueError('positive integer quadrature order required')
        nf=len(faces);required=nf*nf*4*19*8
        if not isinstance(max_pair_bytes,(int,np.integer)) or max_pair_bytes<0:
            raise ValueError('nonnegative integer raw-pair byte limit required')
        if required>max_pair_bytes:
            raise MemoryError(f'raw anchor needs {required} pair bytes; limit {max_pair_bytes}')
        self.points=np.array(points,copy=True);self.faces=np.array(faces,dtype=np.int64,copy=True)
        self.points.setflags(write=False);self.faces.setflags(write=False)
        self.chart_token=chart_token;self.order=int(order);self.pair_bytes=required
        self.device_faces=cp.asarray(self.faces)
        geometry,weights,_,_=_geometry_jacobians(points,faces,order=order)
        self.value=cp.empty((nf,nf,4),dtype=cp.float64)
        self.jacobian=cp.empty((nf,nf,18,4),dtype=cp.float64)
        _build(((nf*nf*6+127)//128,),(128,),
            (*geometry,weights,np.int32(nf),np.int32(len(weights)),self.value,self.jacobian))
        if not bool((cp.all(cp.isfinite(self.value)) & cp.all(cp.isfinite(self.jacobian))).get()):
            raise FloatingPointError('nonfinite anchor pair values/derivatives')
        incidence=np.argsort(self.faces.ravel(),kind='stable').astype(np.int32)
        pointer=np.r_[0,np.cumsum(np.bincount(self.faces.ravel(),minlength=len(points)))].astype(np.int32)
        self.incident=tuple(cp.asarray(a) for a in (pointer,incidence//3,incidence%3))

    def evaluate(self,points,faces,*,chart_token):
        """Guard first externally; reject changed chart/faces or non-SPD V here."""
        points,faces=_validate(points,faces)
        if (chart_token!=self.chart_token or points.shape!=self.points.shape
                or not np.array_equal(faces,self.faces)):
            raise ValueError('raw anchor requires identical semantic chart/faces/winding')
        # Existing analytic geometry owner supplies exact current curls/Jacobian.
        # Its packed pair geometry is temporary, not retained or used by hits.
        geometry,_,curl,curl_jac=_geometry_jacobians(points,faces,order=self.order)
        del geometry
        nf,nv=len(faces),len(points)
        raw=cp.empty((nf,nf),dtype=cp.float64);h=cp.zeros((nf,nv),dtype=cp.float64)
        _update(((nf*nf+127)//128,),(128,),
            (self.value,self.jacobian,cp.asarray(np.ascontiguousarray(points-self.points)),
             self.device_faces,np.int32(nf),np.int32(nv),raw,h))
        v=(raw+raw.T)*.5
        with errstate(linalg='raise'):
            factor=cp.linalg.cholesky(v)
            z=solve_triangular(factor,h,lower=True)
        w=cp.empty((nv,nv),dtype=cp.float64)
        _curl(((nv*nv+127)//128,),(128,),
            (v,cp.asarray(curl),*self.incident,np.int32(nf),np.int32(nv),w))
        w+=z.T@z
        energy=w.get()
        if not np.all(np.isfinite(energy)):
            raise FloatingPointError('nonfinite linearized energy')
        owned=np.array(points,copy=True);owned.setflags(write=False);energy.setflags(write=False)
        diagnostics=dict(pair_bytes=self.pair_bytes,
            h_constant_norm=float(cp.max(cp.abs(h@cp.ones(nv))).get()),
            energy_constant_norm=float(np.max(np.abs(energy@np.ones(nv)))),
            max_displacement=float(np.max(np.abs(points-self.points))))
        return LinearizedAssembly(energy,diagnostics,factor,h,owned,v,curl,curl_jac,self)
