"""Flat full-plane compact P1 WEAK DtN by aperiodic FFT convolution.

Static experimental operator only. It is K=2*B.T*V*B, not M^-1*K or a point
velocity. No physical hull hole, curved surface or hidden-water evolution.
"""
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from numbers import Real,Integral
import numpy as np
from scipy.fft import next_fast_len
from scipy.sparse import coo_matrix
from topwater.rankine import TrianglePanels,QuadratureError
from topwater.numerics import LaplaceControls
from tests.laplace_reference import triangle_samples
from .cuda_operator import cp

_source=Path(__file__).with_name('cuda_pair.cuh').read_text()+'\n'+Path(__file__).with_suffix('.cu').read_text()
_lags=cp.RawKernel(_source,'weak_hat_lags',options=('--std=c++11','--fmad=false'))


def hat(hx,hy):
    """Complete six-triangle BL→TR hat, owned by integer cell/orientation IDs."""
    cells=np.array([[0,0],[0,0],[-1,-1],[-1,-1],[-1,0],[0,-1]],dtype=np.int32)
    orientation=np.array([0,1,0,1,0,1],dtype=np.int32)
    offsets=np.array([[[0,0],[1,0],[1,1]],[[0,0],[1,1],[0,1]]])
    lattice=cells[:,None,:]+offsets[orientation]
    triangles=np.concatenate((lattice*np.array([hx,hy]),np.zeros((6,3,1))),axis=2)
    values=np.all(lattice==0,axis=2).astype(float)
    panels=TrianglePanels(triangles,1)
    gradient=np.linalg.solve(triangles[:,1:,:2]-triangles[:,:1,:2],
                             (values[:,1:]-values[:,:1])[...,None])[...,0]
    return panels,gradient,values,cells,orientation


@dataclass(frozen=True)
class WeakControls:
    relative_tolerance: float=1e-5
    absolute_tolerance: float=1e-10  # Multiplied by min(hx,hy), since K has length units.
    maximum_order: int=64
    maximum_bytes: int=LaplaceControls().max_operator_bytes
    maximum_evaluations: int=500_000_000

    def __post_init__(self):
        for name in ('relative_tolerance','absolute_tolerance'):
            value=getattr(self,name)
            if isinstance(value,bool) or not isinstance(value,Real) or not np.isfinite(value) or not 0<value<1:
                raise ValueError('weak error requests must be finite and between zero and one')
        if (isinstance(self.maximum_order,bool) or not isinstance(self.maximum_order,Integral)
                or self.maximum_order not in (4,8,16,32,64)
                or any(isinstance(v,bool) or not isinstance(v,Integral) or v<1
                       for v in (self.maximum_bytes,self.maximum_evaluations))):
            raise ValueError('bounded integer order and positive storage budget required')


def _storage(shape):
    ny,nx=shape;count=ny*nx;lags=(2*ny-1)*(2*nx-1)
    fft_shape=tuple(next_fast_len(2*n-1) for n in shape)
    nonzeros=count+2*ny*(nx-1)+2*nx*(ny-1)+2*(nx-1)*(ny-1)
    spectrum_bytes=16*fft_shape[0]*(fft_shape[1]//2+1)
    # 64-bit CSR indices bound either SciPy index width. Excludes construction,
    # application scratch, FFT workspace/cache and independent references.
    retained=8*lags+spectrum_bytes+16*nonzeros+8*(count+1)
    return dict(fft_shape=fft_shape,retained_operator_bytes_upper_bound=retained,
        storage_scope='Host lag kernel+consistent sparse mass+device spectrum; excludes scratch/workspace/caches')


def lag_kernel(shape,hx,hy,controls=WeakControls()):
    if (len(shape)!=2 or any(isinstance(n,bool) or not isinstance(n,Integral) or not 1<=n<2**31 for n in shape)
            or not np.isfinite([hx,hy]).all() or min(hx,hy)<=0
            or not isinstance(controls,WeakControls)):
        raise ValueError('positive integer shape, finite spacings and bounded weak controls required')
    shape=tuple(int(n) for n in shape)
    storage=_storage(shape)
    if storage['retained_operator_bytes_upper_bound']>controls.maximum_bytes:
        raise MemoryError('weak retained operators exceed storage budget; excludes scratch/workspace/caches')
    started=perf_counter()
    ny,nx=shape
    # Only central inversion is assumed; diagonal triangulation is not radial.
    canonical=np.array([(dx,dy) for dy in range(ny) for dx in range(-nx+1,nx)
                       if dy>0 or dx>=0],dtype=np.int32)
    panels,gradient,nodal,cells,orientation=hat(hx,hy)
    geometry=[cp.asarray(np.ascontiguousarray(a)) for a in
        (panels.triangles,panels.geometric_normals,panels.edge_outward,panels.tangents,
         panels.edge_lengths,panels.areas,panels.centroids,panels.normal_sign,gradient)]
    own=cp.asarray(cells),cp.asarray(orientation)
    previous=np.full(len(canonical),np.nan);accepted=np.full(len(canonical),np.nan)
    final_order=np.zeros(len(canonical),dtype=int);asymmetry=np.zeros(len(canonical))
    pending=np.arange(len(canonical));work=[];evaluations=0
    order=2
    while order<=controls.maximum_order and len(pending):
        separation=np.maximum(abs(canonical[pending])-2,0)*np.array([hx,hy])
        far_count=int(np.count_nonzero(np.sum(separation**2,axis=1)>=hx*hx+hy*hy))
        near_count=len(pending)-far_count
        near_samples=near_count*72*order**2
        far_terms=far_count*36*order**4
        if evaluations+near_samples+far_terms>controls.maximum_evaluations:
            raise QuadratureError('weak lag integration exceeded its bounded sample/term budget')
        evaluations+=near_samples+far_terms
        points,weights=triangle_samples(panels.triangles,order)
        bary,_=triangle_samples(np.eye(3)[None],order)
        phi=np.einsum('qj,tj->tq',bary[0],nodal)
        hat_weights=(panels.areas[:,None]*weights[None,:]*phi).ravel()
        values=cp.empty(len(pending),dtype=cp.float64);defect=cp.empty_like(values);roundoff=cp.empty_like(values)
        _lags(((len(pending)+63)//64,),(64,),(*geometry,cp.asarray(points),cp.asarray(weights),
            *own,np.int32(len(weights)),cp.asarray(np.ascontiguousarray(points[:,:,:2].reshape(-1,2))),
            cp.asarray(hat_weights),np.int32(len(hat_weights)),cp.asarray(canonical[pending]),
            np.int32(len(pending)),np.float64(hx),np.float64(hy),values,defect,roundoff))
        current,raw_asymmetry=cp.asnumpy(values),cp.asnumpy(defect)
        tolerance=controls.absolute_tolerance*min(hx,hy)+controls.relative_tolerance*abs(current)
        error=abs(current-previous[pending])+raw_asymmetry/2+cp.asnumpy(roundoff)
        good=np.isfinite(current)&np.isfinite(error)&(error<=tolerance)
        selected=pending[good]
        accepted[selected]=current[good];final_order[selected]=order;asymmetry[selected]=raw_asymmetry[good]
        previous[pending]=current
        work.append(dict(order=order,evaluated=len(pending),accepted=int(good.sum()),
            near_source_samples_upper_bound=near_samples,far_inverse_cube_terms=far_terms))
        pending=pending[~good];order*=2
    if len(pending):
        failure=QuadratureError(f'weak lag quadrature unresolved for {len(pending)} offsets at order {controls.maximum_order}')
        failure.weak_diagnostics=dict(quadrature_work=work,unresolved=[dict(
            lag=lag.tolist(),value=float(value),estimated_error=float(e),tolerance=float(t),raw_asymmetry=float(a))
            for lag,value,e,t,a in zip(canonical[pending][:20],current[~good][:20],error[~good][:20],
                                      tolerance[~good][:20],raw_asymmetry[~good][:20],strict=True)])
        raise failure
    kernel=np.empty((2*ny-1,2*nx-1))
    for (dx,dy),value in zip(canonical,accepted,strict=True):
        kernel[dy+ny-1,dx+nx-1]=kernel[-dy+ny-1,-dx+nx-1]=value
    return kernel,dict(assembly_seconds=perf_counter()-started,quadrature_work=work,
        **storage,
        maximum_order=int(final_order.max()),maximum_raw_asymmetry=float(asymmetry.max()),
        sample_term_work_upper_bound=evaluations,maximum_evaluations=controls.maximum_evaluations,
        relative_tolerance=controls.relative_tolerance,absolute_tolerance_scaled=controls.absolute_tolerance*min(hx,hy),
        scope='Target/source symmetric near quadrature; whole-hat separated inverse-cube quadrature; order agreement is not a certified integral bound')


def mass_matrix(shape,hx,hy):
    """Consistent P1 mass on interior nodes with zero outer perimeter."""
    ny,nx=shape
    yy,xx=np.indices(shape);rows=[];cols=[];data=[]
    for dx,dy,weight in [(0,0,.5),(1,0,1/12),(-1,0,1/12),(0,1,1/12),(0,-1,1/12),(1,1,1/12),(-1,-1,1/12)]:
        good=(xx+dx>=0)&(xx+dx<nx)&(yy+dy>=0)&(yy+dy<ny)
        rows.extend((yy[good]*nx+xx[good]).tolist())
        cols.extend(((yy[good]+dy)*nx+xx[good]+dx).tolist())
        data.extend(np.full(int(good.sum()),weight*hx*hy))
    return coo_matrix((data,(rows,cols)),shape=(nx*ny,nx*ny)).tocsr()


class WeakPatch:
    """One finite weak Toeplitz block; apply returns weak loads, not velocities."""
    def __init__(self,shape,hx,hy,controls=WeakControls()):
        kernel,self.diagnostics=lag_kernel(shape,hx,hy,controls)
        self.shape=tuple(int(n) for n in shape);self.kernel=kernel
        kernel.setflags(write=False)
        self.mass=mass_matrix(shape,hx,hy)
        for array in (self.mass.data,self.mass.indices,self.mass.indptr): array.setflags(write=False)
        self.fft_shape=self.diagnostics['fft_shape']
        embed=cp.zeros(self.fft_shape,dtype=cp.float64)
        ny,nx=shape;dy,dx=np.indices(kernel.shape)
        embed[cp.asarray((dy-ny+1)%self.fft_shape[0]),cp.asarray((dx-nx+1)%self.fft_shape[1])]=cp.asarray(kernel)
        self._spectrum=cp.fft.rfft2(embed)
        self.diagnostics.update(kernel_bytes=kernel.nbytes,fft_shape=self.fft_shape,spectrum_bytes=self._spectrum.nbytes,
            mass_nonzeros=self.mass.nnz)

    def apply(self,values):
        values=cp.asarray(values)
        if values.shape!=self.shape or values.dtype!=cp.float64 or not bool(cp.isfinite(values).all()):
            raise ValueError('one finite float64 coefficient per interior patch node required')
        padded=cp.zeros(self.fft_shape,dtype=cp.float64)
        padded[:self.shape[0],:self.shape[1]]=values
        return cp.fft.irfft2(cp.fft.rfft2(padded)*self._spectrum,s=self.fft_shape)[:self.shape[0],:self.shape[1]]
