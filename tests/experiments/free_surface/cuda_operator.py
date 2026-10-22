"""Float64 CUDA version of the existing experimental P1/P0 energy.

Use the optional CuPy CUDA runtime in the baseline command. No single-precision
operator, row correction, physical approximation, or fast-math flag is used.
"""
import ctypes
from importlib.metadata import distribution,distributions,PackageNotFoundError
import os
from pathlib import Path
import numpy as np
from tests.laplace_reference import triangle_samples


# uv's Windows overlay can hide the component-wheel DLLs from path-finder.
# Resolve those installed wheels in this process; leave system CUDA alone when
# no component runtime is installed. Keep directory/library handles alive.
_handles=[]
if os.name=='nt':
    try:
        runtime=distribution('nvidia-cuda-runtime-cu12')
    except PackageNotFoundError:
        runtime=None
    if runtime is not None:
        os.environ['CUDA_PATH']=str(runtime.locate_file('nvidia/cuda_runtime'))
        dlls=[Path(package.locate_file(file)) for package in distributions()
            if package.metadata['Name'].startswith('nvidia-') for file in package.files or ()
            if str(file).endswith('.dll')]
        directories=sorted({str(file.parent) for file in dlls})
        os.environ['PATH']=os.pathsep.join(directories+[os.environ['PATH']])
        _handles.extend(os.add_dll_directory(directory) for directory in directories)
        nvrtc=next((file for file in dlls if file.name.startswith('nvrtc64_')),None)
        if nvrtc is not None:
            _handles.append(ctypes.WinDLL(str(nvrtc)))

import cupy as cp
from cupyx import errstate
from cupyx.scipy.linalg import solve_triangular

_source=(Path(__file__).with_name('cuda_pair.cuh').read_text()+'\n'+
    Path(__file__).with_name('cuda_energy.cu').read_text())
_pairs=cp.RawKernel(_source,'energy_pairs',options=('--std=c++11',))
_curl=cp.RawKernel(_source,'curl_energy',options=('--std=c++11',))


def assemble_energy(faces,panels,order=3):
    points,weights=triangle_samples(panels.triangles,order)
    opposite=np.roll(panels.triangles,-2,axis=1)-np.roll(panels.triangles,-1,axis=1)
    gradient=np.cross(panels.geometric_normals[:,None],opposite)/(2*panels.areas[:,None,None])
    curl=np.cross(panels.geometric_normals[:,None],gradient)
    nf,nv,nq=len(panels),int(faces.max())+1,len(weights)
    geometry=[cp.asarray(np.ascontiguousarray(array)) for array in (
        panels.triangles,panels.geometric_normals,panels.edge_outward,panels.tangents,
        panels.edge_lengths,panels.areas,panels.centroids,panels.normal_sign,
        faces.astype(np.int64),gradient,points,weights)]
    v=cp.empty((nf,nf),dtype=cp.float64);h=cp.zeros((nf,nv),dtype=cp.float64)
    _pairs(((nf*nf+127)//128,),(128,),(*geometry,np.int32(nf),np.int32(nv),np.int32(nq),v,h))
    asymmetry=float((cp.max(cp.abs(v-v.T))/cp.max(cp.abs(v))).get())
    v=(v+v.T)*.5
    order=np.argsort(faces.ravel(),kind='stable').astype(np.int32)
    pointer=np.r_[0,np.cumsum(np.bincount(faces.ravel(),minlength=nv))].astype(np.int32)
    incident=[cp.asarray(array) for array in (pointer,order//3,order%3)]
    w=cp.empty((nv,nv),dtype=cp.float64)
    _curl(((nv*nv+127)//128,),(128,),
        (v,cp.asarray(curl),*incident,np.int32(nf),np.int32(nv),w))
    with errstate(linalg='raise'):
        factor=cp.linalg.cholesky(v)
        z=solve_triangular(factor,h,lower=True)
    w+=z.T@z
    constant=float(cp.max(cp.abs(h@cp.ones(nv))).get())
    return w.get(),dict(v_quadrature_asymmetry=asymmetry,h_constant_norm=constant,
        operator_bytes=dict(single_layer_factor=v.nbytes,mixed_trace=h.nbytes,energy=w.nbytes))
