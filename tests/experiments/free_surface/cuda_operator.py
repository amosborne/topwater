"""Float64 CUDA version of the existing experimental P1/P0 energy.

Use the optional CuPy CUDA runtime in the baseline command. No single-precision
operator, row correction, physical approximation, or fast-math flag is used.
"""
import ctypes
from dataclasses import dataclass
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
_uniform_pairs=_pairs
_curl=cp.RawKernel(_source,'curl_energy',options=('--std=c++11',))


@dataclass(frozen=True, eq=False)
class EnergyAssembly:
    """One CURRENT energy and its matching GPU single-layer factor/trace.

    Retention is optional. This is not a stale-operator cache: force reuse checks
    exact coordinates, face numbering, winding and quadrature before a solve.
    """
    energy: np.ndarray
    diagnostics: dict
    factor: object
    mixed_trace: object
    faces: np.ndarray
    triangles: np.ndarray
    normal_sign: np.ndarray
    quadrature_order: int
    uniform_quadrature: bool

    def auxiliary(self, points, faces, potential, *, order=3):
        if (self.quadrature_order != order or not self.uniform_quadrature
                or potential.shape != (len(self.energy),)
                or not np.array_equal(faces, self.faces)
                or not np.array_equal(points[faces], self.triangles)
                or np.any(self.normal_sign != 1)):
            raise ValueError('force reuse requires matching quadrature, current geometry and winding')
        with errstate(linalg='raise'):
            value = solve_triangular(self.factor, self.mixed_trace@cp.asarray(potential), lower=True)
            return solve_triangular(self.factor.T, value, lower=False)


def prepare_energy(faces,panels,order=3):
    """Assemble once, retaining current V factor/H for matching shape forces."""
    pair_kernel=_pairs
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
    pair_kernel(((nf*nf+127)//128,),(128,),(*geometry,np.int32(nf),np.int32(nv),np.int32(nq),v,h))
    asymmetry=float((cp.max(cp.abs(v-v.T))/cp.max(cp.abs(v))).get())
    v=(v+v.T)*.5
    incident_order=np.argsort(faces.ravel(),kind='stable').astype(np.int32)
    pointer=np.r_[0,np.cumsum(np.bincount(faces.ravel(),minlength=nv))].astype(np.int32)
    incident=[cp.asarray(array) for array in (pointer,incident_order//3,incident_order%3)]
    w=cp.empty((nv,nv),dtype=cp.float64)
    _curl(((nv*nv+127)//128,),(128,),
        (v,cp.asarray(curl),*incident,np.int32(nf),np.int32(nv),w))
    with errstate(linalg='raise'):
        factor=cp.linalg.cholesky(v)
        z=solve_triangular(factor,h,lower=True)
    w+=z.T@z
    constant=float(cp.max(cp.abs(h@cp.ones(nv))).get())
    diagnostics=dict(v_quadrature_asymmetry=asymmetry,h_constant_norm=constant,
        operator_bytes=dict(single_layer_factor=factor.nbytes,mixed_trace=h.nbytes,energy=w.nbytes))
    owned = [np.array(value, copy=True) for value in (faces, panels.triangles, panels.normal_sign)]
    for value in owned:
        value.setflags(write=False)
    energy = w.get()
    energy.setflags(write=False)
    return EnergyAssembly(energy, diagnostics, factor, h, *owned, order,
                          pair_kernel is _uniform_pairs)


def assemble_energy(faces,panels,order=3):
    """Legacy nonretaining interface; use prepare_energy for a field/force pair."""
    current = prepare_energy(faces, panels, order)
    # The legacy interface exposes no retained factor, and keeps its original
    # writable host-array behavior without another matrix copy.
    current.energy.setflags(write=True)
    return current.energy, current.diagnostics
