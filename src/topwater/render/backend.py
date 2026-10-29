"""Shared geometry and selection of native CPU or hardware GPU drawing."""
import numpy as np


class DrawingGeometry:
    fullscreen=np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,-1.],[1.,1.],[-1.,1.]],np.float32)

    def to_clip(self,xy):
        return ((np.asarray(xy)-self.lower)/(self.upper-self.lower)*2-1).astype(np.float32)

    def vertices(self,xy,color):
        xy=np.asarray(xy).reshape(-1,2)
        return np.c_[xy,np.broadcast_to(color,(len(xy),4))].astype(np.float32)


def DrawingBackend(size,lower,upper,*,device='auto'):
    if device not in ('auto','cpu','gpu'):
        raise ValueError("Render device must be 'auto', 'cpu', or 'gpu'")
    if device!='cpu':
        try:
            from .gpu.backend import DrawingBackend as GPU
            return GPU(size,lower,upper)
        except (ImportError,RuntimeError) as error:
            if device=='gpu':
                raise RuntimeError('Hardware WebGPU rendering is unavailable; install topwater[gpu] and a compatible graphics driver') from error
    from .cpu.backend import DrawingBackend as CPU
    return CPU(size,lower,upper)
