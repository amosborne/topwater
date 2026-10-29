"""Optional GPU field kernels; importing does not initialize a GPU runtime."""
from topwater.solve.gpu.potential import OpenCLPotential

__all__ = ['OpenCLPotential']
