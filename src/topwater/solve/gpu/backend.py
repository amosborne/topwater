"""GPU resource construction for the shared reconstruction pipeline."""
from importlib.metadata import version
import sys

import numpy as np

from topwater.solve.field_backend import GPUUnavailableError
from topwater.solve.gpu.potential import OpenCLPotential
from topwater.solve.fields import FieldBackend


class GPUFieldBackend(FieldBackend):
    device = 'gpu'
    # Smallest measured batch matching peak simple-boat throughput on gfx1201.
    radiation_chunk_points = 64

    def __init__(self, hydro, *, device_name=None, cache_sources=True):
        if version('capytaine') != '3.0.0' or hydro.depth != np.inf:
            raise GPUUnavailableError('GPU fields require Capytaine 3.0.0 and infinite water depth')
        initialized = False
        try:
            self.evaluator = OpenCLPotential(hydro.solver.engine.green_function,
                hydro.body.mesh_including_lid, device_name=device_name)
            from topwater.solve.gpu.fourier import GPUFourierWorkspace
            self.fourier = GPUFourierWorkspace(self.evaluator)
            super().__init__(hydro, self.evaluator, self.fourier, cache_sources=cache_sources)
            initialized = True
        except ImportError as error:
            raise GPUUnavailableError('install topwater[gpu] (uv sync --extra gpu) to enable GPU reconstruction') from error
        except (RuntimeError, ValueError) as error:
            raise GPUUnavailableError(f'GPU initialization failed: {error}') from error
        except Exception as error:
            # PyOpenCL's driver errors do not all inherit built-in RuntimeError.
            import pyopencl as cl
            if isinstance(error, cl.Error):
                raise GPUUnavailableError(f'OpenCL initialization failed: {error}') from error
            raise
        finally:
            if not initialized and hasattr(self, 'evaluator'):
                primary = sys.exception()
                try:
                    self.evaluator.close()
                except Exception as error:
                    if primary is not None:
                        primary.add_note(f'GPU initialization cleanup also failed: {error}')
                    else:
                        raise
        self.device_name = self.evaluator.device.name
        self.driver_version = self.evaluator.device.driver_version
