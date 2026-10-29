"""Host numerical operations for the shared reconstruction pipeline."""
from topwater.solve.fields import FieldBackend
from topwater.solve.cpu.fourier import CPUFourierWorkspace
from topwater.solve.cpu.potential import CPUPotential


class CPUFieldBackend(FieldBackend):
    device = 'cpu'
    radiation_chunk_points = 128

    def __init__(self, hydro, *, cache_sources=True, threads=None):
        evaluator = CPUPotential(hydro.solver.engine.green_function, hydro.body.mesh_including_lid, threads=threads)
        super().__init__(hydro, evaluator, CPUFourierWorkspace(), cache_sources=cache_sources)
