from types import SimpleNamespace

import numpy as np
import pytest

from topwater.solve.gpu.backend import GPUFieldBackend
from topwater.solve.gpu.memory import GPUAllocationError
from topwater.solve.gpu.memory import ReusableBuffers
from topwater.solve.settings import OfflineRadiationSettings


def backend_with_limit(limit, error_type=GPUAllocationError):
    attempts, resets = [], []
    def reserve(n, rows, outputs):
        attempts.append((n, rows, outputs))
        if rows > limit:
            raise error_type('injected')
    backend = GPUFieldBackend.__new__(GPUFieldBackend)
    backend.fourier = SimpleNamespace(reserve=reserve, reset_buffers=lambda: resets.append('fft'),
        m=32768, nbytes=0)
    backend.evaluator = SimpleNamespace(reserve_radiation=lambda *args: None,
        reset_batch_buffers=lambda: resets.append('potential'), nbytes=0, allocation_error=GPUAllocationError)
    return backend, attempts, resets


def test_allocation_failure_halves_and_cleans_partial_workspace():
    backend, attempts, resets = backend_with_limit(16)
    result = backend.prepare_radiation(np.arange(13001)*.01, np.arange(1626)*.08,
        np.array([0., 48.]), OfflineRadiationSettings(chunk_points=64), 100)
    assert result == 16
    assert attempts == [(13001, 64, 1626), (13001, 32, 1626), (13001, 16, 1626)]
    assert resets == ['fft', 'potential', 'fft', 'potential']


def test_one_point_failure_stops_and_reports_context():
    backend, attempts, resets = backend_with_limit(0)
    with pytest.raises(RuntimeError, match='cannot fit a one-point batch'):
        backend.prepare_radiation(np.arange(10)*.01, [0., .08, .09], [0., 48.],
            OfflineRadiationSettings(chunk_points=5), 100)
    assert [attempt[1] for attempt in attempts] == [5, 2, 1]
    assert len(resets) == 6


def test_nonallocation_failure_is_not_retried():
    backend, attempts, resets = backend_with_limit(0, ValueError)
    with pytest.raises(ValueError, match='injected'):
        backend.prepare_radiation(np.arange(10)*.01, [0., .09], [0., 48.],
            OfflineRadiationSettings(chunk_points=64), 100)
    assert len(attempts) == 1
    assert resets == []


def test_named_buffers_are_reused_and_respect_allocation_limit():
    events = []
    class Buffer:
        def __init__(self, size):
            self.size = size
            events.append(('allocate', size))
        def release(self):
            events.append(('release', self.size))
    buffers = ReusableBuffers.__new__(ReusableBuffers)
    buffers.queue = SimpleNamespace(device=SimpleNamespace(max_mem_alloc_size=128), finish=lambda: None)
    buffers.values = {}
    buffers.allocator = Buffer
    first = buffers.get('work', 64)
    assert buffers.get('work', 32) is first
    buffers.get('work', 128)
    assert events == [('allocate', 64), ('release', 64), ('allocate', 128)]
    with pytest.raises(GPUAllocationError, match='allocation limit'):
        buffers.get('too_large', 129)
    buffers.clear()
    assert buffers.nbytes == 0


def test_potential_reservation_covers_larger_frequency_batches_at_short_tail():
    from topwater.solve.gpu.potential import OpenCLPotential
    evaluator = OpenCLPotential.__new__(OpenCLPotential)
    evaluator.nf = 32
    evaluator.matrix_budget = lambda: 4096
    sizes = {}
    evaluator.batch_buffers = SimpleNamespace(get=lambda name, size: sizes.update({name: size}))
    evaluator.reserve_radiation(3, 100)
    # Three-point blocks fit two frequencies; a two-point tail fits four.
    assert sizes['matrix'] == 4096
    assert sizes['wavenumbers'] == 100*8
