"""Explicit opt-in regression for the shared CPU/GPU solve and publication."""
import os
from pathlib import Path

import numpy as np
import pytest

from topwater.solve import simulate
from topwater.solve.outputs import save_solve_outputs
from topwater.workflow_data import Boat
from topwater.workflow_data import Result
from topwater.workflow_data import Sea


@pytest.mark.skipif(os.environ.get('TOPWATER_TEST_INTEGRATION') != '1',
                    reason='complete solve requires explicit integration opt-in')
@pytest.mark.parametrize('device', ['cpu', pytest.param('gpu', marks=pytest.mark.skipif(
    os.environ.get('TOPWATER_TEST_GPU') != '1', reason='GPU integration requires hardware opt-in'))])
def test_simple_boat_cold_warm_and_allocation_backoff(tmp_path, monkeypatch, capsys, device):
    from topwater.solve import bem_checks
    checks = []
    def observe(name, original):
        def check(problems):
            checks.append((name, max(float(problem.omega) for problem in problems
                                     if np.isfinite(float(problem.omega)))))
            return original(problems)
        return check
    for name in ('_check_wavelength_and_mesh_resolution', '_check_wavelength_and_irregular_frequencies'):
        monkeypatch.setattr(bem_checks, name, observe(name, getattr(bem_checks, name)))
    root = Path(__file__).resolve().parents[1]
    example = root/'examples/simple-boat-native'
    source_cache = tmp_path/'sources'
    monkeypatch.setenv('TOPWATER_SOURCE_CACHE_DIR', str(source_cache))
    boat = Boat.load(example/'boat.npz')
    sea = Sea.load(example/'sea.npz')
    with np.load(root/'examples/simple-boat-wecsim/result.npz', allow_pickle=False) as archive:
        reference_time = archive['time_s']
    cold = None
    for name in ('cold', 'warm', 'backoff'):
        checks.clear()
        allocation_failures = []
        if name == 'backoff' and device == 'gpu':
            from topwater.solve.gpu.memory import ReusableBuffers
            from topwater.solve.gpu.memory import GPUAllocationError
            original_get = ReusableBuffers.get
            fft_length = 1 << (2*len(reference_time)-2).bit_length()
            def limited(self, buffer_name, size):
                if buffer_name == 'work_b' and size > 16*fft_length*16:
                    allocation_failures.append(size//(fft_length*16))
                    raise GPUAllocationError('injected allocation failure')
                return original_get(self, buffer_name, size)
            monkeypatch.setattr(ReusableBuffers, 'get', limited)
        elif name == 'backoff':
            from topwater.solve.cpu.fourier import CPUFourierWorkspace as FourierWorkspace
            original_reserve = FourierWorkspace.reserve
            def limited(self, n, rows, outputs):
                original_reserve(self, n, rows, outputs)
                if rows > 16:
                    allocation_failures.append(rows)
                    raise MemoryError('injected host allocation failure')
            monkeypatch.setattr(FourierWorkspace, 'reserve', limited)
        result = simulate(boat, sea, 100., device=device)
        maximum=2*np.pi*np.fft.rfftfreq(len(result.time_s),.01)[-1]
        assert checks == [('_check_wavelength_and_mesh_resolution', maximum),
                          ('_check_wavelength_and_irregular_frequencies', maximum)]
        np.testing.assert_array_equal(result.time_s, reference_time)
        assert result.surface_interval_s == .1
        assert result.view_bounds_xy_m.shape == (2, 2)
        if cold is None:
            cold = result
        else:
            for member in Result.members:
                np.testing.assert_allclose(getattr(result, member), getattr(cold, member), rtol=0., atol=1e-10, err_msg=member)
        directory = tmp_path/name
        destination, plot = save_solve_outputs(result, directory, sea=sea)
        restored = Result.load(destination)
        np.testing.assert_array_equal(restored.surface_vertices_m, result.surface_vertices_m)
        assert plot.stat().st_size > 0
        assert list(source_cache.glob('*.npz'))
        output = capsys.readouterr().err
        assert output.count('BEM_PREFLIGHT ') == 1
        assert f'SOURCE_CACHE {"stored" if name == "cold" else "hit"}' in output
        if name == 'backoff':
            assert allocation_failures == ([64, 32] if device == 'gpu' else [128, 64, 32])
            assert 'retrying batch_points=32' in output
            assert 'retrying batch_points=16' in output
            assert 'batch_points=16 reusable_bytes=' in output
