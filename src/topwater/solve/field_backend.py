"""Select per-solve reconstruction resources without global patching."""
import sys


class GPUUnavailableError(RuntimeError):
    """GPU reconstruction cannot be initialized for this solve."""


def _cpu_backend(hydro, reason=None, *, cpu_workers=None):
    from topwater.solve.cpu.backend import CPUFieldBackend
    backend = CPUFieldBackend(hydro, threads=cpu_workers)
    reason = '' if reason is None else f' reason={reason}'
    print(f'FIELD_BACKEND device=cpu{reason} kernel=potential-only-fp64 '
          f'threads={backend.evaluator.threads}', file=sys.stderr)
    return backend


def select_field_backend(hydro, device='auto', *, cpu_workers=None):
    """Select numerical operations for one shared reconstruction pipeline."""
    if not isinstance(device, str) or not device.strip():
        raise ValueError('device must be auto, cpu, gpu, or a GPU name')
    device = device.strip()
    if device == 'cpu':
        return _cpu_backend(hydro, cpu_workers=cpu_workers)
    selector = None if device in ('auto', 'gpu') else device
    from topwater.solve.gpu.backend import GPUFieldBackend
    try:
        backend = GPUFieldBackend(hydro, device_name=selector)
    except GPUUnavailableError as error:
        if device != 'auto':
            raise
        return _cpu_backend(hydro, error, cpu_workers=cpu_workers)
    print(f'FIELD_BACKEND device=gpu name={backend.device_name} '
          f'driver={backend.driver_version}', file=sys.stderr)
    return backend
