"""Fourier grid, cache identity and transform preparation for both devices."""
import numpy as np
from .settings import OfflineRadiationSettings
from .validation import _finite, _grid


def fft_dimensions(sample_count,native_dt,settings):
    """One transform-length policy for workspace allocation and field synthesis."""
    native=settings.period_s is None
    dt=(native_dt if native else native_dt/2) if settings.fft_dt_s is None else float(settings.fft_dt_s)
    if not np.isfinite(dt) or dt<=0:
        raise ValueError('Offline radiation requires a positive finite FFT timestep')
    if native:
        return sample_count,dt
    period=float(settings.period_s)
    if not np.isfinite(period) or period<=0:
        raise ValueError('Explicit padded FFT period must be positive and finite')
    n=int(round(period/dt))
    if n<2:
        raise ValueError('Offline FFT period/grid/band must cover the requested history')
    return n,dt


class FourierPreparation:
    def __init__(self):
        self._key = self._values = None

    def reset_preparation(self):
        self._key = self._values = None

    def prepare(self, w, values, n, dt):
        key = (n, dt, w.shape, w.tobytes())
        if key != self._key or not np.array_equal(values, self._values):
            self.f = 2*np.pi*np.fft.rfftfreq(n, dt)
            if w.shape != self.f.shape or not np.array_equal(w, self.f):
                raise ValueError('Radiation transfer frequencies must equal the velocity FFT bins; no frequency extrapolation or interpolation')
            self.prepare_transform(values, n)
            self._key, self._values = key, np.array(values, copy=True)

    def prepare_transform(self, values, n):
        raise NotImplementedError


def surface_radiation_frequencies(time_s, settings=None):
    """Actual velocity FFT bins, shared by initial preflight and reconstruction."""
    settings = OfflineRadiationSettings() if settings is None else settings
    time = _grid(time_s, 'motion time_s', zero=True)
    n, dt = fft_dimensions(len(time), time[1]-time[0], settings)
    return 2*np.pi*np.fft.rfftfreq(n, dt)


def offline_radiation_response(omega, velocity_transfer, time_s, velocity, acceleration,
                               *, output_time_s=None, settings=None, workspace=None):
    """Full complex H_eta,v applied to the complete saved native velocity.

    H uses exp(-i*w*t); NumPy's inverse uses exp(+i*w*t), hence conj(H).
    By default the transform contains only the saved native velocity samples.
    Explicit period_s settings retain the historical padded/Hermite policy.
    Transfers are evaluated at every velocity FFT bin. These choices are logged,
    NOT a universal causality claim
    or an automatic numerical-error certificate. A*a is already part of full H:
    never add it a second time or force its first sample to zero.
    """
    settings = OfflineRadiationSettings() if settings is None else settings
    stamp = _grid(time_s, 'motion time_s', zero=True)
    w = _grid(omega, 'transfer omega', zero=True)
    if workspace is None:
        from topwater.solve.cpu.fourier import CPUFourierWorkspace
        workspace = CPUFourierWorkspace()
    H = workspace.prepare_transfer(velocity_transfer, len(w))
    v = _finite(velocity, (len(stamp), 6), 'native velocity')
    a = _finite(acceleration, v.shape, 'native acceleration')
    output = stamp if output_time_s is None else _grid(output_time_s, 'output time_s', zero=True)
    native_dt = stamp[1]-stamp[0]
    native = settings.period_s is None
    n, dt = fft_dimensions(len(stamp), native_dt, settings)
    future = settings.future_s
    if (not np.isfinite([dt, future]).all() or dt <= 0 or future < 0
            or stamp[0] != 0 or w[0] != 0 or np.any(v[0] != 0)
            or not np.allclose(np.diff(stamp), native_dt, rtol=0, atol=1e-12)):
        raise ValueError('Offline radiation requires finite settings, uniform zero-start history and exact initial velocity rest')
    if native:
        if abs(dt-native_dt) > 1e-12 or future != 0:
            raise ValueError('Native FFT requires saved timestep and no future continuation')
        zero = 0
        period = n*dt  # DFT period includes one sample interval beyond the last knot.
        selected = np.rint(output/dt).astype(int)
        if output[0] < 0 or output[-1] > stamp[-1]:
            raise ValueError('Offline FFT output must lie within the saved history')
        if not np.allclose(selected*dt, output, rtol=0, atol=1e-10):
            raise ValueError('Output clocks must be exactly represented by the native FFT dt')
        values = v  # No interpolation, oversampling, padding or invented future.
    else:
        period = float(settings.period_s)
        if (n % 4 or abs(n*dt-period) > 1e-10
                or output[0] < 0 or output[-1] > stamp[-1]):
            raise ValueError('Offline FFT period/grid/band must cover the requested history without frequency extrapolation')
        knot_indices = np.rint(stamp/dt).astype(int)
        selected = np.rint(output/dt).astype(int)
        steps = int(round(future/dt))
        if (not np.allclose(knot_indices*dt, stamp, rtol=0, atol=1e-10)
                or not np.allclose(selected*dt, output, rtol=0, atol=1e-10)
                or steps < 0 or abs(steps*dt-future) > 1e-10):
            raise ValueError('Saved motion, output clocks and future duration must be exactly represented by FFT dt')
        zero = n//4
        last = int(knot_indices[-1])
        if zero+last+steps >= n:
            raise ValueError('FFT period does not fit the saved history and declared smooth future')
        values = np.zeros((n, 6))
        inside = np.arange(last+1)*dt
        for mode in range(6):
            values[zero:zero+last+1, mode] = np.interp(inside, stamp, v[:, mode])
        values[zero+knot_indices] = v  # Preserve every native sample exactly.
        if not np.array_equal(values[zero+knot_indices], v):
            raise ValueError('Offline interpolation changed saved velocity samples')
        u = np.arange(1, steps+1)/steps if steps else np.empty(0)
        values[zero+last+1:zero+last+steps+1] = (
            (2*u**3-3*u**2+1)[:, None]*v[-1]
            +future*(u**3-2*u**2+u)[:, None]*a[-1])
    # Accumulate modes into one output spectrum: no FFT-frequency x point x DOF
    # tensor. Spatial chunking by the caller bounds this remaining allocation.
    spectrum = workspace.spectrum(H,w,values,n,dt)
    result,raw_initial,pre_max = workspace.inverse_response(spectrum,n,zero,selected)
    evidence = dict(method='intact full-complex offline FFT(v)',
        period_s=float(period), fft_dt_s=float(dt), n_fft=n, fft_time_origin_s=float(-zero*dt),
        maximum_frequency_hz=float(w[-1]/(2*np.pi)), frequency_spacing_hz=1/(n*dt),
        frequency_count=len(w), transfer_interpolation='none; exact velocity FFT bins',
        saved_motion_interpolation=('none; exact native samples' if native else 'linear with every original native velocity knot exact'),
        prehistory=('none; unpadded periodic transform' if native else 'zero native velocity'), future_duration_s=float(future),
        future_policy=('none' if native else 'C1 Hermite from saved final velocity/acceleration to zero velocity/acceleration'),
        pre_response_max_abs_m=(None if native else pre_max),
        temporal_policy=('native unpadded' if native else 'explicit historical padded'),
        raw_initial_radiation_m=raw_initial, initial_sample_modified=False,
        direct_term_policy='intact H, no separately added acceleration; raw initial value preserved',
        universal_causality_claim=False, numerical_accuracy_qualified=False,
        unresolved_response_above_nyquist=True)
    return result, evidence
