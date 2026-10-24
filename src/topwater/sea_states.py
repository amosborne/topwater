"""Observed sea spectra and synthetic, replayable directional realizations.

Measured band energies and four directional moments remain paired by observation.
Phases are unobserved: generated phases do not reconstruct a buoy elevation trace.
All propagation angles use WaveComponents' east/CCW convention and negative omega*t.
"""
from dataclasses import dataclass
from datetime import datetime, timedelta
from importlib.resources import files
import json
from numbers import Integral
import secrets
from types import MappingProxyType
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike

from ._validation import array
from .inputs import WaveComponents


@dataclass(frozen=True, eq=False)
class ObservedSeaState:
    """One observation: band cosine amplitudes [m], edges [Hz], FROM moments.

    Moments are (a1,b1,a2,b2) for bearings clockwise from true north. They constrain
    a directional reconstruction; they do not uniquely measure its angular PDF.
    Coverage outside frequency_edges_hz is unknown. Tm02 uses flat density within
    each measured/rebinned interval, not a frequency sampled for a realization.
    """
    id: str
    time_utc: str
    frequency_edges_hz: ArrayLike
    band_amplitudes_m: ArrayLike
    direction_moments: ArrayLike
    provenance: Mapping[str, str]

    def __post_init__(self):
        amplitudes = array(self.band_amplitudes_m, (None,), 'band_amplitudes_m')
        edges = array(self.frequency_edges_hz, (len(amplitudes)+1,), 'frequency_edges_hz')
        moments = array(self.direction_moments, (len(amplitudes), 4), 'direction_moments')
        if (not len(amplitudes) or np.any(amplitudes < 0) or not np.any(amplitudes > 0)
                or np.any(edges <= 0) or np.any(np.diff(edges) <= 0) or np.any(abs(moments) > 1)):
            raise ValueError('positive ordered bands, nonnegative nonempty energy and bounded moments required')
        if not isinstance(self.id, str) or not self.id or not isinstance(self.time_utc, str) or not self.time_utc.endswith('Z'):
            raise ValueError('observation identity and UTC timestamp required')
        try:
            stamp = datetime.fromisoformat(self.time_utc.replace('Z', '+00:00'))
            if stamp.utcoffset() != timedelta(0):
                raise ValueError('expected UTC')
        except ValueError as error:
            raise ValueError('valid observation UTC timestamp required') from error
        if not isinstance(self.provenance, Mapping) or not all(
                isinstance(k, str) and isinstance(v, str) for k, v in self.provenance.items()):
            raise ValueError('provenance must map strings to strings')
        for name, value in (('band_amplitudes_m', amplitudes), ('frequency_edges_hz', edges),
                            ('direction_moments', moments)):
            object.__setattr__(self, name, value)
        object.__setattr__(self, 'provenance', MappingProxyType(dict(self.provenance)))

    @property
    def variance_m2(self):
        return float(np.sum(self.band_amplitudes_m**2)/2)

    @property
    def hm0_m(self):
        return 4*np.sqrt(self.variance_m2)

    @property
    def tm02_s(self):
        edges = self.frequency_edges_hz
        mean_f2 = (edges[1:]**3-edges[:-1]**3)/(3*np.diff(edges))
        return float(np.sqrt(self.variance_m2/np.dot(self.band_amplitudes_m**2/2, mean_f2)))

    @property
    def peak_band_period_s(self):
        index = np.argmax(self.band_amplitudes_m**2/np.diff(self.frequency_edges_hz))
        return float(2/np.sum(self.frequency_edges_hz[index:index+2]))


def _catalog_data():
    text = files('topwater').joinpath('data/socal_sea_states.json').read_text(encoding='utf-8')
    data = json.loads(text)
    if data['version'] != 'topwater-observed-seas/1':
        raise ValueError('unsupported observed sea catalog')
    return data


def _observation(data, record):
    provenance = dict(station_id=data['station']['id'], station_name=data['station']['name'],
        source_library_sha256=data['source_library_sha256'],
        source_manifest_sha256=data['source_manifest_sha256'],
        source_files=json.dumps(record['source_files'], separators=(',', ':')),
        directional_measurement='first two Fourier harmonics, quantized by source export to int16/32767',
        coverage='usable measured frequency bands only; unmeasured tails unknown')
    return ObservedSeaState(record['id'], record['time_utc'], data['frequency_edges_hz'],
        record['band_amplitudes_m'], np.asarray(record['direction_moments_int16'])/32767, provenance)


def catalog():
    """Nine selected observed SoCal spectra, not a weighted population sample."""
    data = _catalog_data()
    return tuple(_observation(data, record) for record in data['records'])


def observed(case_id=None):
    """Return an observed record; default is the jointly selected SoCal median."""
    data = _catalog_data()
    case_id = data['baseline_id'] if case_id is None else case_id
    for record in data['records']:
        if record['id'] == case_id:
            return _observation(data, record)
    raise ValueError('unknown observed sea-state id: '+str(case_id))


def _envelope(frequency_hz, beam_m, gravity):
    if (isinstance(beam_m, bool) or isinstance(gravity, bool) or not np.isfinite(beam_m) or beam_m <= 0
            or not np.isfinite(gravity) or gravity <= 0):
        raise ValueError('positive finite beam_m and gravity required')
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    if not np.all(np.isfinite(frequency_hz)) or np.any(frequency_hz <= 0):
        raise ValueError('positive finite frequency support required')
    limit = min(.6, np.sqrt(gravity/(4*np.pi*beam_m)))
    # Equality is admitted; this allowance covers only floating-point roundoff.
    if np.any(frequency_hz > limit*(1+8*np.finfo(float).eps)):
        raise ValueError('observed spectrum exceeds f<=0.6 Hz or wavelength>=2*beam; no bands were removed')



def validate_envelope(waves, *, beam_m, gravity=9.81):
    """Reject incident components outside f<=0.6 Hz and wavelength>=2*beam.

    Zero-amplitude components carry no incident energy. Observed band-edge
    support, when recorded, is checked as well as the realized frequencies.
    """
    if not isinstance(waves, WaveComponents):
        raise TypeError('waves must be WaveComponents')
    _envelope(waves.angular_frequencies[waves.amplitudes > 0]/(2*np.pi), beam_m, gravity)
    support = waves.provenance.get('observed_active_frequency_max_hz')
    if support is not None:
        _envelope(float(support), beam_m, gravity)

def _directional_weights(moment):
    """Positive MEM on 720 nautical FROM bearings; inherited marine-gyro method.

    Uniform mixing regularizes inadmissible moment sequences. The mixing scale
    is exported; this is a reconstruction, not a uniquely observed distribution.
    """
    a1, b1, a2, b2 = moment
    c1, c2 = complex(a1, b1), complex(a2, b2)

    def parameters(scale):
        u, v = scale*c1, scale*c2
        denominator = 1-abs(u)**2
        if denominator <= 1e-6:
            return None
        p1 = (u-v*u.conjugate())/denominator
        p2 = v-u*p1
        sigma = 1-(p1*u.conjugate()+p2*v.conjugate()).real
        return (p1, p2) if sigma > 1e-6 else None

    scale, result = 1., parameters(1.)
    if result is None:
        left, right = 0., 1.
        for _ in range(40):
            middle = (left+right)/2
            if parameters(middle) is None:
                right = middle
            else:
                left = middle
        scale, result = left, parameters(left)
    p1, p2 = result
    theta = np.arange(720)*np.pi/360
    weights = 1/abs(1-p1*np.exp(-1j*theta)-p2*np.exp(-2j*theta))**2
    return weights/weights.sum(), scale


def realize(case=None, *, beam_m, seed=None, components_per_band=8, gravity=9.81):
    """Generate a DEEP-WATER realization of one WHOLE observed spectrum.

    Omitted seed obtains and records a fresh uint32 seed. All frequency, direction
    and phase draws use PCG64 from that seed. Band energies are unchanged for every
    seed. Within-band frequencies are stratified uniform; directional quantiles
    are stratified from a positive MEM reconstruction. A finite realization only
    approximates that directional PDF. Reject incompatible support before drawing.
    Finite-run elevation RMS fluctuates through interference; never normalize it.
    Breaking/slamming are separate physical validity checks, not certified here.
    """
    state = case if isinstance(case, ObservedSeaState) else observed(case)
    if seed is None:
        seed = secrets.randbits(32)
    if isinstance(seed, bool) or not isinstance(seed, Integral) or not 0 <= seed < 2**32:
        raise ValueError('seed must be a uint32 integer or None')
    if (isinstance(components_per_band, bool) or not isinstance(components_per_band, Integral)
            or components_per_band <= 0):
        raise ValueError('components_per_band must be a positive integer')
    active = state.band_amplitudes_m > 0
    maximum = float(state.frequency_edges_hz[1:][active].max())
    _envelope(maximum, beam_m, gravity)
    rng = np.random.Generator(np.random.PCG64(int(seed)))
    count = int(components_per_band)
    frequency, direction, phase, amplitude = [], [], [], []
    scales = []
    for left, right, a, moment in zip(state.frequency_edges_hz[:-1], state.frequency_edges_hz[1:],
                                    state.band_amplitudes_m, state.direction_moments):
        weights, scale = _directional_weights(moment)
        scales.append(scale)
        f = left+(np.arange(count)+rng.random(count))*(right-left)/count
        quantiles = (rng.permutation(count)+rng.random())/count
        cdf = np.cumsum(weights)
        cdf[-1] = 1.
        bearings = np.searchsorted(cdf, quantiles, side='right')*np.pi/360
        frequency.extend(f)
        direction.extend((1.5*np.pi-bearings) % (2*np.pi))
        phase.extend(rng.uniform(0, 2*np.pi, count))
        amplitude.extend(np.full(count, a/np.sqrt(count)))
    provenance = dict(state.provenance, observation_id=state.id, source_time_utc=state.time_utc,
        generator='topwater-observed-sea-realization/1; numpy PCG64',
        phase_distribution='independent uniform [0,2*pi); synthetic, not observed phases',
        direction_model='positive MEM, 0.5-degree FROM grid; finite stratified quantiles',
        regularization_min_scale=repr(float(min(scales))),
        regularization_band_scales=json.dumps(scales, separators=(',', ':')), components_per_band=str(count),
        measured_variance_m2=repr(state.variance_m2), beam_m=repr(float(beam_m)),
        gravity_m_s2=repr(float(gravity)), depth='infinite', observed_active_frequency_max_hz=repr(maximum))
    return WaveComponents(amplitude, 2*np.pi*np.asarray(frequency), direction, phase,
                          seed=int(seed), provenance=provenance)


def dumps_realization(waves):
    """Save actual arrays: replay never depends on regenerating a PRNG sequence."""
    return json.dumps(dict(version='topwater-sea-realization/1',
        waves={name: getattr(waves, name).tolist() for name in
               ('amplitudes', 'angular_frequencies', 'directions', 'phases')},
        seed=waves.seed, provenance=dict(waves.provenance)), allow_nan=False, indent=2)+'\n'


def loads_realization(text, *, beam_m=None):
    """Replay explicit arrays, checking variance closure and current hull envelope."""
    try:
        data = json.loads(text)
        if set(data) != {'version', 'waves', 'seed', 'provenance'} or data['version'] != 'topwater-sea-realization/1':
            raise ValueError('unsupported sea realization')
        if set(data['waves']) != {'amplitudes', 'angular_frequencies', 'directions', 'phases'}:
            raise ValueError('expected four explicit wave arrays')
        waves = WaveComponents(**data['waves'], seed=data['seed'], provenance=data['provenance'])
        source = waves.provenance
        if source['depth'] != 'infinite':
            raise ValueError('only deep-water realizations are supported')
        beam_m = float(source['beam_m']) if beam_m is None else beam_m
        gravity = float(source['gravity_m_s2'])
        _envelope(float(source['observed_active_frequency_max_hz']), beam_m, gravity)
        _envelope(waves.angular_frequencies[waves.amplitudes > 0]/(2*np.pi), beam_m, gravity)
        expected = float(source['measured_variance_m2'])
        if (not np.isfinite(expected) or expected <= 0 or not np.isclose(
                np.sum(waves.amplitudes**2)/2, expected, rtol=8e-15, atol=0)):
            raise ValueError('explicit realization no longer closes measured variance')
        return waves
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise ValueError('invalid sea realization: '+str(error)) from error
