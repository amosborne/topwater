"""Measured-spectrum conservation and replay controls on actual SoCal data."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from topwater import sea_states


def test_baseline_matches_independent_native_noaa_energy_and_direction_row():
    state = sea_states.observed()
    native = json.loads((Path(__file__).parent/'data/socal_baseline_native.json').read_text())
    assert state.id == native['id'] == 'ndbc-46086-20250705T2040Z'
    assert state.time_utc == native['time_utc'] == '2025-07-05T20:40:00Z'
    f, width, density = [np.asarray(native[name]) for name in
                         ('frequency_hz', 'bandwidth_hz', 'density_m2_per_hz')]
    edges = state.frequency_edges_hz
    overlap = np.maximum(0, np.minimum((f+width/2)[:, None], edges[None, 1:])
                         -np.maximum((f-width/2)[:, None], edges[None, :-1]))
    band_variance = density@overlap
    assert_allclose(state.band_amplitudes_m**2/2, band_variance, rtol=2e-14, atol=1e-17)
    # Native frequency bands are partitioned completely; rebinning cannot change m0.
    assert_allclose(state.variance_m2, density@width, rtol=2e-14)
    numerator = np.einsum('n,nk,nb->bk', density, native['direction_moments'], overlap)
    moments = np.divide(numerator, band_variance[:, None], out=np.zeros_like(numerator),
                        where=band_variance[:, None] > 0)
    assert_allclose(state.direction_moments, moments, atol=1.5/32767, rtol=0)
    assert_allclose(state.hm0_m, 1.3922643427165697, rtol=2e-14)
    assert_allclose(state.tm02_s, 6.940525365844814, rtol=2e-14)


def test_catalog_is_nine_distinct_joint_observations_with_readonly_owned_arrays():
    states = sea_states.catalog()
    assert len(states) == len({state.id for state in states}) == len({state.time_utc for state in states}) == 9
    for state in states:
        assert state.provenance['station_id'] == '46086'
        assert len(json.loads(state.provenance['source_files'])) == 5
        assert state.direction_moments.shape == (48, 4)
        for values in (state.frequency_edges_hz, state.band_amplitudes_m, state.direction_moments):
            assert not values.flags.writeable
    state = states[0]
    original = state.band_amplitudes_m.copy()
    replacement = replace(state, band_amplitudes_m=original)
    original[:] = 0
    assert np.any(replacement.band_amplitudes_m > 0)
    with pytest.raises(ValueError, match='UTC timestamp'):
        replace(state, time_utc='garbageZ')


def test_realization_retains_each_observed_bands_variance_for_every_seed():
    state = sea_states.observed()
    for seed in (0, 42, 2**32-1):
        waves = sea_states.realize(state, beam_m=9*.3048, seed=seed)
        amplitudes = waves.amplitudes.reshape(48, 8)
        frequency = waves.angular_frequencies.reshape(48, 8)/(2*np.pi)
        assert_allclose(np.sum(amplitudes**2, axis=1), state.band_amplitudes_m**2, rtol=3e-15)
        assert np.all(frequency >= state.frequency_edges_hz[:-1, None])
        assert np.all(frequency < state.frequency_edges_hz[1:, None])
        assert np.all((waves.phases >= 0)&(waves.phases < 2*np.pi))
        assert waves.seed == seed
        assert waves.provenance['source_time_utc'] == state.time_utc
        assert waves.provenance['depth'] == 'infinite'


def test_direction_convention_and_stratified_realization_retain_observed_harmonics():
    state = sea_states.observed()
    count = 256
    waves = sea_states.realize(state, beam_m=9*.3048, seed=42, components_per_band=count)
    # Convert east/CCW PROPAGATION back to nautical FROM; test against all four
    # measured harmonics, rather than a plot or a mean heading alone.
    theta = (1.5*np.pi-waves.directions.reshape(48, count)) % (2*np.pi)
    samples = np.stack((np.cos(theta), np.sin(theta), np.cos(2*theta), np.sin(2*theta)), axis=-1).mean(axis=1)
    scales = np.asarray(json.loads(waves.provenance['regularization_band_scales']))
    assert scales.shape == (48,)
    assert np.all((scales > 0)&(scales <= 1))
    assert_allclose(samples, state.direction_moments*scales[:, None], rtol=0, atol=.04)
    # Bearings belong to the documented 0.5-degree FROM grid.
    assert_allclose(theta/(np.pi/360), np.rint(theta/(np.pi/360)), atol=2e-12)


def test_optional_seed_records_fresh_seed_and_explicit_arrays_replay(monkeypatch):
    monkeypatch.setattr(sea_states.secrets, 'randbits', lambda bits: 123456)
    generated = sea_states.realize(beam_m=9*.3048)
    repeated = sea_states.realize(beam_m=9*.3048, seed=generated.seed)
    saved = sea_states.dumps_realization(generated)
    replayed = sea_states.loads_realization(saved)
    assert generated.seed == 123456
    for name in ('amplitudes', 'angular_frequencies', 'directions', 'phases'):
        assert_array_equal(getattr(generated, name), getattr(repeated, name))
        assert_array_equal(getattr(generated, name), getattr(replayed, name))
    # Exact replay includes the traveling-wave phase at a location away from the
    # origin, using WaveComponents' k dot x - omega*t + phase convention.
    def elevation(waves):
        k = waves.angular_frequencies**2/9.81
        return np.sum(waves.amplitudes*np.cos(k*(3*np.cos(waves.directions)-2*np.sin(waves.directions))
            -waves.angular_frequencies*7+waves.phases))
    assert elevation(generated) == elevation(replayed)
    changed = sea_states.realize(beam_m=9*.3048, seed=generated.seed+1)
    assert not np.array_equal(generated.phases, changed.phases)
    assert_array_equal(generated.amplitudes, changed.amplitudes)


def test_incompatible_support_rejected_without_changing_observation_or_replay():
    # This actual hour has positive measured energy through 0.495 Hz. The median
    # hour's zero tail ends earlier and need not be rejected for an 11 ft beam.
    state = sea_states.observed('ndbc-46086-20230704T0340Z')
    before = state.band_amplitudes_m.copy()
    with pytest.raises(ValueError, match='no bands were removed'):
        sea_states.realize(state, beam_m=11*.3048, seed=42)
    waves = sea_states.realize(state, beam_m=9*.3048, seed=42)
    with pytest.raises(ValueError, match='no bands were removed'):
        sea_states.loads_realization(sea_states.dumps_realization(waves), beam_m=11*.3048)
    assert_array_equal(state.band_amplitudes_m, before)
    equality_beam = 9.81/(4*np.pi*state.frequency_edges_hz[-1]**2)
    sea_states.realize(state, beam_m=equality_beam, seed=42)
    corrupt = json.loads(sea_states.dumps_realization(waves))
    corrupt['waves']['amplitudes'][int(np.argmax(waves.amplitudes))] *= 2
    with pytest.raises(ValueError, match='variance'):
        sea_states.loads_realization(json.dumps(corrupt))


def test_frequency_ceiling_rejects_full_positive_support_even_for_small_beam():
    state = sea_states.observed()
    outside = replace(state, frequency_edges_hz=state.frequency_edges_hz*2)
    with pytest.raises(ValueError, match='no bands were removed'):
        sea_states.realize(outside, beam_m=.1, seed=42)


@pytest.mark.parametrize('seed', [True, -1, 2**32, 1.5])
def test_invalid_seed_is_not_coerced(seed):
    with pytest.raises(ValueError, match='uint32'):
        sea_states.realize(beam_m=9*.3048, seed=seed)


def test_input_wave_envelope_checks_actual_components_and_observed_support():
    from topwater import WaveComponents

    waves = WaveComponents([.1, 0.], 2*np.pi*np.array([.2, 1.]), [0., 0.], [0., 0.])
    sea_states.validate_envelope(waves, beam_m=3.)
    energetic_short = WaveComponents([.1], [2*np.pi*.6], [0.], [0.])
    with pytest.raises(ValueError, match='wavelength'):
        sea_states.validate_envelope(energetic_short, beam_m=3.)
    observed = WaveComponents([.1], [2*np.pi*.2], [0.], [0.],
        provenance={'observed_active_frequency_max_hz': '.6'})
    with pytest.raises(ValueError, match='no bands were removed'):
        sea_states.validate_envelope(observed, beam_m=3.)
