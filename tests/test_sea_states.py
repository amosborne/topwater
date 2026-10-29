"""Measured-spectrum conservation and replay controls on actual SoCal data."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest

from topwater.sea import states as sea_states


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
        waves = sea_states.realize(state, seed=seed)
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
    waves = sea_states.realize(state, seed=42, components_per_band=count)
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
    generated = sea_states.realize()
    repeated = sea_states.realize(seed=generated.seed)
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
    changed = sea_states.realize(seed=generated.seed+1)
    assert not np.array_equal(generated.phases, changed.phases)
    assert_array_equal(generated.amplitudes, changed.amplitudes)


def test_realization_replay_retains_observation_and_checks_variance():
    state = sea_states.observed('ndbc-46086-20230704T0340Z')
    before = state.band_amplitudes_m.copy()
    waves = sea_states.realize(state, seed=42)
    replayed = sea_states.loads_realization(sea_states.dumps_realization(waves))
    assert_array_equal(state.band_amplitudes_m, before)
    for name in ('amplitudes', 'angular_frequencies', 'directions', 'phases'):
        assert_array_equal(getattr(waves, name), getattr(replayed, name))
    corrupt = json.loads(sea_states.dumps_realization(waves))
    corrupt['waves']['amplitudes'][int(np.argmax(waves.amplitudes))] *= 2
    with pytest.raises(ValueError, match='variance'):
        sea_states.loads_realization(json.dumps(corrupt))


def test_full_positive_high_frequency_support_is_not_clipped():
    state = sea_states.observed()
    high = replace(state, frequency_edges_hz=state.frequency_edges_hz*3)
    original = sea_states.realize(state, seed=42)
    waves = sea_states.realize(high, seed=42)
    repeated = sea_states.realize(high, seed=42)
    replayed = sea_states.loads_realization(sea_states.dumps_realization(waves))
    assert np.max(waves.angular_frequencies[waves.amplitudes > 0]/(2*np.pi)) > .6
    assert_allclose(waves.angular_frequencies, original.angular_frequencies*3, rtol=1e-15)
    assert_allclose(np.sum(waves.amplitudes**2)/2, high.variance_m2, rtol=8e-15)
    for name in ('amplitudes', 'directions', 'phases'):
        assert_array_equal(getattr(waves, name), getattr(original, name))
    for name in ('amplitudes', 'angular_frequencies', 'directions', 'phases'):
        assert_array_equal(getattr(waves, name), getattr(repeated, name))
        assert_array_equal(getattr(waves, name), getattr(replayed, name))


@pytest.mark.parametrize('seed', [True, -1, 2**32, 1.5])
def test_invalid_seed_is_not_coerced(seed):
    with pytest.raises(ValueError, match='uint32'):
        sea_states.realize(seed=seed)


def test_solve_preserves_components_beyond_old_domain_limits(capsys):
    from topwater import WaveComponents
    from topwater.solve.preparation import prepare_case
    from .reference_boat import reference_body

    waves = WaveComponents([.1, .02, 0.], 2*np.pi*np.array([.8, 1.2, 12.]),
        [.2, .3, .4], [.1, .2, .3],
        provenance={'observed_active_frequency_max_hz': '4.0'})
    setup = prepare_case(reference_body(), waves, 2., rest_s=1., ramp_s=1.)
    assert_array_equal(setup.native_waves.amplitudes, waves.amplitudes)
    assert_array_equal(setup.native_waves.angular_frequencies, waves.angular_frequencies)
    assert_array_equal(setup.train.carriers.angular_frequencies, waves.angular_frequencies)
    assert setup.train.omega[-1]/(2*np.pi) >= 1.2
    assert 'packet-sideband-operating-envelope' not in capsys.readouterr().err
