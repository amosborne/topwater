"""Calculation controls affect clocks, packet fitting and numerical operators."""
from dataclasses import replace

import numpy as np
import pytest

from topwater.solve import SolveOptions
from topwater.solve.sampling import surface_sample_times
from topwater.solve.native.radiation import radiation_kernel
from topwater.solve.packets import minimum_cutoff_packet


def test_selected_timestep_and_packet_tolerance_change_the_constructed_packet():
    loose = minimum_cutoff_packet(.43, .7, 1., 1., 4., time_step_s=.02, error_tolerance_m=.03)
    tight = minimum_cutoff_packet(.43, .7, 1., 1., 4., time_step_s=.02, error_tolerance_m=.002)
    assert len(tight.target.time_s) == 400
    np.testing.assert_allclose(np.diff(tight.target.time_s), .02, rtol=0., atol=1e-15)
    assert tight.positive_components > loose.positive_components
    assert tight.diagnostics['certified_envelope_exceedance_bound_m'] <= .002


@pytest.mark.parametrize('step,interval', [(.01, .1), (.02, .2), (.03, .3)])
def test_default_cadence_selects_native_samples_without_restricting_duration(step, interval):
    options = SolveOptions(time_step_s=step)
    assert options.surface_interval_s == pytest.approx(interval)
    time = np.arange(28)*step
    saved = surface_sample_times(time, every=options.render_every)
    assert np.isin(saved, time).all()
    assert saved[-1] == time[-1]


@pytest.mark.parametrize('kwargs', [dict(render_every=0), dict(render_every=1.5), dict(render_every=None), dict(render_every=True), dict(time_step_s=0.),
    dict(packet_tolerance_m=0.), dict(surface_padding_beams=-1.),
    dict(cpu_workers=0), dict(batch_points=0)])
def test_invalid_calculation_controls_fail_before_work(kwargs):
    with pytest.raises(ValueError):
        SolveOptions(**kwargs)


def test_spectral_controls_change_coefficient_and_field_frequency_banks():
    options = SolveOptions(coefficient_step_rad_s=.2,
        radiation_bandwidth_rad_s=4., batch_points=16)
    radiation = options.radiation_frequencies(2.1)
    assert radiation[-1] == pytest.approx(4.)
    np.testing.assert_allclose(np.diff(radiation), .2)
    from topwater.solve.fourier import surface_radiation_frequencies
    settings = options.radiation_settings(128)
    assert settings.chunk_points == 16
    np.testing.assert_array_equal(surface_radiation_frequencies(np.arange(61)*.01, settings), 2*np.pi*np.fft.rfftfreq(61,.01))




def test_result_metadata_roundtrip_and_partial_video_endpoint(tmp_path):
    from .test_result_render import saved_result
    from topwater.workflow_data import Result
    from topwater.render.geometry import render_clock
    original = saved_result()
    source = replace(original, surface_interval_s=.1, view_bounds_xy_m=[[-2., -2.], [2., 2.]], video_size_px=[642, 642])
    destination = tmp_path/'result.npz'
    source.save(destination)
    restored = Result.load(destination)
    np.testing.assert_array_equal(restored.view_bounds_xy_m, source.view_bounds_xy_m)
    assert restored.surface_interval_s == .1
    np.testing.assert_array_equal(restored.video_size_px, [642, 642])
    indices, fps = render_clock(restored)
    np.testing.assert_array_equal(indices, [0])
    assert fps == 10.
    # The off-cadence endpoint is preserved in the numeric result.
    assert restored.surface_time_s[-1] == original.surface_time_s[-1] == .08


def test_real_movie_uses_solver_view_and_omits_off_cadence_endpoint(tmp_path):
    from .test_result_render import saved_result
    from topwater.render.workflow import render_result
    from topwater.render.scene import RenderSettings
    source = replace(saved_result(), surface_interval_s=.1,
        view_bounds_xy_m=[[-2., -2.], [2., 2.]])
    summary = render_result(source, tmp_path/'video.mp4', settings=RenderSettings(framing=.5))
    assert summary['frames'] == 1 and summary['fps'] == 10.
    assert summary['source_end_s'] == 0.
    assert source.surface_time_s[-1] == .08


@pytest.mark.parametrize('changes', [dict(time_step_s=.02)])
def test_wecsim_rejects_unsupported_motion_controls_before_discovery(monkeypatch, changes):
    from topwater.solve.wecsim import engine as wecsim_backend
    monkeypatch.setattr(wecsim_backend, 'preflight', lambda: pytest.fail('MATLAB discovery must not run'))
    with pytest.raises(ValueError, match='requires the native solver'):
        from topwater.solve import simulate
        simulate(None, None, 100., solver='wecsim', options=SolveOptions(**changes))


@pytest.mark.parametrize('every,count', [(1, 28), (3, 10), (7, 5), (100, 2)])
def test_render_stride_controls_saved_samples_and_keeps_endpoint(every, count):
    options = SolveOptions(time_step_s=.02, render_every=every)
    time = np.arange(28)*options.time_step_s
    saved = surface_sample_times(time, every=options.render_every)
    assert len(saved) == count
    assert saved[-1] == time[-1]
    np.testing.assert_array_equal(saved[:-1], time[::every][:len(saved)-1])
    assert options.surface_interval_s == pytest.approx(.02*every)



def test_renderer_uses_archive_resolution_instead_of_catalogue_default(tmp_path):
    from .test_result_render import saved_result
    from topwater.render.workflow import render_result
    source = replace(saved_result(), surface_interval_s=.1,
        view_bounds_xy_m=[[-2., -2.], [2., 2.]], video_size_px=[642, 642])
    summary = render_result(source, tmp_path/'video.mp4')
    assert summary['size'] == (642, 642)
    assert summary['frames'] == 1


@pytest.mark.parametrize('size', [[0, 800], [801, 800], [800, -2]])
def test_archive_rejects_invalid_video_dimensions(size):
    from .test_result_render import saved_result
    with pytest.raises(ValueError, match='positive even'):
        replace(saved_result(), video_size_px=size)
