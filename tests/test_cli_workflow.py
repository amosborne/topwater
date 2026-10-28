"""Public directory defaults, deterministic generation and failure preservation."""
from pathlib import Path

import numpy as np
import pytest

from topwater.cli import main, parser
from topwater.workflow_data import Sea


@pytest.mark.parametrize('argv,directory,dwell', [
    (['solve','60'], Path('.'),60.),
    (['solve','example with spaces','60','--solver=wecsim'],Path('example with spaces'),60.),
])
def test_optional_directory_does_not_consume_required_dwell(argv,directory,dwell):
    parsed=parser().parse_args(argv)
    assert parsed.directory==directory and parsed.dwell==dwell
    assert parsed.solver == 'wecsim'


def test_retired_native_solver_is_not_a_cli_choice():
    with pytest.raises(SystemExit) as error:
        parser().parse_args(['solve', '60', '--solver=native'])
    assert error.value.code == 2


def test_generate_sea_defaults_to_cwd_and_preserves_exact_band_energy(tmp_path,monkeypatch):
    from topwater import sea_states
    monkeypatch.chdir(tmp_path)
    observed=sea_states.observed()
    arguments=['generate-sea','--record',observed.id,'--seed','42']
    assert main(arguments)==0
    first=Sea.load(tmp_path/'sea.npz')
    np.testing.assert_allclose(np.sum(first.amplitudes_m.reshape(48,8)**2,axis=1),
                               observed.band_amplitudes_m**2,rtol=3e-15)
    assert main(arguments)==0
    second=Sea.load(tmp_path/'sea.npz')
    for name in Sea.members:
        np.testing.assert_array_equal(getattr(first,name),getattr(second,name))
    assert set(p.name for p in tmp_path.iterdir())=={'sea.npz','sea.log'}
    assert 'status=success' in (tmp_path/'sea.log').read_text()


def test_generation_failure_preserves_existing_sea(tmp_path,capsys):
    destination=tmp_path/'sea.npz'
    destination.write_bytes(b'previous sea')
    assert main(['generate-sea',str(tmp_path),'--record','unknown'])==2
    assert destination.read_bytes()==b'previous sea'
    assert 'unknown observed' in capsys.readouterr().err
    assert 'unknown observed' in (tmp_path/'sea.log').read_text()
    assert 'status=failed' in (tmp_path/'sea.log').read_text()


@pytest.mark.parametrize('duration',['60','nan','-1'])
def test_missing_inputs_or_bad_duration_never_replaces_result(tmp_path,duration,capsys):
    (tmp_path/'result.npz').write_bytes(b'previous result')
    (tmp_path/'video.mp4').write_bytes(b'previous video')
    (tmp_path/'result.png').write_bytes(b'previous plot')
    assert main(['solve',str(tmp_path),duration])==2
    assert (tmp_path/'result.npz').read_bytes()==b'previous result'
    assert (tmp_path/'video.mp4').read_bytes()==b'previous video'
    assert (tmp_path/'result.png').read_bytes()==b'previous plot'
    assert capsys.readouterr().err.startswith('topwater solve:')
    assert 'status=failed' in (tmp_path/'result.log').read_text()
    assert 'output_published=false' in (tmp_path/'result.log').read_text()


def test_stronger_observation_is_packaged_without_changing_legacy_catalog():
    from topwater import sea_states
    assert len(sea_states.catalog())==9
    observed=sea_states.observed('ndbc-46086-20231223T0040Z')
    np.testing.assert_allclose(observed.hm0_m,1.6574679484080534,rtol=1e-14)
    np.testing.assert_allclose(observed.tm02_s,7.23656714864799,rtol=1e-14)


def test_render_missing_result_is_clear_failure(tmp_path,capsys):
    assert main(['render',str(tmp_path)])==2
    assert 'topwater render:' in capsys.readouterr().err


def test_corrupt_result_is_clear_fatal_and_preserves_video(tmp_path,capsys):
    (tmp_path/'result.npz').write_bytes(b'not a ZIP archive')
    (tmp_path/'video.mp4').write_bytes(b'previous video')
    assert main(['render',str(tmp_path)])==2
    assert 'invalid or truncated Result NPZ' in capsys.readouterr().err
    assert (tmp_path/'video.mp4').read_bytes()==b'previous video'


@pytest.mark.parametrize('flags,timing,total,midpoint', [
    ([], (100., 10., 5.), 130., 65.),
    (['--solver=wecsim', '--rest', '2', '--ramp', '3'], (100., 2., 3.), 110., 55.),
])
def test_wecsim_cli_publishes_result_and_motion_plot(tmp_path, monkeypatch, flags, timing, total, midpoint):
    from topwater import solve_outputs, wecsim_backend
    from .test_workflow_data import boat, result
    from topwater.workflow_data import Result
    solved = result()
    boat().save(tmp_path/'boat.npz')
    original_sea = Sea([.1, .2], [1., 2.], [.3, .6], [.4, .8])
    original_sea.save(tmp_path/'sea.npz')
    original_inputs = {name: (tmp_path/name).read_bytes() for name in ('boat.npz', 'sea.npz')}
    calls = []
    plotted_seas = []
    publish_pair = solve_outputs.save_solve_outputs
    def capture_plot_sea(data, directory, *, sea=None):
        plotted_seas.append(sea)
        return publish_pair(data, directory, sea=sea)
    def ready_motion(body, waves, dwell_s, *, rest_s, ramp_s):
        calls.append((dwell_s, rest_s, ramp_s))
        return solved
    monkeypatch.setattr(wecsim_backend, 'solve', ready_motion)
    monkeypatch.setattr(solve_outputs, 'save_solve_outputs', capture_plot_sea)
    assert main(['solve', str(tmp_path), '100', *flags]) == 0
    assert calls == [timing]
    assert len(plotted_seas) == 1
    for name in Sea.members:
        np.testing.assert_array_equal(getattr(plotted_seas[0], name), getattr(original_sea, name))
    for name, original in original_inputs.items():
        assert (tmp_path/name).read_bytes() == original
    np.testing.assert_array_equal(Result.load(tmp_path/'result.npz').position_m, solved.position_m)
    with np.load(tmp_path/'result.npz', allow_pickle=False) as archive:
        assert set(archive.files) == set(Result.members)
        assert not any(name.startswith('sea_') for name in archive.files)
    assert (tmp_path/'result.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    log = (tmp_path/'result.log').read_text()
    assert (f'PACKET dwell_s={timing[0]:g} rest_each_s={timing[1]:g} ramp_each_s={timing[2]:g} '
            f'total_s={total:g} dwell_midpoint_s={midpoint:g} time_step_s=0.01') in log
    assert 'DIRECTLY at dwell midpoint' in log
    assert 'phase zero means a target carrier crest' in log
    assert 'PUBLISHED ' + str(tmp_path/'result.npz') in log
    assert 'PUBLISHED ' + str(tmp_path/'result.png') in log
    assert 'END status=success' in log


def test_real_render_cli_records_validated_publication_in_video_log(tmp_path):
    from .test_result_render import saved_result
    saved_result().save(tmp_path/'result.npz')
    assert main(['render',str(tmp_path)])==0
    assert (tmp_path/'video.mp4').stat().st_size>0
    report=(tmp_path/'video.log').read_text()
    assert 'STAGE load self-contained result' in report
    assert 'status=success' in report and 'output_published=true' in report
    assert 'BUDGET render maximum_seconds=240 surface_frames=2' in report
    assert 'VIDEO frames=2 fps=12.5 size=800x800' in report
    assert 'source_end_s=0.080 encoded_duration_s=0.160 verified=true' in report
    assert 'TIMING render preparation_s=' in report
    assert 'verification_s=' in report
    assert 'timings={' not in report


@pytest.mark.parametrize('argv,rest,ramp', [
    (['solve', '100', '--solver=wecsim'], 10., 5.),
    (['solve', 'example', '100', '--rest', '2', '--ramp', '3', '--solver=wecsim'], 2., 3.),
])
def test_solve_exposes_dwell_and_explicit_per_side_timing(argv, rest, ramp):
    args = parser().parse_args(argv)
    assert args.dwell == 100.
    assert args.rest == rest and args.ramp == ramp


@pytest.mark.parametrize('dwell,flags,message', [
    ('nan', [], 'dwell must be positive and finite'),
    ('0', [], 'dwell must be positive and finite'),
    ('100', ['--rest', '-1'], 'rest must be nonnegative and finite'),
    ('100', ['--ramp', '0'], 'ramp must be positive and finite'),
    ('100', ['--rest', '1.001'], 'rest must align to the fixed 0.01s'),
    ('100', ['--ramp', '5.005'], 'ramp must align to the fixed 0.01s'),
    ('100.001', [], 'dwell must align to the fixed 0.01s'),
])
def test_bad_packet_timing_is_fatal_before_backend_or_publication(tmp_path, monkeypatch, capsys,
                                                               dwell, flags, message):
    from topwater import wecsim_backend
    monkeypatch.setattr(wecsim_backend, 'solve', lambda *args, **kwargs: pytest.fail('backend must not run'))
    preserved = {name: name.encode() for name in ('boat.npz', 'sea.npz', 'result.npz', 'result.png', 'video.mp4')}
    for name, data in preserved.items():
        (tmp_path/name).write_bytes(data)
    assert main(['solve', str(tmp_path), dwell, '--solver=wecsim', *flags]) == 2
    assert message in capsys.readouterr().err
    for name, data in preserved.items():
        assert (tmp_path/name).read_bytes() == data
    log = (tmp_path/'result.log').read_text()
    assert 'END status=failed' in log and 'output_published=false' in log
    assert 'PUBLISHED ' not in log


@pytest.mark.parametrize('motion_count,surface_count,budget', [
    (13001, 1626, 1626.),
    (13001, 2, 240.),
    (2, 751, 751.),
])
def test_render_budget_uses_saved_water_frames_not_motion_samples(tmp_path, monkeypatch,
                                                               motion_count, surface_count, budget):
    from types import SimpleNamespace
    from topwater import result_render
    from topwater.visualization import RenderSettings
    from topwater.workflow_data import Result
    saved = SimpleNamespace(time_s=np.linspace(0., 130., motion_count),
                            surface_time_s=np.linspace(0., 130., surface_count))
    monkeypatch.setattr(Result, 'load', classmethod(lambda cls, path: saved))
    calls = []
    def render(result, destination, settings=None):
        assert result is saved and destination == tmp_path/'video.mp4'
        calls.append(settings.maximum_seconds)
        destination.write_bytes(b'mocked validated video')
        timings = dict(frames=surface_count, fps=12.5, size=(800, 800), source_end_s=130.,
            video_duration_s=surface_count/12.5)
        timings.update({key:0. for key in ('preparation_seconds', 'setup_seconds', 'geometry_seconds',
            'draw_seconds', 'encoder_start_seconds', 'encoder_send_seconds', 'encoder_close_seconds',
            'verification_seconds', 'elapsed_s')})
        return timings
    monkeypatch.setattr(result_render, 'render_result', render)
    assert main(['render', str(tmp_path)]) == 0
    assert calls == [budget]
    assert RenderSettings().maximum_seconds == 240.
    log = (tmp_path/'video.log').read_text()
    assert f'BUDGET render maximum_seconds={budget:g} surface_frames={surface_count}' in log
    assert 'END status=success' in log and 'output_published=true' in log
