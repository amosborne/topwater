"""Public directory defaults, deterministic generation and failure preservation."""
from pathlib import Path

import numpy as np
import pytest

from topwater.cli import main
from topwater.cli import parser
from topwater.workflow_data import Sea


@pytest.mark.parametrize('argv,directory,dwell', [
    (['solve'], Path('.'),100.),
    (['solve','--dwell','60'], Path('.'),60.),
    (['solve','example with spaces','--dwell','60','--solver=wecsim'],Path('example with spaces'),60.),
])
def test_named_dwell_and_optional_directory_have_independent_defaults(argv,directory,dwell):
    parsed=parser().parse_args(argv)
    assert parsed.directory==directory and parsed.dwell==dwell
    assert parsed.solver == ('wecsim' if '--solver=wecsim' in argv else 'native')


def test_native_solver_is_the_default_and_explicit_choice():
    assert parser().parse_args(['solve', '--dwell', '60']).solver == 'native'
    assert parser().parse_args(['solve', '--dwell', '60', '--solver=native']).solver == 'native'


def test_generate_sea_defaults_to_cwd_and_preserves_exact_band_energy(tmp_path,monkeypatch):
    from topwater.sea import states as sea_states
    monkeypatch.chdir(tmp_path)
    observed=sea_states.observed()
    arguments=['sea','--record',observed.id,'--seed','42']
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
    assert main(['sea',str(tmp_path),'--record','unknown'])==2
    assert destination.read_bytes()==b'previous sea'
    assert 'unknown observed' in capsys.readouterr().err
    assert 'unknown observed' in (tmp_path/'sea.log').read_text()
    assert 'status=failed' in (tmp_path/'sea.log').read_text()


@pytest.mark.parametrize('duration',['60','nan','-1'])
def test_missing_inputs_or_bad_duration_never_replaces_result(tmp_path,duration,capsys):
    (tmp_path/'result.npz').write_bytes(b'previous result')
    (tmp_path/'video.mp4').write_bytes(b'previous video')
    (tmp_path/'result.png').write_bytes(b'previous plot')
    assert main(['solve',str(tmp_path),'--dwell',duration])==2
    assert (tmp_path/'result.npz').read_bytes()==b'previous result'
    assert (tmp_path/'video.mp4').read_bytes()==b'previous video'
    assert (tmp_path/'result.png').read_bytes()==b'previous plot'
    assert capsys.readouterr().err.startswith('topwater solve:')
    assert 'status=failed' in (tmp_path/'result.log').read_text()
    assert 'output_published=false' in (tmp_path/'result.log').read_text()


def test_stronger_observation_is_packaged_without_changing_legacy_catalog():
    from topwater.sea import states as sea_states
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
    from topwater.solve import outputs as solve_outputs
    from topwater.solve.wecsim import engine as wecsim_backend
    from topwater import solve as native
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
    def ready_motion(body, waves, dwell_s, *, rest_s, ramp_s, **device_options):
        calls.append((dwell_s, rest_s, ramp_s))
        return solved
    monkeypatch.setattr(native, 'simulate', ready_motion)
    monkeypatch.setattr(solve_outputs, 'save_solve_outputs', capture_plot_sea)
    assert main(['solve', str(tmp_path), *flags]) == 0
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
    assert 'BUDGET' not in report
    assert 'VIDEO frames=2 fps=12.5 size=800x800' in report
    assert 'source_end_s=0.080 encoded_duration_s=0.160 verified=true' in report
    assert 'TIMING render preparation_s=' in report
    assert 'verification_s=' in report
    assert 'timings={' not in report


@pytest.mark.parametrize('argv,rest,ramp', [
    (['solve', '--solver=wecsim'], 10., 5.),
    (['solve', 'example', '--dwell', '100', '--rest', '2', '--ramp', '3', '--solver=wecsim'], 2., 3.),
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
    ('100', ['--rest', '1.001'], 'rest must align to the selected 0.01s'),
    ('100', ['--ramp', '5.005'], 'ramp must align to the selected 0.01s'),
    ('100.001', [], 'dwell must align to the selected 0.01s'),
])
def test_bad_packet_timing_is_fatal_before_backend_or_publication(tmp_path, monkeypatch, capsys,
                                                               dwell, flags, message):
    from topwater import solve as solve_api
    monkeypatch.setattr(solve_api, 'simulate', lambda *args, **kwargs: pytest.fail('backend must not run'))
    preserved = {name: name.encode() for name in ('boat.npz', 'sea.npz', 'result.npz', 'result.png', 'video.mp4')}
    for name, data in preserved.items():
        (tmp_path/name).write_bytes(data)
    assert main(['solve', str(tmp_path), '--dwell', dwell, '--solver=wecsim', *flags]) == 2
    assert message in capsys.readouterr().err
    for name, data in preserved.items():
        assert (tmp_path/name).read_bytes() == data
    log = (tmp_path/'result.log').read_text()
    assert 'END status=failed' in log and 'output_published=false' in log
    assert 'PUBLISHED ' not in log


def test_solve_cli_forwards_meaningful_controls_and_accepts_partial_video_interval(tmp_path, monkeypatch):
    from topwater import solve as native
    from topwater.solve import outputs as solve_outputs
    from .test_workflow_data import boat, result
    from topwater.workflow_data import Sea
    boat().save(tmp_path/'boat.npz')
    Sea([.1], [1.], [0.], [0.]).save(tmp_path/'sea.npz')
    selected = []
    def solve(*args, **kwargs):
        selected.append(kwargs['options'])
        return result()
    monkeypatch.setattr(native, 'simulate', solve)
    monkeypatch.setattr(solve_outputs, 'save_solve_outputs',
        lambda *args, **kw: (tmp_path/'result.npz', tmp_path/'result.png'))
    assert main(['solve', str(tmp_path), '--dwell', '100.02', '--time-step', '.02',
        '--packet-tolerance', '.005']) == 0
    from topwater.solve import SolveOptions
    assert selected == [SolveOptions(time_step_s=.02, packet_tolerance_m=.005)]
    assert 'time_step_s=0.02' in (tmp_path/'result.log').read_text()


def test_development_controls_are_not_exposed_as_cli_options():
    options = parser().parse_args(['solve'])
    assert options.dwell == 100.
    assert not hasattr(options, 'surface_spacing_m')
    assert not hasattr(options, 'fps') and not hasattr(options, 'render_every')
    assert not hasattr(options, 'surface_padding_beams')
    assert not hasattr(options, 'radiation_memory_s')


def test_random_seed_is_logged_without_console_output_and_can_replay(tmp_path, capsys):
    import re
    from topwater.sea import states as sea_states
    record = sea_states.observed().id
    assert main(['sea', str(tmp_path), '--record', record]) == 0
    first = Sea.load(tmp_path/'sea.npz')
    console = capsys.readouterr()
    assert 'seed=' not in console.out+console.err
    log = (tmp_path/'sea.log').read_text()
    seed = re.search(r'SEA record=.* seed=(\d+)', log).group(1)
    assert main(['sea', str(tmp_path), '--record', record, '--seed', seed]) == 0
    second = Sea.load(tmp_path/'sea.npz')
    for member in Sea.members:
        np.testing.assert_array_equal(getattr(first, member), getattr(second, member))
    console = capsys.readouterr()
    assert 'seed=' not in console.out+console.err


def test_gpu_name_is_a_direct_device_argument():
    selected = parser().parse_args(['solve', '--device', 'My GPU'])
    assert selected.device == 'My GPU'
    assert not hasattr(selected, 'gpu_device')



@pytest.mark.parametrize('flag', ['--radiation-memory', '--surface-spacing'])
def test_automatic_history_and_spacing_have_no_cli_override(flag):
    with pytest.raises(SystemExit):
        parser().parse_args(['solve', flag, '10'])
