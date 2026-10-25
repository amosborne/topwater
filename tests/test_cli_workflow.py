"""Public directory defaults, deterministic generation and failure preservation."""
from pathlib import Path

import numpy as np
import pytest

from topwater.cli import main, parser
from topwater.workflow_data import Sea


@pytest.mark.parametrize('argv,directory,duration', [
    (['solve','60'], Path('.'),60.),
    (['solve','example with spaces','60','--solver=wecsim'],Path('example with spaces'),60.),
])
def test_optional_directory_does_not_consume_required_duration(argv,directory,duration):
    parsed=parser().parse_args(argv)
    assert parsed.directory==directory and parsed.duration==duration


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
def test_unavailable_native_or_bad_duration_never_replaces_result(tmp_path,duration,capsys):
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


def test_wecsim_cli_publishes_result_and_motion_plot(tmp_path, monkeypatch):
    from topwater import wecsim_backend
    from .test_workflow_data import boat, result
    from topwater.workflow_data import Result
    solved = result()
    boat().save(tmp_path/'boat.npz')
    Sea(np.array([]), np.array([]), np.array([]), np.array([])).save(tmp_path/'sea.npz')
    calls = []
    def ready_motion(body, waves, duration):
        calls.append(duration)
        return solved
    monkeypatch.setattr(wecsim_backend, 'solve', ready_motion)
    assert main(['solve', str(tmp_path), '1', '--solver=wecsim']) == 0
    assert calls == [1.]
    np.testing.assert_array_equal(Result.load(tmp_path/'result.npz').position_m, solved.position_m)
    assert (tmp_path/'result.png').read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    log = (tmp_path/'result.log').read_text()
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
    assert 'VIDEO frames=2 fps=12.5 size=800x800' in report
    assert 'source_end_s=0.080 encoded_duration_s=0.160 verified=true' in report
    assert 'TIMING render preparation_s=' in report
    assert 'verification_s=' in report
    assert 'timings={' not in report
