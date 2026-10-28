"""Profile the unchanged simple-boat case in private storage outside the checkout.

uv run --group wecsim python scripts/profile_workflow.py --stage render
Use --stage all for boat/sea generation, the full 130 s solve, and all video
frames. Profiles and outputs are retained in the printed OS-temporary directory.
"""
import argparse
import cProfile
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import pstats
import shutil
import subprocess
import tempfile
import time

from topwater.cli import main


def profile_matlab(output):
    """Instrument only the benchmark's native stage; never enable this in the CLI."""
    import matlab.engine
    from topwater import wecsim_backend
    original = wecsim_backend._run_matlab
    def run(workspace, fork):
        engine = matlab.engine.connect_matlab(wecsim_backend.ENGINE_NAME)
        state = engine.eval("profile('status')", nargout=1)
        if state['ProfilerStatus'] != 'off':
            raise RuntimeError('MATLAB profiler is already active; refusing to replace its data')
        engine.eval('profile clear; profile on;', nargout=0)
        # MATLAB Engine permits only one attachment per Python process/session.
        # Leave profiling active in MATLAB, release this client, then let the
        # unmodified production adapter own its ordinary connection and cleanup.
        del engine
        try:
            return original(workspace, fork)
        finally:
            engine = matlab.engine.connect_matlab(wecsim_backend.ENGINE_NAME)
            engine.eval('profile off; topwaterPerformanceInfo=profile("info");', nargout=0)
            try:
                data = engine.eval('jsonencode(topwaterPerformanceInfo.FunctionTable)', nargout=1)
                (output/'matlab_profile.json').write_text(data+'\n', encoding='utf-8')
            finally:
                engine.eval('profile clear; clear topwaterPerformanceInfo', nargout=0)
                del engine
    wecsim_backend._run_matlab = run


def profile(stage, arguments, directory):
    profiler = cProfile.Profile()
    started = time.perf_counter()
    try:
        status = profiler.runcall(main, arguments)
    finally:
        elapsed = time.perf_counter()-started
        profiler.dump_stats(str(directory/f'{stage}.prof'))
    stats = pstats.Stats(profiler)
    rows = [dict(file=file, line=line, function=function, primitive_calls=cc,
                 calls=nc, self_seconds=tt, cumulative_seconds=ct)
            for (file, line, function), (cc, nc, tt, ct, _) in stats.stats.items()]
    print(f'PROFILE {stage}: status={status} wall_seconds={elapsed:.3f}', flush=True)
    return dict(status=status, wall_seconds=elapsed,
                top_cumulative=sorted(rows, key=lambda r: r['cumulative_seconds'], reverse=True)[:40],
                top_self=sorted(rows, key=lambda r: r['self_seconds'], reverse=True)[:40])


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=('solve', 'render', 'all'), default='all')
    parser.add_argument('--output', type=Path, help='Profile directory outside the checkout; default: unique OS temp')
    parser.add_argument('--matlab-profile', action='store_true', help='Also collect the shared MATLAB profiler during solve')
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[1]
    output = (args.output if args.output else Path(tempfile.mkdtemp(prefix='topwater-profile-'))).resolve()
    if output.is_relative_to(project):
        parser.error('--output must be outside the checkout')
    output.mkdir(parents=True, exist_ok=True)
    case = output/'case'
    case.mkdir(exist_ok=True)
    names = ['boat.npz', 'sea.npz']+(['result.npz'] if args.stage == 'render' else [])
    inputs = {}
    for name in names:
        source = project/'examples/simple-boat'/name
        destination = case/name
        if destination.exists():
            parser.error(f'profile input already exists: {destination}; use a fresh output directory')
        # Saves publish by replacement, so profiling never writes through these links.
        try:
            os.link(source, destination)
        except OSError:
            shutil.copyfile(source, destination)
        with source.open('rb') as stream:
            inputs[name] = hashlib.file_digest(stream, 'sha256').hexdigest()
    revision = subprocess.run(['git', '-C', str(project), 'rev-parse', 'HEAD'],
                              capture_output=True, text=True, check=True).stdout.strip()
    sources = sorted((project/'src/topwater').rglob('*.py'))+sorted((project/'src/topwater').rglob('*.m'))
    report = dict(revision=revision, python=platform.python_version(), platform=platform.platform(),
                  cpu_count=os.cpu_count(), input_sha256=inputs,
                  source_sha256={p.relative_to(project).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sources},
                  dependencies={name: version(name) for name in ('numpy', 'scipy', 'numba', 'matplotlib')},
                  matlab_profile=args.matlab_profile, stages={})
    if args.matlab_profile:
        profile_matlab(output)
    commands = {}
    if args.stage == 'all':
        commands['boat'] = ['boat', str(case), '--length', '7.3152', '--beam', '2.7432',
                            '--deadrise', '20', '--mass', '2500']
        commands['sea'] = ['generate-sea', str(case), '--record', 'ndbc-46086-20231223T0040Z', '--seed', '42']
    if args.stage in ('solve', 'all'):
        commands['solve'] = ['solve', str(case), '100', '--rest', '10', '--ramp', '5']
    if args.stage in ('render', 'all'):
        commands['render'] = ['render', str(case)]
    print(f'PROFILE directory={output}', flush=True)
    for stage, arguments in commands.items():
        report['stages'][stage] = profile(stage, arguments, output)
        (output/'profile.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        if report['stages'][stage]['status']:
            return report['stages'][stage]['status']
    return 0


if __name__ == '__main__':
    raise SystemExit(run())
