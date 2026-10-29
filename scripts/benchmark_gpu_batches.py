"""Sweep complete GPU solves in fresh processes; retain logs, hashes and timings.

Archives are hashed member-by-member and removed after verification to bound
temporary disk use. The default selection is the smallest batch within 1% of
the fastest median, using three shuffled rounds after one warm-up solve.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import statistics
import subprocess
import tempfile
import zipfile


def archive_hashes(path):
    with zipfile.ZipFile(path) as archive:
        result = {}
        for name in sorted(archive.namelist()):
            with archive.open(name) as stream:
                result[name] = hashlib.file_digest(stream, 'sha256').hexdigest()
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batches', nargs='+', type=int, default=[32, 64, 128, 256, 512, 1024])
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1 or min(args.batches) < 1:
        parser.error('batch sizes and repeats must be positive')
    project = Path(__file__).resolve().parents[1]
    output = Path(tempfile.mkdtemp(prefix='topwater-batch-sweep-')).resolve()
    report = dict(scope='Complete unprofiled simple-boat GPU solve, including publication, excluding rendering',
                  output_directory=str(output), repeats=args.repeats, batches=sorted(set(args.batches)),
                  selection='Smallest batch within 1% of fastest median', seed=42, runs=[])
    report['benchmark_source_sha256'] = {str(path.relative_to(project)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (Path(__file__).resolve(), project/'scripts/profile_workflow.py')}
    reference = None
    print(f'SWEEP directory={output}', flush=True)
    randomizer = random.Random(42)
    schedule = [(0, 128)]
    for repeat in range(1, args.repeats+1):
        batches = report['batches'].copy()
        randomizer.shuffle(batches)
        schedule.extend((repeat, batch) for batch in batches)
    for index, (repeat, batch) in enumerate(schedule):
        directory = output/f'{index:02d}-round{repeat}-batch{batch}'
        log = output/f'{index:02d}-round{repeat}-batch{batch}.log'
        command = ['uv', 'run', '--no-sync', 'python', 'scripts/profile_workflow.py',
                   '--stage', 'solve', '--device', 'gpu', '--no-profile',
                   '--radiation-chunk-points', str(batch), '--output', str(directory)]
        with log.open('w', encoding='utf-8') as stream:
            completed = subprocess.run(command, cwd=project, stdout=stream, stderr=subprocess.STDOUT)
        if completed.returncode:
            raise RuntimeError(f'Batch {batch} failed; inspect {log}')
        profile = json.loads((directory/'profile.json').read_text(encoding='utf-8'))
        if reference is None:
            report['provenance'] = {key: value for key, value in profile.items()
                                    if key not in ('stages', 'radiation_chunk_points', 'solve_backend_events')}
        else:
            for key in ('source_sha256', 'dependencies', 'input_sha256'):
                assert profile[key] == report['provenance'][key], f'Changed benchmark inputs: {key}'
            assert any('SOURCE_CACHE hit' in event for event in profile['solve_backend_events'])
        archive = (directory/'case/result.npz').resolve()
        assert archive.is_relative_to(output)
        hashes = archive_hashes(archive)
        if reference is None:
            reference = hashes
            report['reference_member_sha256'] = hashes
        assert hashes == reference, f'Batch {batch} changed output; inspect {archive}'
        archive.unlink()
        run = dict(repeat=repeat, batch=batch, warmup=repeat == 0,
                   complete_solve_seconds=profile['stages']['solve']['wall_seconds'],
                   output_identical=True, directory=str(directory),
                   backend_events=profile['solve_backend_events'])
        report['runs'].append(run)
        (output/'sweep.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print(f'SWEEP round={repeat} batch={batch} complete_solve_s={run["complete_solve_seconds"]:.3f} '
              f'outputs=identical warmup={repeat == 0}', flush=True)
    report['medians_seconds'] = {str(batch): statistics.median(
        run['complete_solve_seconds'] for run in report['runs'] if not run['warmup'] and run['batch'] == batch)
        for batch in report['batches']}
    best = min(report['medians_seconds'].values())
    report['selected_batch'] = min(int(batch) for batch, seconds in report['medians_seconds'].items()
                                   if seconds <= best*1.01)
    (output/'sweep.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print('SWEEP summary='+json.dumps(dict(medians=report['medians_seconds'],
                                          selected_batch=report['selected_batch'])), flush=True)


if __name__ == '__main__':
    main()
