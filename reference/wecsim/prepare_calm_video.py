"""Bookkeeping only: bind accepted inputs and archive ONE prospective render job.

No field evaluation, renderer import, model, engine, lease or process launch.
Run only after the coordinator has written the genuine accepted-motion receipt.
"""
import argparse
import hashlib
import json
from importlib.metadata import distribution, version
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'artifacts/reference/calm-buoy-20261024'


def read(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def key(path):
    resolved = path.resolve()
    try:
        return resolved.relative_to(ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--stage', choices=['preview', 'encode'], required=True)
parser.add_argument('--receipt', type=Path)
parser.add_argument('--preview-review', type=Path)
parser.add_argument('--owner', default='/root/calm_video')
args = parser.parse_args()
controls_path = BASE/'render-ready-controls.json'
output = BASE/('render-'+args.stage+'-attempt1')
archive = BASE/('render-'+args.stage+'-attempt1-source')
assert not output.exists() and not archive.exists(), 'Preserve prior outputs/source archives'

if args.stage == 'preview':
    assert args.receipt and not controls_path.exists()
    receipt_path = args.receipt.resolve()
    receipt = read(receipt_path)
    assert receipt['decision'] == 'accepted_for_visualization' and receipt['reviewed_by']
    assert receipt['observation_id'] == 'ndbc-46086-20230319T0840Z'
    assert receipt['trajectory_status'] in ('provisional_model_trajectory', 'accepted_reference_case')
    controls = read(ROOT/'reference/wecsim/calm-video-controls.json')
    assert controls['status'] == 'source_only_pending_accepted_case'
    controls['motion_acceptance_receipt'] = key(receipt_path)
    required = ['physical_case', 'heading_file', 'packet_npz', 'incident_report',
                'motion_mat', 'summary_json', 'native_output_mat', 'motion_acceptance_receipt']
    input_paths = [ROOT/controls[name] for name in required]
    input_paths += [ROOT/'reference/wecsim'/name for name in
                   ('render_calm_buoy_reference.py', 'render_generated_reference.py', 'calm_buoy_train.py')]
    controls['input_sha256'] = {key(path):sha(path) for path in input_paths}
    for name, receipt_name in [('motion_mat', 'motion_sha256'), ('summary_json', 'summary_sha256'),
                              ('native_output_mat', 'native_output_sha256'), ('heading_file', 'heading_sha256'),
                              ('packet_npz', 'packet_sha256'), ('incident_report', 'incident_report_sha256')]:
        assert controls['input_sha256'][controls[name]] == receipt[receipt_name], name+' receipt mismatch'
    report = read(ROOT/controls['incident_report'])
    summary = read(ROOT/controls['summary_json'])
    assert report['complete'] and report['accepted']
    assert summary['complete'] and summary['samples'] == 3001 and summary['duration_s'] == 60
    controls['status'] = 'ready_for_coordinator_render_allocation'
    write(controls_path, controls)
else:
    assert args.preview_review and controls_path.is_file()
    controls = read(controls_path)
    assert controls['status'] == 'ready_for_coordinator_render_allocation'

frozen = dict(controls['input_sha256'])
assert all(sha(ROOT/path) == digest for path, digest in frozen.items())
frozen[key(controls_path)] = sha(controls_path)
frozen[key(ROOT/'reference/capytaine/geometry_execution_guard.py')] = sha(ROOT/'reference/capytaine/geometry_execution_guard.py')
frozen[key(Path(__file__))] = sha(Path(__file__))
command = ['reference/wecsim/render_calm_buoy_reference.py', '--config', key(controls_path),
           '--output', key(output), '--mode', args.stage]
if args.stage == 'encode':
    review_path = args.preview_review.resolve()
    review = read(review_path)
    assert review['decision'] == 'preview_pass' and review['reviewed_by']
    preflight_path = ROOT/review['preflight_json']
    cache_path = ROOT/review['surface_cache_npz']
    assert sha(preflight_path) == review['preflight_sha256']
    assert sha(cache_path) == review['surface_cache_sha256']
    preflight = read(preflight_path)
    assert preflight['complete'] and preflight['mode'] == 'preview'
    expected_identity = dict(source_sha256=sha(ROOT/'reference/wecsim/render_calm_buoy_reference.py'),
                             config_sha256=sha(controls_path), input_sha256=controls['input_sha256'])
    assert preflight['identity'] == expected_identity
    assert preflight['surface_cache_sha256'] == review['surface_cache_sha256']
    for path in (review_path, preflight_path, cache_path):
        frozen[key(path)] = sha(path)
    command += ['--preview-review', key(review_path)]

ffmpeg = shutil.which('ffmpeg')
assert ffmpeg, 'FFmpeg executable unavailable'
frozen[key(Path(ffmpeg))] = sha(Path(ffmpeg))
deps = {}
for name in ('numpy', 'scipy', 'matplotlib', 'threadpoolctl', 'pillow', 'contourpy', 'kiwisolver',
             'cycler', 'fonttools', 'packaging', 'pyparsing', 'python-dateutil', 'six'):
    package = distribution(name)
    deps[name] = dict(version=version(name), file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py', '.pyd', '.dll', '.so', '.ttf', '.afm', '.mplstyle'))})

archive.mkdir(parents=True)
for name, digest in frozen.items():
    path = ROOT/name
    # External binary identity is frozen, but no installed dependency is modified/copied.
    try:
        relative = path.resolve().relative_to(ROOT)
    except ValueError:
        continue
    destination = archive/'frozen'/relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, destination)
    assert sha(destination) == digest
proposal = dict(owner=args.owner, status='ready_for_coordinator_execution_allocation',
    assignment='ONE acceptedmildbuoy60s '+args.stage+' rendering; no physics/field changes',
    output=key(output), producer_argv=command, input_sha256=frozen, dependency_identity=deps,
    verification_receipt_name='verification.json', cpu_threads=1, matlab_engines=0,
    maximum_elapsed_seconds=controls['preview_maximum_seconds' if args.stage == 'preview' else 'encode_maximum_seconds'],
    expected_elapsed_seconds=[30, 300] if args.stage == 'preview' else [120, 600],
    stopping_point='Preview requires independent literal/provenance review; encoder separately allocated; preserve partial outputs',
    numerical_input='Exact accepted immutable train2/native60s motion; visual cache never used for forcing')
write(archive/'proposal.json', proposal)
print(key(archive/'proposal.json'))
