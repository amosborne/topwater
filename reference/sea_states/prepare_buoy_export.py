"""Freeze an unchanged 2023 observed-buoy exporter at fresh caller-selected paths."""
import argparse
from importlib.metadata import distribution, version
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--config', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
config_path = (ROOT/args.config).resolve()
output = (ROOT/args.output).resolve()
assert config_path.is_relative_to(ROOT)
assert output.is_relative_to(ROOT/'artifacts/reference')
archive = output.with_name(output.name+'-source')
assert not output.exists() and not archive.exists()
config = json.loads(config_path.read_text(encoding='utf-8-sig'))
# The unchanged exporter deliberately attaches the existing five2023 raw files.
assert config['time_utc'].startswith('2023-')
exporter = ROOT/'reference/sea_states/export_calm_buoy.py'
assert sha(exporter) == '22bcf2c0e95f6f021eb7ffb9c1d3b5b34ffd0e1b954bbcd4e1d21ddde38e9621'
paths = [config_path, Path(__file__).resolve(), exporter,
    ROOT/'reference/capytaine/geometry_execution_guard.py', ROOT/config['physical_input'],
    ROOT/config['catalog_provenance'], Path(config['source_library'])]
paths += [ROOT/'src/topwater'/name for name in ('__init__.py', '_validation.py', 'inputs.py', 'sea_states.py')]
archive.mkdir(parents=True)
frozen = {}
for path in paths:
    assert path.is_file()
    key = path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()
    frozen[key] = sha(path)
    target = archive/(path.relative_to(ROOT) if path.is_relative_to(ROOT) else Path('external')/path.name)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
dependencies = {}
for name in ('numpy', 'scipy'):
    package = distribution(name)
    dependencies[name] = dict(version=version(name), file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py', '.pyd', '.dll', '.so'))})
proposal = dict(owner='/root/calm_buoy_input', status='ready_for_coordinator_execution_allocation',
    assignment='ONE actual stronger observation export through unchanged existing generator, no field/BEM/dynamics',
    output=output.relative_to(ROOT).as_posix(), input_sha256=frozen, dependency_identity=dependencies,
    producer_argv=[exporter.relative_to(ROOT).as_posix(), '--config',
        config_path.relative_to(ROOT).as_posix(), '--output', output.relative_to(ROOT).as_posix()],
    verification_receipt_name='verification.json', maximum_elapsed_seconds=60,
    expected_elapsed_seconds=[1, 10], cpu_threads=1, matlab_engines=0,
    observation=config['observation_id'], source_library_sha256=config['source_library_sha256'],
    seed=config['seed'], components_per_band=config['components_per_band'],
    amplitudes_and_phases='Existing generator, unchanged/no scaling or fitting',
    stopping_point='Neutral stationary carriers only; separate physical incoming train/calm assessment and native body/video allocation')
(archive/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n', encoding='utf-8')
print(archive/'proposal.json')
