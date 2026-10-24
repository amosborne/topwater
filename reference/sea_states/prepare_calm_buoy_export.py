"""Freeze one generator-export proposal; never generate waves or run a solver."""
from importlib.metadata import distribution, version
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


config_path = ROOT/'reference/sea_states/calm-buoy-export.json'
config = json.loads(config_path.read_text())
archive = ROOT/'artifacts/reference/calm-buoy-20261024/export-attempt1-source'
output = ROOT/'artifacts/reference/calm-buoy-20261024/export-attempt1'
assert not archive.exists() and not output.exists()
paths = [config_path, Path(__file__), ROOT/'reference/sea_states/export_calm_buoy.py',
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
    assignment='ONE actual mild observation export through existing generator, no field/BEM/dynamics',
    output=output.relative_to(ROOT).as_posix(), input_sha256=frozen, dependency_identity=dependencies,
    producer_argv=['reference/sea_states/export_calm_buoy.py', '--config',
        config_path.relative_to(ROOT).as_posix(), '--output', output.relative_to(ROOT).as_posix()],
    verification_receipt_name='verification.json', maximum_elapsed_seconds=60,
    expected_elapsed_seconds=[1, 10], cpu_threads=1, matlab_engines=0,
    observation=config['observation_id'], source_library_sha256=config['source_library_sha256'],
    seed=42, components_per_band=8, amplitudes_and_phases='Existing generator, unchanged/no scaling or fitting',
    stopping_point='Neutral stationary carriers only; separate physical incoming finite-train preflight and native body/video allocation')
(archive/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n', encoding='utf-8')
print(archive/'proposal.json')
