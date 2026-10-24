"""Archive exact identities for an already allocated consumer job; never launch it."""
import argparse
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


p = argparse.ArgumentParser()
p.add_argument('stage', choices=['bemio', 'harmonic', 'compare', 'jan90'])
args = p.parse_args()
base = Path('artifacts/reference/generated-bank-consumer-20261024')
archive = base/(args.stage+'-attempt1-source')
archive.mkdir(parents=True, exist_ok=False)
native = Path('artifacts/reference/generated-bank-native-20261024/attempt1')
paths = [path for path in native.iterdir() if path.is_file()]
paths += [Path('reference/wecsim/prepare_consumer_manifest.py'),
          Path('reference/capytaine/geometry_execution_guard.py')]
if args.stage == 'bemio':
    paths += [Path('reference/wecsim/prepare_generated_bemio.m'), Path('reference/wecsim/run_matlab.py')]
    paths += [path for path in Path('artifacts/dependencies/topwater-wec-sim/source').rglob('*') if path.is_file()]
    argv = ['reference/wecsim/run_matlab.py', 'reference/wecsim/prepare_generated_bemio.m',
            '--log', str(base/'bemio-attempt1.log')]
elif args.stage == 'harmonic':
    paths += [Path('reference/wecsim/screen_generated_banks.py'), Path('reference/wecsim/check_generated_h5.m'),
              Path('reference/wecsim/run_matlab.py'), Path('reference/wecsim/run_generated_harmonic_job.py')]
    paths += [path for path in (base/'bemio-attempt1').iterdir() if path.is_file()]
    argv = ['reference/wecsim/run_generated_harmonic_job.py']
elif args.stage == 'compare':
    paths += [Path('reference/wecsim/screen_generated_banks.py'), Path('reference/wecsim/compare_generated_radiation.py')]
    paths += [path for path in (base/'bemio-attempt1').iterdir() if path.is_file()]
    incident = Path('artifacts/reference/incident-input-diagnostics-20261024')
    paths += list(incident.glob('*-offset-0deg.npz'))
    paths += list(incident.glob('*-offset-0deg.json'))+[incident/'manifest.json']
    argv = ['reference/wecsim/compare_generated_radiation.py', '--output', str(base/'compare-attempt1')]
else:
    paths += [Path('reference/wecsim/compare_jan90.py'), Path('reference/wecsim/screen_generated_banks.py')]
    paths += [path for path in (base/'bemio-attempt1').iterdir() if path.is_file()]
    scout = Path('artifacts/reference/radiation-heading-scout-20261024/attempt1')
    paths += [path for path in scout.iterdir() if path.is_file()]
    incident = Path('artifacts/reference/incident-input-diagnostics-20261024')
    paths += [incident/'ndbc-46086-20250105T1140Z-offset-90deg.npz', incident/'ndbc-46086-20250105T1140Z-offset-90deg.json']
    bank = Path('artifacts/reference/generated-sea-candidates-20261024/ndbc-46086-20250105T1140Z')
    paths += [bank/'heading-offset-90deg.json', bank/'realization-global.json']
    argv = ['reference/wecsim/compare_jan90.py', '--scout', str(scout), '--output', str(base/'jan90-attempt1')]
frozen = {}
for path in paths:
    key = path.as_posix()
    frozen[key] = sha(path)
    target = archive/key
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
deps = {}
for name in ['numpy', 'scipy', 'xarray', 'threadpoolctl', 'matlabengine']:
    package = distribution(name)
    deps[name] = dict(version=version(name), file_sha256={str(path):sha(package.locate_file(path))
        for path in package.files if str(path).endswith(('.py', '.pyd', '.dll', '.so'))})
proposal = dict(owner='/root/wecsim_reference', status='ready_for_coordinator_execution_allocation' if args.stage in ('compare', 'jan90')
    else 'ready_for_coordinator_execution_allocation_granted',
    assignment='REF-002 allocated '+args.stage+' ONLY; no body trajectory or custom integrator',
    output=str(base/(args.stage+'-attempt1')), producer_argv=argv,
    input_sha256=frozen, dependency_identity=deps, cpu_threads=1,
    maximum_elapsed_seconds=600 if args.stage == 'bemio' else (180 if args.stage == 'harmonic' else 120),
    verification_receipt_name='verification.json',
    expected='OfficialBEMIO40s/phase/exactmetadata pass' if args.stage == 'bemio'
        else ('All3 stationary harmonic banks/force signs; radiationdiagnostics, no qualification' if args.stage == 'harmonic'
              else 'Storeddata native-vs40sK all6particular histories, finitebandKK/normalization and wave-relative proxies; no qualification'))
(archive/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n', encoding='utf-8')
print(archive/'proposal.json')
