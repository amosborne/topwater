"""Archive the one prospective initialization fixture; never launch an engine/job."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


base = Path('artifacts/reference/generated-bank-consumer-20261024')
archive = base/'seed-initialization-fixture-attempt3-source'
archive.mkdir(parents=True, exist_ok=False)
names = ['check_seeded_initialization.m', 'cic_seed_table.m', 'cic_seed_lookup.m',
         'harmonic_seed_state.m', 'make_harmonic_seed.m', 'set_initial_joint_state.m',
         'apply_cic_history_adapter.m', 'run_matlab.py', 'prepare_seed_fixture_manifest.py']
paths = [Path('reference/wecsim')/name for name in names]
paths += [Path('reference/capytaine/geometry_execution_guard.py'), Path('docs/experts/wecsim-fork.json'),
          Path('docs/experts/wecsim_expert/seed-initialization-fixture-proposal-2026-10-24.md'),
          Path('artifacts/reference/generated-bank-native-20261024/attempt1/boat.stl')]
paths += [p for p in Path('artifacts/dependencies/topwater-wec-sim/source').rglob('*') if p.is_file()]
frozen = {}
for path in paths:
    frozen[path.as_posix()] = sha(path)
    target = archive/path
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
package = distribution('matlabengine')
dependencies = {'matlabengine':dict(version=version('matlabengine'), file_sha256={str(p):sha(package.locate_file(p))
    for p in package.files if str(p).endswith(('.py', '.pyd', '.dll', '.so'))})}
proposal = dict(owner='/root/wecsim_reference', status='ready_for_coordinator_execution_allocation',
    assignment='ONE nativeCIC/prehistory/PCHIP and joint-initialization fixture; no canonicalbodytrajectory',
    output=str(base/'seed-initialization-fixture-attempt3'),
    producer_argv=['reference/wecsim/run_matlab.py', 'reference/wecsim/check_seeded_initialization.m',
                   '--log', str(base/'seed-initialization-fixture-attempt3.log')],
    input_sha256=frozen, dependency_identity=dependencies, cpu_threads=1, matlab_engines=1,
    maximum_elapsed_seconds=600, expected_elapsed_seconds=[10, 60], verification_receipt_name='verification.json',
    expected='Directindexedfullpast6x6force rows/nativePCHIP queries agree; correctionvanishes; actualnativejointCGpose/v/worldomega agrees',
    exclusions=['No 60s body trajectory, BEM, coefficient/input/fork change or custom ODE',
                'No acceleration-delay interface/state change', 'No automatic retry or next case'])
(archive/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n')
print(archive/'proposal.json')
