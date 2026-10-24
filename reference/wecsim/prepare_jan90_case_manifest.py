"""Archive prospective native case identities; no engine connection or numerical launch."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


root = Path.cwd()
base = Path('artifacts/reference/generated-bank-consumer-20261024')
archive = base/'case-jan90-attempt1-source'
archive.mkdir(parents=True, exist_ok=False)
sources = ['run_jan90_case.m', 'prepare_jan90_case_inputs.m', 'wecSimInputFile_jan90.m',
           'build_generated_model.m', 'apply_jan90_force_table.m', 'run_matlab.py',
           'prepare_jan90_case_manifest.py']
paths = [Path('reference/wecsim')/name for name in sources]
paths += [Path('reference/capytaine/geometry_execution_guard.py'), Path('reference/capytaine/canonical-boat.json'),
          Path('docs/experts/wecsim-fork.json')]
fork = Path('artifacts/dependencies/topwater-wec-sim/source')
paths += [p for p in fork.rglob('*') if p.is_file()]
for folder in ['generated-bank-native-20261024/attempt1', 'radiation-heading-scout-20261024/attempt1',
               'generated-bank-consumer-20261024/bemio-attempt1', 'generated-bank-consumer-20261024/jan90-attempt1']:
    paths += [p for p in (Path('artifacts/reference')/folder).iterdir() if p.is_file()]
bank = Path('artifacts/reference/generated-sea-candidates-20261024/ndbc-46086-20250105T1140Z')
paths += [bank/'realization-global.json', bank/'heading-offset-90deg.json']
frozen = {}
for path in dict.fromkeys(paths):
    frozen[path.as_posix()] = sha(path)
    target = archive/path
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(path, target)
dependencies = {}
for name in ['matlabengine', 'numpy', 'scipy', 'xarray', 'threadpoolctl']:
    package = distribution(name)
    dependencies[name] = dict(version=version(name), file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py', '.pyd', '.dll', '.so'))})
proposal = dict(owner='/root/wecsim_reference', status='ready_for_coordinator_execution_allocation',
    assignment='ONE provisional17panel Jan05offset90 native held-release60s; no automatic rerun/refinement',
    output=str(base/'case-jan90-attempt1'),
    producer_argv=['reference/wecsim/run_matlab.py', 'reference/wecsim/run_jan90_case.m',
                   '--log', str(base/'case-jan90-attempt1.log')],
    input_sha256=frozen, dependency_identity=dependencies, cpu_threads=1, matlab_engines=1,
    maximum_elapsed_seconds=1800, expected_elapsed_seconds=[30, 180], verification_receipt_name='verification.json',
    fork_commit='f962103a348d63101ce669c2611f7d7c109755c3',
    initial_state='held before0; q0=[0,0,eta0,atan(sy/sqrt(1+sx^2)),-atan(sx),0]; v0=0 and radiationhistory0',
    radiation='existing native40sCIC dt.02; no coefficient fitting or scout substitution',
    forcing='full384 unscaled native pairH, generatedlocalphase, exp(-iwt); stationary sea present at0; no ramp',
    expected='3001raw60columnnative samples/fullsigned6forcebus vs exacttable/nonzeroq0/v0/memory/frame/mass readbacks',
    stopping_point='Exact process terminal and runtime evidence; no automatic renderer or refined case',
    physical_limits='Tangentplane pose not curvedsurface equilibrium; heldpast uses meanbody diffraction assumption; mesh/radiation/geometry accuracy unqualified')
(archive/'proposal.json').write_text(json.dumps(proposal, indent=2)+'\n', encoding='utf-8')
print(archive/'proposal.json')
