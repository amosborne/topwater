"""Prepare ONE prospective official CIC PTO-only control; no simulation or model editing."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


fork = Path('artifacts/dependencies/topwater-wec-sim')
baseline = Path('artifacts/reference/official-examples-20261024-ready2/cic')
old = json.loads((baseline/'proposal.json').read_text())
base = Path('artifacts/reference/official-pto-zero-20261024/cic')
assert not base.exists()
tree = base/'official'
shutil.copytree(fork/'examples/RM3',tree/'examples/RM3')
folder = Path('tests/RegressionTests/RegularWaves/regularCIC')
shutil.copytree(fork/folder,tree/folder)
paths = []
for path in tree.rglob('*'):
    if path.is_file():
        relative = path.relative_to(tree).as_posix()
        assert sha(path) == old['official_asset_proof'][relative]['checkout_sha256']
        paths.append(path)
(base/'case.json').write_text(json.dumps(dict(case_directory=(tree/folder).as_posix(),
    golden_available=False,duration_s=200),indent=2)+'\n')
paths += [base/'case.json',baseline/'proposal.json',baseline/'official'/folder/'audit.mat']
paths += [Path('reference/wecsim')/p for p in ['run_official_pto_zero.m',
    'audit_official_pto_energy.m','prepare_official_pto_zero.py','run_matlab.py']]
paths += [Path('reference/capytaine/geometry_execution_guard.py')]
paths += list((fork/'source').rglob('*'))
archive = base/'source'
frozen = {}
for path in paths:
    if not path.is_file():
        continue
    frozen[path.as_posix()] = sha(path)
    dest = archive/path
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(path,dest)
package = distribution('matlabengine')
deps = {'matlabengine':dict(version=version('matlabengine'),file_sha256={str(p):sha(package.locate_file(p))
    for p in package.files if str(p).endswith(('.py','.pyd','.dll','.so'))})}
proposal = dict(owner='/root/wecsim_reference',status='ready_for_coordinator_execution_allocation',
    assignment='ONE official provided regularCIC control: sole physical change PTO damping1200000 to0',
    output=(tree/folder).as_posix(),input_sha256=frozen,dependency_identity=deps,
    producer_argv=['reference/wecsim/run_matlab.py','reference/wecsim/run_official_pto_zero.m',
        '--log',(base/'engine.log').as_posix()],
    verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
    expected_elapsed_seconds=[30,240],cpu_threads=1,matlab_engines=1,
    unchanged='Official2heave bodies/model/H5/M/K/Ainf/ramp100/CIC20/dt.1/200s/allotherphysicalsettings',
    override='OfficialnativeMCR pto(1).damping=0 plus alreadydeclared explorerOFF/displayonly',
    comparison='Against passed officialCIC baseline. Modifiedtrajectory has NO suppliedgolden and NO copied1e-5 expectation.',
    energy='ActualSolid masssplit plus restoring storage; signed Fe/Frad/APPLIEDAr/PTO work; report openbudget residuals/scales, not conservation assertion',
    stopping_point='Exactterminal/idleunderlease/separateverify/finalrelease; no automatic nextstep/rerun/model/freeDOF changes')
(base/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(base/'proposal.json')
