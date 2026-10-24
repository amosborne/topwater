"""Prepare byte-identical official assets and prospective manifests; never launch MATLAB."""
import hashlib
import json
import shutil
import subprocess
from importlib.metadata import distribution, version
from pathlib import Path

ROOT = Path.cwd()
FORK = ROOT/'artifacts/dependencies/topwater-wec-sim'
BASE = 'a23d2993807410fe78cc63546c1ba9d81912024f'
PIN = 'f962103a348d63101ce669c2611f7d7c109755c3'
OUT = ROOT/'artifacts/reference/official-examples-20261024-ready2'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')


assert subprocess.check_output(['git','-C',str(FORK),'rev-parse','HEAD'],text=True).strip() == PIN
assert not subprocess.check_output(['git','-C',str(FORK),'status','--porcelain'],text=True).strip()
OUT.mkdir(parents=True, exist_ok=False)
for key, folder, duration in [('standard','examples/RM3',400),
                              ('cic','tests/RegressionTests/RegularWaves/regularCIC',200)]:
    case_root = OUT/key
    tree = case_root/'official'
    shutil.copytree(FORK/'examples/RM3', tree/'examples/RM3')
    if key == 'cic':
        shutil.copytree(FORK/folder, tree/folder)
    copied = {}
    official_proof = {}
    for path in tree.rglob('*'):
        if not path.is_file():
            continue
        rel = path.relative_to(tree).as_posix()
        blob = subprocess.check_output(['git','-C',str(FORK),'show',f'{BASE}:{rel}'])
        # Checkout CRLF conversion is allowed; compare the repository blob identity.
        tracked = subprocess.check_output(['git','-C',str(FORK),'hash-object','--path',rel,str(FORK/rel)],text=True).strip()
        original = subprocess.check_output(['git','-C',str(FORK),'rev-parse',f'{BASE}:{rel}'],text=True).strip()
        assert tracked == original, rel
        assert sha(path) == sha(FORK/rel)
        copied[path.relative_to(ROOT).as_posix()] = sha(path)
        official_proof[rel] = {'base_git_blob':original,'base_blob_sha256':hashlib.sha256(blob).hexdigest(),
                               'checkout_sha256':sha(path),'bytes':path.stat().st_size}
    config = dict(case_directory=(tree/folder).relative_to(ROOT).as_posix(),
                  golden_available=key == 'cic', duration_s=duration)
    write(case_root/'case.json',config)
    entry = case_root/'entry.m'
    entry.write_text("addpath('C:/Users/alex/topwater/reference/wecsim');\n"
                     f"run_official_example('{key}');\n",encoding='utf-8')
    frozen_paths = [ROOT/'reference/wecsim/run_official_example.m',ROOT/'reference/wecsim/run_matlab.py',
                    ROOT/'reference/wecsim/prepare_official_examples.py',ROOT/'reference/capytaine/geometry_execution_guard.py',
                    case_root/'case.json',entry]
    frozen_paths += [p for p in (FORK/'source').rglob('*') if p.is_file()]
    for path in frozen_paths:
        copied[path.relative_to(ROOT).as_posix()] = sha(path)
    archive = case_root/'source'
    for rel in copied:
        dest = archive/rel
        dest.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(ROOT/rel,dest)
    package = distribution('matlabengine')
    deps = {'matlabengine':dict(version=version('matlabengine'),file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py','.pyd','.dll','.so'))})}
    proposal = dict(owner='/root/wecsim_reference',status='ready_for_coordinator_execution_allocation',
        assignment=f'ONE official {key} RM3 baseline; original input/model/H5 physics unchanged',
        output=(tree/folder).relative_to(ROOT).as_posix(),input_sha256=copied,dependency_identity=deps,
        producer_argv=['reference/wecsim/run_matlab.py',entry.relative_to(ROOT).as_posix(),
                       '--log',(case_root/'engine.log').relative_to(ROOT).as_posix()],
        verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
        expected_elapsed_seconds=[30,240],cpu_threads=1,matlab_engines=1,
        official_base=BASE,maintained_pin=PIN,official_asset_proof=official_proof,
        display_only_override='native MCR simu.explorer=off; defaultFigureVisible=off',
        stopping_point='exact terminal, frozen identity, idle under lease/separate readback/final release; no retries',
        acceptance='finite complete official grid; zero initial velocity/radiation; official CIC golden1e-5m heave',
        limits='Open system with incident forcing/PTO/constraints. Not custom boat BEM validation or passivity proof.')
    write(case_root/'proposal.json',proposal)
    print(case_root/'proposal.json')
