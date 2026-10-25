"""Freeze fixed physicaltrain preflight only; no field computation or numerical launch."""
import argparse
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument('--base',type=Path,default=Path('artifacts/reference/calm-buoy-20261024'))
parser.add_argument('--train-name',default='train-attempt2')
parser.add_argument('--controls',type=Path)
parser.add_argument('--entry',default='reference/wecsim/run_calm_buoy_case.m')
args = parser.parse_args()
base = args.base
out = base/args.train_name
archive = base/(args.train_name+'-source')
assert not out.exists()
archive.mkdir(parents=True,exist_ok=False)
paths = [Path('reference/wecsim')/p for p in ['calm_buoy_train.py','prepare_calm_buoy_train.py',
    'calm-buoy-train-controls.json','prepare_calm_buoy_train_manifest.py']]
paths += [Path('reference/capytaine/canonical-boat.json'),Path('reference/capytaine/geometry_execution_guard.py')]
paths += [base/'export-attempt1'/p for p in ['components.npz','case.json','observation.json']]
paths += [Path('artifacts/reference/calm-buoy-20261024')/'excitation-attempt1/excitation.npz',Path('artifacts/reference/calm-buoy-20261024')/'excitation-attempt1/manifest.json']
if args.controls:
    paths.remove(Path('reference/wecsim/calm-buoy-train-controls.json'))
    paths.append(args.controls)
    paths += [Path('artifacts/reference/calm-buoy-20261024/export-attempt1')/p for p in ['case.json','components.npz']]
frozen = {}
for path in paths:
    frozen[path.as_posix()] = sha(path)
    dest = archive/path
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(path,dest)
deps = {}
for name in ['numpy','scipy','threadpoolctl']:
    package = distribution(name)
    deps[name] = dict(version=version(name),file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py','.pyd','.dll','.so'))})
proposal = dict(owner='/root/wecsim_reference',status='ready_for_coordinator_execution_allocation',
    assignment='ONE data-derived buoy flat-top spatialtrain preflight; no temporalramp/bodyhistoryseed',
    output=out.as_posix(),input_sha256=frozen,dependency_identity=deps,
    producer_argv=['reference/wecsim/prepare_calm_buoy_train.py','--output',out.as_posix()]+(['--controls',args.controls.as_posix()] if args.controls else []),
    verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
    expected_elapsed_seconds=[15,180],cpu_threads=1,matlab_engines=0,
    physical_construction='W(s/cg+declared tc), plateau15/support25/quintictransition; exact signedκ dispersion; generatorcarriersunscaled/unfitted',
    support='sharedOmega0:.005:lambda>=2beamcap endpoint; signedκ±; analyticpressureH0/quadraticOmega²bridge below.1; upperremovedsupport/reversesides reported/no renorm',
    stage_acceptance='With separate_assessment true exit0 means artifacts complete ONLY, never calm acceptance; subsequent engineering assessment required',
    order='coefficients/support/actual60PSD/cheapall6F0 FIRST; efficient percomponentinterval rectangle bounds ONLYifF0pass',
    accuracy='NewfullheadingH linearinterpolation and fixed17mesh unqualified, originalradiationH5/CIC10 unchanged',
    stopping_point='Exactterminal/allfrozen/idleUNDERlease/separatereadback/FINALrelease; no automaticbody/retiming/search/encoder')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(archive/'proposal.json')
