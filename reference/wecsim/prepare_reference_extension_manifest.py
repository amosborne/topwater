"""Incident-extension readiness archive only; does not evaluate the field or launch a job."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


base = Path('artifacts/reference/reference-extension-20261024')
out = base/'incident-attempt1'
archive = base/'incident-attempt1-source'
assert not out.exists()
archive.mkdir(parents=True,exist_ok=False)
paths = [Path('reference/wecsim')/p for p in ['extend_accepted_reference_packet.py',
    'reference-extension-controls.json','prepare_reference_extension_manifest.py',
    'incident_packet.py','packet-controls.json']]
paths += [Path('reference/capytaine/canonical-boat.json'),Path('reference/capytaine/geometry_execution_guard.py')]
original = Path('artifacts/reference/minimum-20261024')
paths += [original/'packet-attempt3'/p for p in ['packet.npz','diagnostic.json',
    'startup-forces.npz','wave_1.txt','wave_2.txt','wave_3.txt']]
paths += [original/'packet-full-planform/diagnostic.json',original/'consumer/bemio.mat',
          original/'expected-retained-force.mat']
frozen = {}
for path in paths:
    frozen[path.as_posix()] = sha(path)
    dest = archive/path
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(path,dest)
deps = {}
for name in ['numpy','scipy','xarray','threadpoolctl']:
    package = distribution(name)
    deps[name] = dict(version=version(name),file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py','.pyd','.dll','.so'))})
proposal = dict(owner='/root/wecsim_reference',status='ready_for_coordinator_execution_allocation',
    assignment='ONE exactacceptedREF0013packet physicalfield extension/prefix preflight ONLY',
    output=out.as_posix(),input_sha256=frozen,dependency_identity=deps,
    producer_argv=['reference/wecsim/extend_accepted_reference_packet.py','--output',out.as_posix()],
    verification_receipt_name='verification.json',cpu_threads=1,matlab_engines=0,
    maximum_elapsed_seconds=120,expected_elapsed_seconds=[2,20],
    acceptance='Exactoriginalfieldcoefficients/waveprefix→inheritfullplanformcalmbounds; first501×6retainedforceidentity/F0; full60futurepadding',
    physics='NOamplitude/frequency/phase/heading/focus/bandwidth change; ORIGINALH5/CIC10 proposed for laterbody',
    stopping_point='IncidentONLY; exactterminal/idleunderlease/separateverify/finalrelease; no automaticbody/retry')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(archive/'proposal.json')
