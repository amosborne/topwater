"""Archive read-only train2 reassessment sources/inputs; no field computation or launch."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


base = Path('artifacts/reference/calm-buoy-20261024')
out = base/'assessment-attempt1'
archive = base/'assessment-attempt1-source'
assert not out.exists()
archive.mkdir(parents=True,exist_ok=False)
paths = [Path('reference/wecsim')/p for p in ['assess_calm_buoy_train.py',
    'calm-buoy-assessment-controls.json','prepare_calm_buoy_assessment_manifest.py','calm_buoy_train.py']]
paths += [Path('reference/capytaine/canonical-boat.json'),Path('reference/capytaine/geometry_execution_guard.py')]
paths += [p for p in (base/'train-attempt2').iterdir() if p.is_file()]
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
    assignment='ONE read-only reassessment immutablemildtrain2 with explicitREF002engineeringcalmpolicy +fullhullbounds',
    output=out.as_posix(),input_sha256=frozen,dependency_identity=deps,
    producer_argv=['reference/wecsim/assess_calm_buoy_train.py','--output',out.as_posix()],
    verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
    expected_elapsed_seconds=[30,180],cpu_threads=1,matlab_engines=0,
    physics='NO train/parent/C/force/H regeneration/window/retime/historyseed; old.1Nrejects remain immutable',
    policy='eps.001 normF<=epsmg/normM<=epsmgLchar ANDperaxis<=epsS samefixednativebasis; unchangedeta/u/arrival',
    field_bound='Combined p6Taylor beforeabs, rigorous7th-orderexp remainder;1mXY/.2mz boxes, batch8, unresolvedsubdivision max3levels within480s',
    interpretation='Actualcenterviolation vs conservativeboundunresolved/budget distinguished; engineeringassumption not accuracyproof',
    stopping_point='Exactterminal/allfrozen/idleUNDERlease/separateverify/FINALrelease; no automaticbody/encoder/fieldchange')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(archive/'proposal.json')
