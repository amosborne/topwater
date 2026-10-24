"""Freeze one future native mildbuoy case after accepted immutable field assessment."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


base = Path('artifacts/reference/calm-buoy-20261024')
packet = base/'train-attempt2'
assessment = json.loads((base/'assessment-attempt1/diagnostic.json').read_text())
assert assessment['complete'] and assessment['accepted'] and assessment['field_status'] == 'proved_calm'
for path,digest in assessment['immutable_input_sha256'].items():
    assert sha(path) == digest,path
out = base/'case-attempt1'
archive = base/'case-attempt1-source-ready2'
assert not out.exists()
archive.mkdir(parents=True,exist_ok=False)
names = ['run_calm_buoy_case.m','wecSimInputFile_calm_buoy.m','prepare_calm_buoy_case_inputs.m',
    'apply_jan90_force_table.m','log_mild_packet_cic.m','verify_reference_extension_cic.m',
    'verify_seeded_cic_output.m','build_reference_model.m','run_matlab.py','prepare_calm_buoy_case_manifest.py']
paths = [Path('reference/wecsim')/p for p in names]
paths += list(Path('artifacts/dependencies/topwater-wec-sim/source').rglob('*'))
paths += [Path('reference/capytaine/geometry_execution_guard.py'),Path('reference/capytaine/canonical-boat.json'),
          Path('docs/experts/wecsim-fork.json')]
old = Path('artifacts/reference/minimum-20261024')
paths += [old/'native-gl2/native.mat',old/'native-gl2/boat.stl',old/'consumer/canonical.h5']
paths += list(packet.iterdir())
paths += list((base/'assessment-attempt1').iterdir())
paths += [base/'export-attempt1/case.json',base/'export-attempt1/components.npz',
          base/'excitation-attempt1/excitation.npz',base/'excitation-attempt1/manifest.json']
frozen = {}
for path in dict.fromkeys(paths):
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
    assignment='ONE60s mildbuoy physicalincoming train2 nativecase afterexplicitREF002fieldassessmentPASS',
    output=out.as_posix(),input_sha256=frozen,dependency_identity=deps,
    producer_argv=['reference/wecsim/run_matlab.py','reference/wecsim/run_calm_buoy_case.m',
                   '--log',(base/'case-attempt1.log').as_posix()],
    verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
    expected_elapsed_seconds=[30,180],cpu_threads=1,matlab_engines=1,
    operator='ORIGINALminimumH5/nativeGL2M/K/Ainf/STL/CIC10/dt.02/nativeAr unchanged',
    frame='NativeWECmodelworld=meanboat; sourceHtheta/force/M unchanged. GlobalstationR201.308/t embedding metadata only, apply exactlyonce toCG/geometry/vectors/video.',
    initial_state='Knownneutralq0/allv0/CICmemory0; ramp0/viscosity0/all6free; no history adapter/prerun',
    excitation='Unchangedacceptedtrain2 spectralHforce onceviareviewedcase-onlyhook; scalarzeroelevation nativeimport compatibilityONLY',
    checks='3001rawfinite60channels/full6Fe/actualRT-PCHIP-CIC/NE/Ar/H5fullM-K-Ainf-CG/inertia/headinglineage',
    limits='REF002engineeringcalmpolicy explicit, old.1Nrejects preserved; numericalaccuracyunqualified',
    stopping_point='Exactterminal/idleUNDERlease/separateverify/FINALrelease; no automaticrerun/encoder/fieldchange')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(archive/'proposal.json')
