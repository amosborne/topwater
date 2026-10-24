"""Freeze prospective ORIGINAL REF001 H5/CIC10 body duration extension; never launch."""
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


base = Path('artifacts/reference/reference-extension-20261024')
packet = base/'incident-attempt1'
report = json.loads((packet/'diagnostic.json').read_text())
assert report['accepted'] and report['unchanged_physical_field']
assert report['first501_retained_force_bit_identical'] and report['old_wave_prefix_bytes_identical']
report_hashes = {key.replace('\\','/'):value for key,value in report['input_sha256'].items()}
controls = Path('reference/wecsim/reference-extension-controls.json')
assert sha(controls) == report_hashes[controls.as_posix()]
out = base/'case-attempt1'
archive = base/'case-attempt1-source'
assert not out.exists()
archive.mkdir(parents=True,exist_ok=False)
names = ['run_reference_extension_case.m','wecSimInputFile_reference_extension.m',
    'prepare_reference_extension_inputs.m','inspect_reference_extension_dispatch.m',
    'verify_reference_extension_cic.m','log_mild_packet_cic.m','verify_seeded_cic_output.m',
    'build_reference_model.m','run_matlab.py','prepare_reference_extension_body_manifest.py',
    'reference-extension-controls.json']
paths = [Path('reference/wecsim')/p for p in names]
paths += list(Path('artifacts/dependencies/topwater-wec-sim/source').rglob('*'))
paths += [Path('reference/capytaine/geometry_execution_guard.py'),Path('reference/capytaine/canonical-boat.json'),
          Path('docs/experts/wecsim-fork.json')]
old = Path('artifacts/reference/minimum-20261024')
paths += [old/'native-gl2/native.mat',old/'native-gl2/boat.stl',old/'consumer/canonical.h5',
          old/'case-attempt1/motion.mat',old/'expected-retained-force.mat']
paths += list(packet.iterdir())
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
    assignment='ONE ORIGINAL REF001 unchangedpacket/H5/CIC10 native60s duration extension',
    output=out.as_posix(),input_sha256=frozen,dependency_identity=deps,
    producer_argv=['reference/wecsim/run_matlab.py','reference/wecsim/run_reference_extension_case.m',
                   '--log',(base/'case-attempt1.log').as_posix()],
    verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
    expected_elapsed_seconds=[30,180],cpu_threads=1,matlab_engines=1,
    operator='ORIGINAL minimum-20261024/consumer/canonical.h5/native-gl2 M/K/Ainf/STL/CIC10/dt.02',
    actual_excitation='Native3import histories/preprocessing/onepositiveTOTALbranch; NOforcehookoverride',
    initial_state='Knownneutralq0/allv0/CICmemory0; ramp0/viscosity0/all6free; no history adapter/prerun',
    checks='3001rawfinite60channels/full6prepared+appliedforce/actualRT-PCHIP-CIC/NE/Ar/first501q-v repeatability REPORT',
    limits='Existing17panel coefficients/temporalaccuracy UNQUALIFIED; boundedmotion not passivityproof',
    stopping_point='Exactterminal/idleunderlease/separateverify/finalrelease; no automaticrerun/encoder/newBEM')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(archive/'proposal.json')
