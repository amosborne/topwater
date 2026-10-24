"""Archive one prospective incident/body job; no numerical work or engine connection."""
import argparse
import hashlib
import json
import shutil
from importlib.metadata import distribution, version
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


parser = argparse.ArgumentParser()
parser.add_argument('--stage', choices=['incident','body'], required=True)
args = parser.parse_args()
base = Path('artifacts/reference/mild-packet-20261024')
name = 'incident-attempt1' if args.stage == 'incident' else 'case-attempt1'
out = base/name
archive = base/(name+'-source-ready2')
assert not out.exists() and not archive.exists()
sources = ['mild-packet-controls.json','prepare_mild_packet_manifest.py']
if args.stage == 'incident':
    sources += ['prepare_mild_packet.py','incident_packet.py']
else:
    sources += ['run_mild_packet_case.m','prepare_mild_packet_inputs.m','wecSimInputFile_mild_packet.m',
        'verify_mild_packet_cic.m','verify_seeded_cic_output.m','log_mild_packet_cic.m',
        'build_generated_model.m','apply_jan90_force_table.m','run_matlab.py']
paths = [Path('reference/wecsim')/s for s in sources]
paths += [Path('reference/capytaine/canonical-boat.json'),Path('reference/capytaine/geometry_execution_guard.py')]
paths += [Path('artifacts/reference/minimum-20261024/native-gl2/canonical.nc'),
    Path('artifacts/reference/minimum-20261024/consumer/bemio.mat')]
if args.stage == 'body':
    paths += list(Path('artifacts/dependencies/topwater-wec-sim/source').rglob('*'))
    paths += [Path('artifacts/reference/generated-bank-native-20261024/attempt1/native.mat'),
        Path('artifacts/reference/generated-bank-native-20261024/attempt1/boat.stl'),
        Path('artifacts/reference/generated-bank-consumer-20261024/bemio-attempt1/canonical.h5')]
    packet = base/'incident-attempt1'
    report = json.loads((packet/'diagnostic.json').read_text())
    assert report['complete'] and report['accepted']
    report_hashes = {key.replace('\\','/'):value for key,value in report['input_sha256'].items()}
    assert report_hashes['reference/wecsim/mild-packet-controls.json'] == sha('reference/wecsim/mild-packet-controls.json')
    paths += list(packet.glob('*'))
archive.mkdir(parents=True)
frozen = {}
for path in dict.fromkeys(paths):
    if not path.is_file():
        continue
    frozen[path.as_posix()] = sha(path)
    dest = archive/path
    dest.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(path,dest)
deps = {}
for name in (['numpy','scipy','xarray','threadpoolctl'] if args.stage == 'incident' else ['matlabengine']):
    package = distribution(name)
    deps[name] = dict(version=version(name), file_sha256={str(p):sha(package.locate_file(p))
        for p in package.files if str(p).endswith(('.py','.pyd','.dll','.so'))})
command = ['reference/wecsim/prepare_mild_packet.py','--output',out.as_posix()]
if args.stage == 'body':
    command = ['reference/wecsim/run_matlab.py','reference/wecsim/run_mild_packet_case.m',
               '--log',(base/'case-attempt1.log').as_posix()]
proposal = dict(owner='/root/wecsim_reference',status='ready_for_coordinator_execution_allocation',
    assignment='ONE mild single incoming packet '+args.stage+' diagnostic; no finalbuoysea or BEM accuracy claim',
    output=out.as_posix(),producer_argv=command,input_sha256=frozen,dependency_identity=deps,
    verification_receipt_name='verification.json',maximum_elapsed_seconds=600,
    expected_elapsed_seconds=[120,480] if args.stage == 'incident' else [30,180],
    cpu_threads=1,matlab_engines=0 if args.stage == 'incident' else 1,
    stopping_point='Separate incident/body allocations; preserve failed preflight; no retiming/rerun/encoder',
    inputs='a=.015m,f=.25Hz,theta45deg,phase0,B=.5rad/s,focus40s; twozero-energy evaluator slots',
    startup='Fullhullplanform calm bounds+six direct/retainedF0 before body; body q/v/CIC0, ramp0, viscosity0',
    operator='Existing17panelH5/CIC40 UNQUALIFIED; REF001 +/-20 excitationIRF and futurepadding through80s')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n')
print(archive/'proposal.json')
