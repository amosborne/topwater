"""Freeze telemetry augmentation + shared default Matplotlib CLI, no execution."""
import hashlib,json,shutil
from pathlib import Path
from importlib.metadata import distributions
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def key(p):return p.resolve().relative_to(ROOT).as_posix()
base=ROOT/'artifacts/reference/stronger-buoy-20261024';archive=base/'matplotlib-clip-attempt1-source'
assert not archive.exists();cpath=ROOT/'reference/wecsim/reference-2d-controls.json';c=read(cpath)
assert not (ROOT/c['output']).exists() and not (ROOT/c['augmented_scene']).exists() and not (ROOT/c['qa_directory']).exists()
assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
paths=[ROOT/p for p in c['input_sha256']]
paths += [ROOT/'reference/wecsim'/p for p in ['augment_scene_telemetry.py','run_reference_2d.py','reference-2d-controls.json','prepare_reference_2d_manifest.py']]
paths += list((ROOT/'src/topwater').rglob('*.py'))+[ROOT/'pyproject.toml',ROOT/'tests/test_visualization.py',ROOT/'tests/test_package.py',ROOT/'reference/capytaine/geometry_execution_guard.py']
archive.mkdir();frozen={}
for path in dict.fromkeys(paths):
 frozen[key(path)]=sha(path);destination=archive/'frozen'/key(path);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,destination)
for name in ['uv','ffmpeg','ffprobe']:
 path=Path(shutil.which(name)).resolve();frozen[str(path)]=sha(path)
deps={}
for package in distributions():
 deps[package.metadata['Name']]=dict(version=package.version,file_sha256={str(p):sha(package.locate_file(p)) for p in package.files or [] if str(p).endswith(('.py','.pyd','.dll','.so','.exe','.pth','.ttf','.afm','.mplstyle'))})
proposal=dict(owner='/root/visualization',status='ready_for_coordinator_execution_allocation',assignment='ONE exactcachedgeometry telemetryaugmentation +cleanPUBLICdefaultMatplotlibtopdown63frameproof; no fieldprep/physics',output=c['output'],producer_argv=['reference/wecsim/run_reference_2d.py','--config',key(cpath)],input_sha256=frozen,dependency_identity=deps,verification_receipt_name='verification.json',maximum_elapsed_seconds=90,expected_elapsed_seconds=[10,60],cpu_threads=1,matlab_engines=0,python_executable='artifacts/dependencies/pyvista-venv/Scripts/python.exe',geometry_policy='Originalscene1079ae90 kept; newNPZ arrays time/hull/water/connectivity BITEXACT, only optionaltelemetryfromexactacceptednativeq; no rezero/resample/headings/scale/waves',public_cli='uv run --no-sync topwater render INPUT.npz OUTPUT.mp4; defaultMatplotlib2D, fixed800square,fullNPZ,MP4onlypublicoutput',checks='Perarrayidentity +exactnative-time telemetry, full63decode800x800/12.5fps/5.04s andexplicitreferenceQApeakposter',stopping_point='Exactterminal/children/frozenhashes, idleUNDERlease/SEPARATEreadback/FINALrelease. ROOTliteralreview; no autocorrection/rerun/fieldchanges')
p=archive/'proposal.json';p.write_text(json.dumps(proposal,indent=2)+'\n');print(key(p));print(sha(p))
