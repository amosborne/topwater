"""Pure immutable freeze for one cached public-CLI lighting rerender."""
import hashlib,json,shutil
from pathlib import Path
from importlib.metadata import distribution,version
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def key(p):return p.resolve().relative_to(ROOT).as_posix()
base=ROOT/'artifacts/reference/stronger-buoy-20261024';archive=base/'pyvista-clip-attempt3-source'
assert not archive.exists();archive.mkdir()
control=ROOT/'reference/wecsim/pyvista-cached-controls.json';c=read(control)
assert not (ROOT/c['output']).exists() and not (ROOT/c['qa_directory']).exists()
paths=[ROOT/p for p in c['input_sha256']]
assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
paths += [ROOT/'reference/wecsim'/p for p in ['run_pyvista_cached.py','pyvista-cached-controls.json','pyvista-cached-settings.json','prepare_pyvista_cached_manifest.py']]
paths += list((ROOT/'src/topwater').rglob('*.py'))+[ROOT/'pyproject.toml',ROOT/'reference/capytaine/geometry_execution_guard.py']
frozen={}
for path in dict.fromkeys(paths):
 frozen[key(path)]=sha(path);dest=archive/'frozen'/key(path);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
for name in ['uv','ffmpeg','ffprobe']:
 p=Path(shutil.which(name)).resolve();frozen[str(p)]=sha(p)
previous=read(base/'pyvista-clip-attempt2-source/proposal.json');deps=previous['dependency_identity']
for name,d in deps.items():
 package=distribution(name);assert version(name)==d['version']
 assert all(sha(package.locate_file(p))==v for p,v in d['file_sha256'].items()),name
proposal=dict(owner='/root/visualization',status='ready_for_coordinator_execution_allocation',assignment='ONE exactcached63frame publicPyVista rerender shadows/waterlineOFF continuouswaterpalette; no waveprep/tests/install/physics',output=c['output'],producer_argv=['reference/wecsim/run_pyvista_cached.py','--config',key(control)],input_sha256=frozen,dependency_identity=deps,verification_receipt_name='verification.json',maximum_elapsed_seconds=90,expected_elapsed_seconds=[20,60],cpu_threads=1,matlab_engines=0,python_executable='artifacts/dependencies/pyvista-venv/Scripts/python.exe',scene_sha256=c['input_sha256'][c['scene']],qa_directory=c['qa_directory'],checks='StandardCLI exactMP4 path, actualbackend/timing, full63frame ffprobe decode800square12.5fps5.04s, decodedpeakposter; unchangedNPZ before/after',stopping_point='Exactterminal/children/frozeninputs, idleUNDERlease/SEPARATEreadback/FINALrelease; ROOTliteralquality, no autoretry/newfield/renderfeatures')
p=archive/'proposal.json';p.write_text(json.dumps(proposal,indent=2)+'\n');print(key(p));print(sha(p))
