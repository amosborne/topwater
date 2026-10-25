"""Freeze one cached shared-CLI rerender with external/process phase timing."""
import hashlib,json,shutil
from pathlib import Path
from importlib.metadata import distributions
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def key(p):return p.resolve().relative_to(ROOT).as_posix()
base=ROOT/'artifacts/reference/stronger-buoy-20261024';cpath=ROOT/'reference/wecsim/render-cached-controls.json';c=read(cpath);archive=base/(Path(c['output']).stem+'-source')
assert not archive.exists() and not (ROOT/c['output']).exists() and not (ROOT/c['qa_directory']).exists()
assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
paths=[ROOT/p for p in c['input_sha256']]+[ROOT/'reference/wecsim'/p for p in ['run_render_cached.py','render-cached-controls.json','prepare_render_cached_manifest.py']]
paths += list((ROOT/'src/topwater').rglob('*.py'))+[ROOT/'pyproject.toml',ROOT/'tests/test_visualization.py',ROOT/'reference/capytaine/geometry_execution_guard.py']
archive.mkdir();frozen={}
for path in dict.fromkeys(paths):
 frozen[key(path)]=sha(path);dest=archive/'frozen'/key(path);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
for name in ['uv','ffmpeg','ffprobe']:
 path=Path(shutil.which(name)).resolve();frozen[str(path)]=sha(path)
deps={}
for package in distributions():
 deps[package.metadata['Name']]=dict(version=package.version,file_sha256={str(p):sha(package.locate_file(p)) for p in package.files or [] if str(p).endswith(('.py','.pyd','.dll','.so','.exe','.pth','.ttf','.afm','.mplstyle'))})
proposal=dict(owner='/root/visualization',status='ready_for_coordinator_execution_allocation',assignment=c.get('purpose','ONE exactcachedtelemetryNPZ defaultMatplotlib style/timing63frameproof; no augmentation/waveprep/physics/benchmark'),output=c['output'],producer_argv=['reference/wecsim/run_render_cached.py','--config',key(cpath)],input_sha256=frozen,dependency_identity=deps,verification_receipt_name='verification.json',maximum_elapsed_seconds=90,expected_elapsed_seconds=[10,60],cpu_threads=1,matlab_engines=0,python_executable='artifacts/dependencies/pyvista-venv/Scripts/python.exe',scene_sha256=c['input_sha256'][c['scene']],timing='External Popen-to-communicate CLIwall includesuv/interpreterload; preserve returned phase JSON and explicit QA/runnerwall',checks='Scenehashbefore/after; full63decode800x800/12.5fps/5.04s and peakPNG explicitreferenceQA; focusedtests enabled onlyifcontrollerrequests' if c.get('run_focused_tests',True) else 'No testrepeat forfont-only; scenehashbefore/after,full63decode800x800/12.5fps/5.04s,peakPNGexplicitreferenceQA',stopping_point='Exactterminal/children/frozenhashes idleUNDERlease/SEPARATEreadback/FINALrelease; ROOTliteral/timinganalysis no extra benchmarks or physics')
p=archive/'proposal.json';p.write_text(json.dumps(proposal,indent=2)+'\n');print(key(p));print(sha(p))
