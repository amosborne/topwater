"""Pure freeze of one PyVista hardware/poster/63frame+synthetic smoke job."""
import hashlib
import json
from pathlib import Path
from importlib.metadata import distributions
import shutil
ROOT=Path(__file__).resolve().parents[2]
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def key(p):return p.resolve().relative_to(ROOT).as_posix()
base=ROOT/'artifacts/reference/stronger-buoy-20261024';archive=base/'pyvista-clip-attempt2-source';out=base/'pyvista-clip-attempt2.mp4'
assert not archive.exists() and not out.exists()
control=ROOT/'reference/wecsim/pyvista-reference-controls.json';config=read(control)
paths=[ROOT/name for name in config['input_sha256']]
assert all(sha(ROOT/name)==digest for name,digest in config['input_sha256'].items())
paths += [ROOT/'reference/wecsim'/name for name in ['prepare_pyvista_scene.py','run_pyvista_reference.py','pyvista-reference-controls.json','pyvista-visualization-requirements.txt','prepare_pyvista_reference_manifest.py','calm_buoy_train.py']]
paths += [ROOT/'src/topwater/__init__.py',ROOT/'src/topwater/visualization.py',ROOT/'src/topwater/cli.py',ROOT/'src/topwater/__main__.py',ROOT/'tests/test_visualization.py',ROOT/'tests/test_package.py',ROOT/'pyproject.toml',ROOT/'reference/capytaine/geometry_execution_guard.py']
paths += list((ROOT/'src/topwater').rglob('*.py'))
archive.mkdir(parents=True)
frozen={}
for path in dict.fromkeys(paths):
 assert path.is_file(),path
 frozen[key(path)]=sha(path);dest=archive/'frozen'/key(path);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,dest)
deps={}
for package in distributions():
 name=package.metadata['Name'];deps[name]=dict(version=package.version,file_sha256={str(p):sha(package.locate_file(p)) for p in package.files or [] if str(p).endswith(('.py','.pyd','.dll','.so','.ttf','.afm','.mplstyle','.exe','.pth'))})
uv=Path(shutil.which('uv')).resolve();frozen[str(uv)]=sha(uv)
ffmpeg=Path(shutil.which('ffmpeg')).resolve();frozen[str(ffmpeg)]=sha(ffmpeg)
proposal=dict(owner='/root/visualization',status='ready_for_coordinator_execution_allocation',assignment='ONE PyVista/VTK hardware+poster+63frame50-54.96s squareclip andsyntheticconvexmesh smoke; no simulation/customrenderer',output=key(out),producer_argv=['reference/wecsim/run_pyvista_reference.py','--config',key(control)],input_sha256=frozen,dependency_identity=deps,verification_receipt_name='verification.json',cpu_threads=1,matlab_engines=0,maximum_elapsed_seconds=300,expected_elapsed_seconds=[30,120],python_executable='artifacts/dependencies/pyvista-venv/Scripts/python.exe',prepared_scene=config['prepared_scene'],prepared_scene_policy='Same accepted physicalcoefficients sampledwiderDISPLAYONLY for63times, centralacceptedcache bitidentical; reuse once forfuturestylechanges',rendering='Standard generic topwater.visualization PyVista/VTK settings/convexrequirement; backendGPUunproven untilactualcontextreport',checks='Selectedgenericunit tests; mainposter first, all63squareframes/cornerwatercoverage/actualcontextbackend, secondsyntheticboxposter; exactnativecache/inputlineage',stopping_point='Exactterminal+ffmpegchildproof, input/dependencyhashes, idleUNDERlease/SEPARATEreadback/FINALrelease. No automaticfull60s render/physics/installchange.')
(archive/'proposal.json').write_text(json.dumps(proposal,indent=2)+'\n');print(key(archive/'proposal.json'));print(sha(archive/'proposal.json'))
