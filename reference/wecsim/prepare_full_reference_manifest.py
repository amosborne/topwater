"""Freeze full751 accepted-cache/minimal-border producer and standard 2D CLI."""
import hashlib,json,shutil
from pathlib import Path
from importlib.metadata import distributions
ROOT=Path(__file__).resolve().parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def key(p):return p.resolve().relative_to(ROOT).as_posix()
base=ROOT/'artifacts/reference/stronger-buoy-20261024';archive=base/'matplotlib-full60-attempt1-source';control=ROOT/'reference/wecsim/full-reference-controls.json';c=read(control)
assert not archive.exists() and not (ROOT/c['output']).exists() and not (ROOT/c['scene_output']).exists() and not (ROOT/c['qa_directory']).exists()
assert all(sha(ROOT/p)==v for p,v in c['input_sha256'].items())
paths=[ROOT/p for p in c['input_sha256']]+[ROOT/'reference/wecsim'/p for p in ['prepare_full_reference_scene.py','run_full_reference.py','full-reference-controls.json','prepare_full_reference_manifest.py','calm_buoy_train.py']]
paths+=list((ROOT/'src/topwater').rglob('*.py'))+[ROOT/'pyproject.toml',ROOT/'reference/capytaine/geometry_execution_guard.py']
archive.mkdir();frozen={}
for p in dict.fromkeys(paths):
 frozen[key(p)]=sha(p);destination=archive/'frozen'/key(p);destination.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,destination)
for name in ['uv','ffmpeg','ffprobe']:
 p=Path(shutil.which(name)).resolve();frozen[str(p)]=sha(p)
deps={}
for package in distributions():
 deps[package.metadata['Name']]=dict(version=package.version,file_sha256={str(p):sha(package.locate_file(p)) for p in package.files or [] if str(p).endswith(('.py','.pyd','.dll','.so','.exe','.pth','.ttf','.afm','.mplstyle'))})
proposal=dict(owner='/root/visualization',status='ready_for_coordinator_execution_allocation',assignment='ONE full60s751frameacceptedcache/minimal1550borderpoints producer+unchanged2DpublicCLI; no sea/physics changes',output=c['output'],producer_argv=['reference/wecsim/run_full_reference.py','--config',key(control)],input_sha256=frozen,dependency_identity=deps,verification_receipt_name='verification.json',maximum_elapsed_seconds=600,expected_elapsed_seconds=[120,240],cpu_threads=1,matlab_engines=0,python_executable='artifacts/dependencies/pyvista-venv/Scripts/python.exe',data_policy='All751acceptedhull/time/nativeHUD; central65x65waterXYZBITEXACT; SAMEfrozentrain evaluatedONLYnewminimumoriginalspacingborder nodes75x77; no96mgrid/generator/physics',checks='Exactframing1.45/allhulls/fixedviewcornercoverage; scene reloadidentity/nativeqindices0:4:3000; actualcleanPUBLICCLI751/fixed800/12.5fps/60.08s/full751decode+0/30/51.44/60PNGs; preparation/CLIwall+phaseprofile',stopping_point='Exactterminal/children/frozenidentity, idleUNDERlease/SEPARATEreadback/FINALrelease. ROOTliteral finalreview; no automaticrerun/fieldtuning/rendererchanges')
p=archive/'proposal.json';p.write_text(json.dumps(proposal,indent=2)+'\n');print(key(p));print(sha(p))
