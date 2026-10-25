"""Pure presentation-preview archive preparation: no renderer import or evaluation."""
import argparse
import hashlib
import json
from pathlib import Path
from importlib.metadata import distribution,version
import shutil
ROOT=Path(__file__).resolve().parents[2]
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(p):return p.resolve().relative_to(ROOT).as_posix()
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['preview','encode'],default='preview');parser.add_argument('--preview-review',type=Path);args=parser.parse_args()
base=ROOT/'artifacts/reference/stronger-buoy-20261024'
out=base/('presentation-preview-attempt3' if args.stage=='preview' else 'presentation-encode-attempt2');archive=base/(out.name+'-source')
assert not out.exists() and not archive.exists()
controls_source=ROOT/'reference/wecsim/presentation-preview-controls.json'
controls=read(controls_source);assert controls['status']=='source_only'
accepted=read(ROOT/controls['accepted_controls']);review=read(ROOT/controls['accepted_preview_review'])
assert review['decision']=='preview_pass' and review['reviewed_by']
paths=[ROOT/p for p in accepted['input_sha256']]
assert all(sha(ROOT/p)==v for p,v in accepted['input_sha256'].items())
paths += [ROOT/controls['accepted_controls'],ROOT/controls['accepted_preview_review'],ROOT/review['preflight_json'],ROOT/review['surface_cache_npz']]
assert sha(ROOT/review['preflight_json'])==review['preflight_sha256']
assert sha(ROOT/review['surface_cache_npz'])==review['surface_cache_sha256']
paths += [ROOT/'reference/wecsim/render_reference_presentation.py',ROOT/'reference/wecsim/prepare_presentation_preview.py',ROOT/'reference/wecsim/presentation-preview-controls.json',ROOT/'reference/capytaine/geometry_execution_guard.py']
controls['input_sha256']={key(p):sha(p) for p in dict.fromkeys(paths)}
controls['status']='ready_for_coordinator_presentation_preview'
archive.mkdir(parents=True)
if args.stage=='preview':
 ready=archive/'presentation-ready-controls.json';write(ready,controls)
else:
 ready=base/'presentation-preview-attempt3-source/presentation-ready-controls.json'
 assert read(ready)==controls
 assert args.preview_review
 approved=read(args.preview_review);assert approved['decision']=='preview_pass' and approved['reviewed_by']
 report=ROOT/approved['preflight_json'];assert sha(report)==approved['preflight_sha256']
 controls['input_sha256']=dict(controls['input_sha256'])
frozen=dict(controls['input_sha256']);frozen[key(ready)]=sha(ready)
if args.stage=='encode':
 frozen[key(args.preview_review)]=sha(args.preview_review);frozen[key(report)]=sha(report)
 ffmpeg=Path(shutil.which('ffmpeg'));frozen[str(ffmpeg.resolve())]=sha(ffmpeg)
for name,digest in frozen.items():
 p=ROOT/name
 try:relative=p.resolve().relative_to(ROOT)
 except ValueError:continue
 dest=archive/'frozen'/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest);assert sha(dest)==digest
packages=['numpy','scipy','matplotlib','threadpoolctl','pillow','contourpy','kiwisolver','cycler','fonttools','packaging','pyparsing','python-dateutil','six']
deps={}
for name in packages:
 package=distribution(name)
 deps[name]=dict(version=version(name),file_sha256={str(p):sha(package.locate_file(p)) for p in package.files if str(p).endswith(('.py','.pyd','.dll','.so','.ttf','.afm','.mplstyle'))})
command=['reference/wecsim/render_reference_presentation.py','--config',key(ready),'--output',key(out),'--mode',args.stage]
if args.stage=='encode':command+=['--preview-review',key(args.preview_review)]
proposal=dict(owner='/root/wecsim_reference',status='ready_for_coordinator_execution_allocation',assignment='ONE selected24 presentation '+args.stage+'; no physics or field regeneration',output=key(out),producer_argv=command,input_sha256=frozen,dependency_identity=deps,verification_receipt_name='verification.json',cpu_threads=1,matlab_engines=0,maximum_elapsed_seconds=120 if args.stage=='preview' else 1500,expected_elapsed_seconds=[5,30] if args.stage=='preview' else [120,1400],stopping_point='Preserve exact terminal/hash/idle-underlease/separateverify/finalrelease. ROOT literal view/contactselection before anyfullvideo; no auto appearance rerun.',display='Selected24/meanbow-110/10m/equalmetres, originalhull depthbuffer+vectorizedsamecachedwatertriangle rayhit, exact2px contact, no axes/footer/fieldregen/amplification',physical_limits='Incident-surface contact visualization ONLY, no flood/wet-body/newphysics claim. Closed display convexvolume mustmatchactualphysicalsurface or failclosed.')
write(archive/'proposal.json',proposal);print(key(archive/'proposal.json'));print(sha(archive/'proposal.json'))
