"""Reference adapter runner; all drawing is delegated to topwater.visualization."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VTK_SMP_MAX_THREADS'):
    os.environ[key]='1'
import argparse
from pathlib import Path
import sys
import time
import json
import subprocess
import shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from prepare_pyvista_scene import prepare_scene,read,write,sha


def main():
    from topwater.visualization import SceneSequence
    import pyvista as pv
    parser=argparse.ArgumentParser();parser.add_argument('--config',type=Path,required=True);args=parser.parse_args()
    config=read(args.config);started=time.monotonic()
    qa=ROOT/config['qa_directory'];qa.mkdir(parents=True,exist_ok=False)
    env=dict(os.environ,UV_PROJECT_ENVIRONMENT=str(Path(sys.executable).resolve().parents[1]),VIRTUAL_ENV=str(Path(sys.executable).resolve().parents[1]),PYTHONPATH=str(ROOT/'src'),UV_CACHE_DIR=str(ROOT/'artifacts/tools/uv-cache'))
    test=subprocess.run([sys.executable,'-m','pytest','-q','--basetemp',str(qa/'pytest-temp'),'tests/test_visualization.py','tests/test_package.py'],cwd=ROOT,env=env,check=True)
    uv=shutil.which('uv');assert uv,'Public uv entrypoint unavailable'
    assert all(sha(ROOT/p)==v for p,v in config['input_sha256'].items())
    prepared=ROOT/config['prepared_scene'];path=prepared/'scene.npz'
    if path.exists():
        scene=SceneSequence.load(path)
        assert scene.metadata['input_sha256']==config['input_sha256'],'Prepared scene belongs to another physical input'
        assert scene.metadata['preparation_controls']=={k:config[k] for k in scene.metadata['preparation_controls']},'Prepared scene belongs to another clip/grid'
    else:scene=prepare_scene(config,prepared)
    settings=dict(config['render_settings'])
    settings['azimuth_deg']=scene.metadata['mean_heading_deg']+config['camera_azimuth_from_heading_deg']
    output=ROOT/config['output']
    settings['maximum_seconds']=max(1,config['maximum_job_seconds']-(time.monotonic()-started)-5)
    settings_path=prepared/'render-settings.json';write(settings_path,settings)
    def public_render(scene_path,destination,setting_path,label):
        command=[uv,'run','--no-sync','topwater','render',str(scene_path),str(destination),'--settings',str(setting_path)]
        process=subprocess.Popen(command,cwd=ROOT,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True);print(dict(public_cli=command,cli_pid=process.pid,parent_pid=os.getpid()),flush=True)
        stdout,stderr=process.communicate();(qa/f'{label}-stdout.txt').write_text(stdout);(qa/f'{label}-stderr.txt').write_text(stderr)
        assert process.returncode==0,f'Public renderer CLI failed with exit{process.returncode}: {stderr}'
        result=json.loads(stdout);write(qa/f'{label}-report.json',result);return result
    report=public_render(path,output,settings_path,'main')
    assert time.monotonic()-started<config['maximum_job_seconds'],'Job budget exhausted'
    # A closed synthetic box exercises the SAME public renderer with other topology.
    dims=np.asarray(config['smoke_dimensions_m']);xy=(scene.hull_points[0,:,:2].min(axis=0)+scene.hull_points[0,:,:2].max(axis=0))/2
    center=np.r_[xy,config['smoke_center_z_m']]
    box=pv.Box(bounds=np.c_[center-dims/2,center+dims/2].ravel()).triangulate()
    smoke=SceneSequence(time_s=scene.time_s[:2],hull_points=np.repeat(np.asarray(box.points)[None],2,axis=0),hull_faces=np.asarray(box.faces).reshape(-1,4)[:,1:],water_points=scene.water_points[:2],water_faces=scene.water_faces,metadata=dict(units='m',coordinate_system='right-handed-z-up',source='Synthetic convex topology/framing smoke only; not physical response'))
    smoke_settings=dict(settings,maximum_seconds=max(1,config['maximum_job_seconds']-(time.monotonic()-started)))
    smoke.save(prepared/'synthetic-box.npz');write(prepared/'smoke-settings.json',smoke_settings)
    smoke_report=public_render(prepared/'synthetic-box.npz',ROOT/config['synthetic_smoke_output'],prepared/'smoke-settings.json','synthetic')
    ffmpeg=shutil.which('ffmpeg');assert ffmpeg
    for movie,index,name in [(output,18,'poster.png'),(ROOT/config['synthetic_smoke_output'],0,'synthetic-poster.png')]:
        subprocess.run([ffmpeg,'-v','error','-i',str(movie),'-vf',f'select=eq(n\\,{index})','-frames:v','1',str(qa/name)],check=True)
    assert all(sha(ROOT/p)==v for p,v in config['input_sha256'].items())
    write(qa/'reference-job.json',dict(complete=True,elapsed_s=time.monotonic()-started,scene_sha256=sha(path),input_sha256=config['input_sha256'],render_report=report,synthetic_smoke_report=smoke_report,limits='Same fixed physical trajectory/train; visual convex closure, not flooding/scatteredfluid/nonlinear accuracy'))


if __name__=='__main__':main()
